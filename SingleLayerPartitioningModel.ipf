#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


// Simple numerical model solving the transport of a chemical species with a given C* along a PFA Teflon
// tube under laminar flow, as published in Pagonis et al., AMT, 2017.
// Link to AMTD version: https://www.atmos-meas-tech-discuss.net/amt-2017-279/
// See paper for ranges of flow rates, diameters and C* tested experimentally
// May not work well outside those bounds
// Questions? <demetrios.pagonis@colorado.edu> and <jose.jimenez@colorado.edu>

//Launch panel: TubingModelPanel()

// Master function
	// Enter the parameters for the model run
		//c* (saturation concentration) at the temperature of interest
		//total time
		//Analyte diffusion coefficient in air at the pressure of interest
		//tube diameter
		//flow rate
		//tube length
		//bin length
		//Chamber SAV
		//Chamber Cw
		//step size

//Sample initialization function
	//Set up the time series for the gas being sampled

//v2.1 -	fixes a bug in the convolution of the instrument response function that affected cases
//			where concentration values exceeded 1
//v2.2 - 	fixes a bug where the user's Cw input for 'other' tubing was overwritten
//			conductive silicone tubing has been added as a tubing material


Function ModelTubing(CStar, Cw_ugm3, TotalTime_s, Q_Lmin, Dgas_cm2s, Dtube_cm, Ltube_cm, SampleFlag)

	//-------------------------------------------------------------------------------------------------------------------
	//user-defined parameters

	Variable CStar			// compound c*, ug/m3
	Variable Cw_ugm3			// tubing Cw, ug/m3, for 0.47 cm ID tubing
	Variable TotalTime_s	//Total time, seconds - must be an integer number of seconds
	Variable Q_Lmin
	Variable Dgas_cm2s  	//Gas-phase diffusion coeff, cm^2/s, 0.067 is used in the paper
	Variable Dtube_cm 		//Tube diameter, cm 0.47= 3/16" ID (1/4" OD), 0.15875 cm = 1/16" ID (1/8" OD)
	Variable Ltube_cm		//Tube length, cm
	Variable SampleFlag		//Flag for time sample time series


	Variable Lbin_cm = 2	//bin length, cm
	Variable SAV_quarterOD = 8.51064 //Surface area to volume ratio, 0.47 cm ID tubing
	Variable StepSize = 	0.001	//Euler method time step, sec, recommend to use 0.001 - must divide evenly into 1 second
	Variable Q_cm3s = Q_Lmin * 1000 / 60 // Volumetric flow rate, cm^3/s

	//-------------------------------------------------------------------------------------------------------------------
	//Rate constant & equilibrium constant calculations


	Variable NPTS = TotalTime_s	//capture model output at 1-second intervals

	If( TotalTime_s < 60)	// capture model output at a higher time resolution when the model run is short
		NPTS *= 25
	EndIf

	Variable UpdateInt = round(NPTS/10)

	Variable Vgas_cm3 = 0.25 * pi* Dtube_cm^2 * Ltube_cm		// Volume of gas in the whole length of tube under study
	Variable SAV_tube = (pi*Dtube_cm*Ltube_cm)/Vgas_cm3		// Calculate the surface area to volume ratio for the tube

	Variable Kgw_cham = Cw_ugm3 / Cstar									//[W]/[G]
	Variable Kgw_tube = Kgw_cham * SAV_tube / SAV_quarterOD		//scale molar equilibrium constants by SA / V ratios



	Variable StepInt = (totalTime_s / Stepsize) / Npts		//number of calculation steps to wait between recording outputs

	Variable kdgw = (2 * Dgas_cm2s) / (Dtube_cm/2)^2			//rate constant for gas-wall transport, kgw - Eq. 4 in Pagonis et al.
	Variable kdwg = kdgw / Kgw_Tube								//rate constant for wall-gas transport, kwg - Eq. 5 in Pagonis et al.
	Variable v_cms = Q_cm3s / ( pi * (Dtube_cm/2)^2 )	//Flow speed assuming plug flow
	Variable kflow = v_cms / Lbin_cm								//rate constant for perfectly mixed flow between bins, kf - Eq. 7 in Pagonis et al.


	Variable Nbins = floor(Ltube_cm / Lbin_cm)


	//-------------------------------------------------------------------------------------------------------------------
	// wave creation / initialization

	Make /N=(Npts,Nbins) /O wGas, wWall				//Arrays that hold concentration of sample in the gas & walls over bin number and time
	Make /N=(Npts) /O wTime, wSample, wLast, wInstrument
	Make /N=(Nbins) /O /D wg, ww, wdg, wdw			//waves to capture the concentration & 1st derivative of sample in gas & wall for each bin
	wLast=NaN												//the derivatives are used to update the concentration waves in the euler method below
	wInstrument = NaN
	DoUpdate

	InitializeSampleWave(wSample, CStar, SampleFlag)			//Function populates the sample wave, defines what the analyte concentration is entering the tube

	wg = 0
	ww = 0


	//-------------------------------------------------------------------------------------------------------------------
	//bookkeeping before the model gets run
	Variable CurrentTime
	Variable iStep = 0
	Variable iRow = 0
	Variable iCol
	Variable iBin
	Variable CurrentSample
	Variable BinW, BinG


	//-------------------------------------------------------------------------------------------------------------------
	// run the model using an euler approximation

	wTime = StepSize * p * StepInt

	Do
		//update the time and samlple variables
		CurrentTime += StepSize
		CurrentSample = interp(CurrentTime, wTime, wSample )

		// Calculate d[G]/dt and d[W]/dt for each bin - Eq. 8 & 9 in Pagonis et al. AMT 2017
		wdg[0] = (kflow * CurrentSample) - (kflow * wg[0]) - (kdgw * wg[0]) + (kdwg * ww[0])
		wdw[0] = kdgw * wg[0] - kdwg * ww[0]

		For( iBin = 1; iBin < NBins ; iBin += 1)
			Binw = ww[iBin]
			Bing = wg[iBin]
			// Eq. 8 & 9 in Pagonis et al. AMT 2017
			wdg[iBin] = (kflow * wg[iBin-1]) - (kflow * Bing) - (kdgw * Bing) + (kdwg * Binw)
			wdw[iBin] = kdgw*Bing - kdwg * Binw

		EndFor

		//Euler method - solving the whole wave (at all positions) at the same time
		//Xi = X(i-1) + dX * dt
		wg += wdg * StepSize
		ww += wdw * StepSize

		//save concentrations to output waves at the specified interval (ie run the model at 1 ms and record output at 1 s)
		if( mod( iStep, StepInt) == 0 && iRow < Npts )

			For( iCol = 0 ; iCol < NBins ; iCol += 1)
				wGas[iRow][iCol] = wg[iCol]
				wWall[iRow][iCol] = ww[iCol]
			EndFor

			If(mod(iRow,UpdateInt)==0)
				wLast[0,iRow]=wGas[p][(NBins-1)]
				DoUpdate
			EndIf

			iRow += 1
		endif


		iStep += 1
	While( CurrentTime < TotalTime_s)


	KillWaves /Z wg,ww,wdg,wdw

	//make a wave for the last bin in the tube
	Make /N=(NPts) /O wLast = wGas[p][(NBins-1)]

	//Convolve tubing model with instrument timescales
	NVAR gA1
	NVAR gTau1
	NVAR gA2
	NVAR gTau2
	Make /N=(NPts) /O wInstrument
	DP_ConvolveXYDbleExp(wTime, wLast, gTau1, gA1, gTau2, gA2, wInstrument)


	Return 0
End


//====================================================================================================================


Static Function InitializeSampleWave(w, CStar, SampleFlag)

	Wave w
	Variable CStar
	Variable SampleFlag

	//--------------------------------------------------------------------------

	If( SampleFlag == 0 )
		w=1		// use this to model step-function increases in measured concentration
		w[0]=0
	ElseIf (SampleFlag == 1)

	//--------------------------------------------------------------------------
		w=0
		w[10,69] = 1
	ElseIf(SampleFlag ==2)

	//--------------------------------------------------------------------------------------
	// this section models the 10-second oxidation experiments done in Krechmer et al 2016

		Variable Cw_ugm3
		If( Cstar < 1)
			Cw_ugm3 = 16
		ElseIf( Cstar < 6e3)
			Cw_ugm3 = 16*(Cstar^0.6)
		Else
			Cw_ugm3 = 20*1000
		Endif 		//from Krechmer 16
		Variable Kgw_cham = Cw_ugm3 / Cstar	//[W]/[G]
		Variable Fg = 1 / (1+Kgw_cham) 	//observed K in Z group chambers
		Variable Fw = 1-Fg
		w = Fw*exp( -1 * (p - 9) / 600 ) + Fg
		w[0,9] = (p+1)/10

	ElseIf(SampleFlag == 3)
		Wave wCustom=root:wCustom
		w = wCustom

	EndIf

	Return 0
End

//-------------------------------------------------------


Function DP_ConvolveXY( wx, wy, sW, TimeConst )
	//convolves XY data with a single-exponential defined by TimeConst
	//creates/overwrites a wave sW for output
	// function zero-pads start of the wave, one-pads the end and normalizes using those regions

	//Convolution:
	// (f*g)(t) = integral( f(tau)g(t-tau)dtau) over all time

	Wave wx,wy
	String sW
	Variable TimeConst

	Variable tau
	Variable deltatau = 0.05


	Variable NPTS = numpnts(wy)
	Variable NZeroPad = 20
	Variable NOnePad  =  (TimeConst < 10) ? 60 : 10*ceil(TimeConst)
	Variable NPTS_redim = NPTS + NZeroPad + NOnePad
	Variable iPt
	Variable IntegralSum
	Variable xi
	Variable OnePadVal = wavemax(wy)

	Variable Term1, Term2

	//Redimension wX and wY, zero- and one-pad
	InsertPoints NPTS, NOnePad, wx, wy
	wx[NPTS,(NPTS+NOnePad-1)] = wX[NPTS-1]+p-NPTS+1
	wy[NPTS,(NPTS+NOnePad-1)] = OnePadVal
	InsertPoints 0, NZeroPad, wx, wy
	wx += NZeroPad
	wx[0,NZeroPad-1] = p
	wy[0,NZeroPad-1] = 0

	Make /N=(NPTS_redim) /O $sW = NaN
	Wave wConvolvedData = $sW


	//Do the convolution
	For( iPt = 0 ; iPt < NPTS_redim ; iPt += 1)	//for each point in the output wave (f*g)(t)
		xi = wx[iPt]
		IntegralSum = 0

		For( tau = 0 ; tau <= 15*TimeConst ; tau += deltatau )		//integrate across tau

			If(tau <= xi )
				Term1 = interp((xi-tau), wx, wy)		//g(t-tau) term
				Term2 = exp(-1*(tau)/TimeConst)		//f(tau) term
				IntegralSum += Term1*Term2
			EndIf

		EndFor

		wConvolvedData[iPt] = IntegralSum
	EndFor



	//Normalize to zero- and one-padded regions
	Variable ZeroVal = mean(wConvolvedData,0,(NZeroPad-1))
	Variable OneVal = mean(wConvolvedData,NPTS_redim-11,NPTS_redim-1)
	wConvolvedData -= ZeroVal
	wConvolvedData *= OnePadVal / OneVal


	//Undo redimensioning
	wx -= NZeroPad
	DeletePoints 0, NZeroPad, wx, wy, wConvolvedData
	DeletePoints NPTS, NOnePad, wx, wy, wConvolvedData

	Return 0
End

//-------------------------------------------------------

Function DP_ConvolveXYDbleExp(wX, wY, Tau1, Int1, Tau2, Int2, wDest)

	//Convolve XY data with a double exponential function with parameters
	// tau 1 and 2, relative intensities Int1 and 2
	// create/overwrite destination wave sDest

	Wave wX, wY
	Variable Tau1
	Variable Tau2
	Variable Int1
	Variable Int2
	Wave wDest

	Variable NPTS = numpnts(wX)
	Variable NAdded = (Tau1 < Tau2) ? round(5*Tau2) : round(5*tau1)
	Redimension /N=(NPTS+NAdded) wx, wy
	Variable DelX = wX[1]-wX[0]
	wx[NPTS,(NPTS+NAdded-1)] = wX[p-1] + delX
	wY[NPTS,(NPTS+NAdded-1)] = 1

	Make /N=(numpnts(wX)) /O wConv1, wConv2

	DP_ConvolveXY( wx, wy, "wConv1", Tau1 )
	DP_ConvolveXY( wx, wy, "wConv2", Tau2 )

	wDest = Int1*wConv1 + Int2*wConv2

	Redimension /N=(NPTS) wX,wY

	Killwaves /Z wConv1, wConv2

	Return 0
End

//-------------------------------------------------------




Function ListBoxProc(lba) : ListBoxControl
	STRUCT WMListboxAction &lba

	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave

	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			break
		case 3: // double click
			break
		case 4: // cell selection

			NVAR gTotalTime_s
			NVAR gCStar

			Variable NPTS = gTotalTime_s
			Variable Stepsize = 0.001
			Variable StepInt = (gtotalTime_s / Stepsize) / Npts
			If( gTotalTime_s < 60)	// capture model output at a higher time resolution when the model run is short
				NPTS *= 25
			EndIf
			Make /N=(Npts) /O wTime, wSample
			wTime = StepSize * p * StepInt

			NVAR gSampleFlag
			gSampleFlag = row

			InitializeSampleWave( wSample, gCStar, row )
			SetCw()

		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
			break
	endswitch


	return 0
End


Function SetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval

			Wave wTime
			Wave wSample
			Wave wCustom
			Redimension /N=(round(dval)) wTime, wSample, wCustom
			wTime=p
			SetCw()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function ButtonProc_1(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			BrowseURL/Z "https://docs.google.com/document/d/1CCYiW3cfPjPdT69LoTXVvjDk9uaWVLwRxzDg1YstL94/edit?usp=sharing"

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Static Function SetCw()

	NVAR gCw_ugm3
	NVAR giMat
	NVAR gCStar
	Wave wCw = wTubingCw_ugm3

	switch(giMat)	// numeric switch
		case 0:// PFA
			gCw_ugm3 = (gCStar >=1e4) ? wCw[giMat] : 9e4
			gCw_ugm3 = (gCStar >= 40) ? gCw_ugm3 : (1.16e4*gCStar^0.56)
			gCw_ugm3 = (gCStar >= 5) ?  gCw_ugm3 : 2.7e4
			break
		case 1://cPFA
			gCw_ugm3 = (gCStar >=1e4) ? wCw[giMat] : 6.2e5
			gCw_ugm3 = (gCStar >= 1e2) ? gCw_ugm3 : (1.02e4*gCStar^0.85)
			gCw_ugm3 = (gCStar >= 5) ?  gCw_ugm3 : 4.4e4
			break
		case 2: //FEP
			gCw_ugm3 = (gCStar >=1e4) ? wCw[giMat] : 4.9e5
			gCw_ugm3 = (gCStar >= 1e2) ? gCw_ugm3 : (8.46e4*gCStar^0.38)
			gCw_ugm3 = (gCStar >= 5) ?  gCw_ugm3 : 1.2e5
			break
		case 8: //Other material - let user input
			wCw[giMat] = gCw_ugm3
		default:
			gCw_ugm3 = wCw[giMat]		// all other tubings
	endswitch


	Return 0
End


Function ListBoxProc_1(lba) : ListBoxControl
	STRUCT WMListboxAction &lba

	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave

	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			break
		case 3: // double click
			break
		case 4: // cell selection
			Wave wCw = root:wTubingCw_ugm3
			NVAR giMat
			giMat = row
			SetCw()
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
		case 13: // checkbox clicked (Igor 6.2 or later)
			break
	endswitch

	return 0
End

Function CheckRe()
	NVAR gQ_Lmin
	NVAR gTubeLength_cm
	NVAR gDtube_cm

	Variable Q_cm3s=gQ_Lmin*1000/60
	Variable Re = Q_cm3s*gDtube_cm/(0.1511*pi*(gDtube_cm*0.5)^2)
	If(Re>1500)
		TextBox/W=TubingModelPanel#G0/C/N=Re/F=0/A=MC/X=10.00/Y=10.00 "\\K(65535,0,0)\\JCRe > 1500\r\\Z10Model not validated for turbulent flow"
	Else
		TextBox/W=TubingModelPanel#G0/C/N=Re/F=0/A=MC/X=10.00/Y=10.00 " "
	EndIf
	Return 0
End



Function SetVarProc_1(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
			SetCw()
		case 2: // Enter key
			SetCw()
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			SetCw()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetVarProc_2(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			CheckRe()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ButtonProc(ba) : ButtonControl
	// "Run Model" button

	STRUCT WMButtonAction &ba

	NVAR gCStar, gCw_ugm3, gDgas_cm2s, gDtube_cm, gQ_Lmin, gTotalTime_s, gTubeLength_cm, gSampleFlag

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			SetCw()
			CheckRe()
			ModelTubing(gCStar, gCw_ugm3, gTotalTime_s, gQ_Lmin, gDgas_cm2s, gDtube_cm, gTubeLength_cm, gSampleFlag)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetVarProc_3(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			CheckRe()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function CreateTubingModelGlobals()
	dfref dfr = getdatafolderdfr()
	cd root:
	Variable /G gA1 = 0.5
	Variable /G gA2 = 0.5
	Variable /G gTau1 = 2
	Variable /G gTau2 = 30
	Variable /G gCStar = 1e6
	Variable /G gCw_ugm3 = 8e5
	Variable /G gDgas_cm2s = 0.055
	Variable /G gDtube_cm = 0.47
	Variable /G giMat = 0
	Variable /G gQ_Lmin = 1
	Variable /G gSampleFlag = 1
	Variable /G gTotalTime_s = 300
	Variable /G gTubeLength_cm = 100
	Make /N=(gTotalTime_s) /O wCustom, wInstrument, wLast, wSample, wTime = 0
	Make /N=(gTotalTime_s,floor(gTubeLength_cm/2)) /O wGas, wWall = 0
	Make /N=4 /O wSampleFlags=p
	Make /O /T wSampleNames = {"Step function increase","60 second pulse","Krechmer 16 experiments","Custom (root:wCustom)"}
	Make /N=7 /O wTubingCw_ugm3 = {800000,1.3e+06,2e+06,8e+06,1.2e+07,1.6e+07,2.6e11,NaN}
	Make /O /T wTubingMaterials = {"PFA Teflon","Conductive PFA","FEP Teflon","PEEK","PTFE","Conductive PTFE","Conductive Silicone","Other"}
	wTime=p
	InitializeSampleWave( wSample, gCStar, gSampleFlag )
	SetCw()
	cd dfr
	Return 0
End

///////////////////////////////////////////////////////////////////////////////////
// Button prcedures from Donna Sueper
// The popped graph will have the name 'yy_zz'  // yy and zz can each be strings of length at least 1, can be as long as you like
// There needs to exist a function called yy_populate_zz
// If zz is a string with the chars "graph" in it, the code will know that the popped entity is a graph and not a table
// This entire system can be redone without requiring 3 different strings, xx, yy and zz but am doing it this way to be compatible with my old stuff

Function A_Populate_tubingGraph()  // named yy_populate_zz  use of 'graph' in the name is deliberate
	wave  wSample,wLast, wInstrument,wTime  // for now assumes no data folders and these waves exist in root

	AppendToGraph wSample,wLast vs wTime
	AppendToGraph wInstrument vs wTime
	ModifyGraph margin(left)=72,margin(bottom)=72,margin(top)=21,margin(right)=21,gfSize=16
	ModifyGraph frameStyle=1
	ModifyGraph lSize=3
	ModifyGraph rgb(wSample)=(0,0,0),rgb(wLast)=(1,4,52428),rgb(wInstrument)=(65535,43690,0)
	ModifyGraph mirror=2
	ModifyGraph lblMargin(left)=7,lblMargin(bottom)=8
	Label left "Concentration"
	Label bottom "Time (sec)"
	Legend/C/N=text0/J/A=MC/X=29.72/Y=34.39 "\\s(wSample) Sample\r\\s(wLast) Tubing outlet\r\\s(wInstrument) Instrument"
End


// Generic button control for "pop" buttons.  It assumes that the button control is named appropriately, in xx_yy_zz type of format.
// Will pop tables or graphs from a panel
Function gen_butt_popThis(ctrlName) : ButtonControl
	String ctrlName

	string WindowStr="", WindowPrefixStr="", FuncStr="", NameOfPoppedWindow=""

	NameOfPoppedWindow = removeListItem(0, ctrlName, "_")  // presumes the button control is named like  xx_yy_zz.  This removes the xx_ prefix must be 1 or more characters
	WindowPrefixStr = stringFromList(0,NameOfPoppedWindow, "_") // yy
	WindowStr = removeFromList(WindowPrefixStr, NameOfPoppedWindow, "_")   // zz

	if (strlen(NameOfPoppedWindow)==0 || strlen(WindowPrefixStr)==0 || strlen(WindowStr)==0 )
		Abort "something is wrong with the string length of either NameOfPoppedWindow:"+NameOfPoppedWindow +" or WindowPrefixStr "+WindowPrefixStr+" or WindowStr:"+WindowStr
	endif

	DoWindow/F  $NameOfPoppedWindow  // if graph or table already exists, bring to the front instead of making a new one

	if (V_flag==0)
		if (strsearch(lowerStr(NameOfPoppedWindow), "graph", 0)>=0)
			Display/W=(50,50,550,450)/N=$NameOfPoppedWindow as NameOfPoppedWindow // a blank graph window, set arbitrary size & location
			ShowInfo			//1.51Q
		else
			Edit/W=(50,70,450,350)/N=$NameOfPoppedWindow as NameOfPoppedWindow	// a blank table window, set arbitrary size & location
		endif

		FuncStr = WindowPrefixStr+"_populate_"+WindowStr  // yy_populate_zz
		funcref gen_NoParamProtoFunction populateMe = $FuncStr

		if (Exists(FuncStr))
			populateMe()
		else
			print WindowPrefixStr+"_populate_"+WindowStr
			Abort "Could not find the function named "+WindowPrefixStr+"_populate_"+WindowStr+" to populate the new window. Window named passed is "+NameOfPoppedWindow
		endif
	endif

End


// prototype needed for funcref
Function gen_NoParamProtoFunction()
End


Window TubingModelPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	CreateTubingModelGlobals()
	NewPanel /W=(316,108,1073,515)
	ModifyPanel cbRGB=(61166,61166,61166)
	SetDrawLayer UserBack
	SetDrawEnv linethick= 2,fillfgc= (62708,62708,62708)
	DrawRect 24,63,370,341
	SetDrawEnv fsize= 16,fstyle= 1
	DrawText 24,26,"Partitioning delay model for tubing and instruments v2"
	SetDrawEnv fsize= 10
	DrawText 25,56,"Demetrios Pagonis and Jose-Luis Jimenez\rLiu et al. 2019, \\f02AMT\\f00; Deming et al. 2019, \\f02AMT\\f00; Pagonis et al. 2017, \\f02AMT\\f00"
	SetDrawEnv fsize= 10,fstyle= 1
	DrawText 49,83,"Define tubing conditions:"
	SetDrawEnv fsize= 10,fstyle= 1
	DrawText 215,83,"Define sample:"
	SetDrawEnv fsize= 10,fstyle= 1
	DrawText 45,307,"Define instrument response function:  \\f02R(t) = A\\B1\\M∙e\\S-t/τ1\\M+A\\B2\\M∙e\\S-t/τ2"
	DrawPICT 530,2,0.5,0.5,PICT_2
	SetDrawEnv fname= "Segoe UI",fsize= 10
	DrawText 49,186,"Tubing material:"
	SetDrawEnv fname= "Segoe UI",fsize= 10
	DrawText 215,186,"Time series:"
	SetDrawEnv fsize= 9
	DrawText 202,285,"C\\Bw\\M = f(C*) for PFA,\rcPFA, and FEP"
	SetVariable CStar,pos={261.00,90.00},size={94.00,20.00},bodyWidth=40,proc=SetVarProc_1,title="C* (μg m\\S-3\\M)"
	SetVariable CStar,help={"Saturation vapor concentration 298 K"},fSize=10
	SetVariable CStar,limits={0,inf,0},value= gCStar
	SetVariable Cw,pos={84.00,262.00},size={105.00,23.00},bodyWidth=50,title="C\\Bw\\M (μg m\\S-3\\M)"
	SetVariable Cw,help={"Equivalent absorbing mass\rof the tubing walls"},fSize=10
	SetVariable Cw,limits={-inf,inf,0},value= gCw_ugm3
	SetVariable Dgas,pos={259.00,117.00},size={96.00,23.00},bodyWidth=40,title="D\\Bg\\M (cm\\S2\\M s\\S-1\\M)"
	SetVariable Dgas,help={"Analyte diffusivity in air"},fSize=10
	SetVariable Dgas,limits={-inf,inf,0},value= gDgas_cm2s
	SetVariable dtube,pos={58.00,146.00},size={131.00,15.00},bodyWidth=40,proc=SetVarProc_3,title="Inner diameter (cm)"
	SetVariable dtube,fSize=10,limits={-inf,inf,0},value= gDtube_cm
	SetVariable Q,pos={65.00,90.00},size={124.00,20.00},bodyWidth=40,proc=SetVarProc_2,title="Flow rate (L min\\S-1\\M)"
	SetVariable Q,fSize=10,limits={-inf,inf,0},value= gQ_Lmin
	SetVariable time,pos={267.00,146.00},size={88.00,15.00},bodyWidth=40,proc=SetVarProc,title="Time (sec)"
	SetVariable time,help={"Length of model run"},fSize=10
	SetVariable time,limits={-inf,inf,0},value= gTotalTime_s
	SetVariable length,pos={69.00,120.00},size={120.00,15.00},bodyWidth=40,title="Tube length (cm)"
	SetVariable length,fSize=10,limits={-inf,inf,0},value= gTubeLength_cm
	Button RunModel,pos={91.00,351.00},size={200.00,45.00},proc=ButtonProc,title="Run Model"
	Button RunModel,fSize=10,fColor=(65535,65535,0)
	ListBox states,pos={215.00,190.00},size={140.00,60.00},proc=ListBoxProc,fSize=10
	ListBox states,listWave=root:wSampleNames,mode= 2,selRow= 1,editStyle= 1
	SetVariable dtube1,pos={58.00,313.00},size={54.00,18.00},bodyWidth=40,title="A\\B1"
	SetVariable dtube1,help={"Fraction of signal responding with timescale tau1"}
	SetVariable dtube1,fSize=10,limits={-inf,inf,0},value= gA1
	SetVariable dtube2,pos={126.00,313.00},size={66.00,18.00},bodyWidth=40,title="τ\\B1\\M (s)"
	SetVariable dtube2,fSize=10,limits={-inf,inf,0},value= gTau1
	SetVariable dtube3,pos={207.00,313.00},size={54.00,18.00},bodyWidth=40,title="A\\B2"
	SetVariable dtube3,help={"Fraction of signal responding with timescale tau2"}
	SetVariable dtube3,fSize=10,limits={-inf,inf,0},value= gA2
	SetVariable dtube4,pos={275.00,313.00},size={66.00,18.00},bodyWidth=40,title="τ\\B2\\M (s)"
	SetVariable dtube4,fSize=10,limits={-inf,inf,0},value= gTau2
	ListBox materials,pos={49.00,190.00},size={140.00,60.00},proc=ListBoxProc_1
	ListBox materials,fSize=10,listWave=root:wTubingMaterials,mode= 2,selRow= 0
	ListBox materials,editStyle= 1
	Button readme,pos={446.00,6.00},size={40.00,20.00},proc=ButtonProc_1,title="Help"
	Button readme,fSize=10
	Button blah_A_tubinggraph,pos={385.00,34.00},size={65.00,20.00},proc=gen_butt_popThis,title="Pop graph"
	Button blah_A_tubinggraph,fSize=10
	Execute/Q/Z "SetWindow kwTopWin sizeLimit={33,48,inf,inf}" // sizeLimit requires Igor 7 or later
	Display/W=(386,62,750,394)/HOST=#  wSample,wLast vs wTime
	AppendToGraph wInstrument vs wTime
	ModifyGraph margin(left)=65,margin(bottom)=58,margin(top)=21,margin(right)=21,gfSize=12
	ModifyGraph frameStyle=1
	ModifyGraph lSize=3
	ModifyGraph rgb(wSample)=(0,0,0),rgb(wLast)=(1,4,52428),rgb(wInstrument)=(65535,43690,0)
	ModifyGraph mirror=2
	ModifyGraph lblMargin(left)=7,lblMargin(bottom)=8
	Label left "Concentration"
	Label bottom "Time (sec)"
	Legend/C/N=text0/J/A=MC/X=25.00/Y=36.76 "\\s(wSample) Sample\r\\s(wLast) Tubing outlet\r\\s(wInstrument) Instrument"
	TextBox/C/N=Re/F=0/A=MC/X=10.00/Y=10.00 " "
	RenameWindow #,G0
	SetActiveSubwindow ##
EndMacro
