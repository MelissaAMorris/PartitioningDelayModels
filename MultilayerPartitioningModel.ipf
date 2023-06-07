#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Multilayer tubing model, written by Melissa Morris, adapted from Lucas Algrim, 2020
//
// This model is used to simulate transmission of VOCs through a tube made of an absorbing material
// There are many manual parameters that require entry in the "INPUT PARAMETERS LIST"
// 
// For more details on the model, please see the publication: ___________
//
// Version 1.0 finished on 14-Mar-2023
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Function MultiLayerTubingModel(str)

	string str // name for the run
	
	// ============================================================================
	// ============================================================================
	// START OF INPUT PARAMETERS LIST
	
	// The following parameters pertain to the VOC and absorption medium, and require manual entry.
	// Values are pre-loaded to simulate a C8 2-ketone sampled through polyethylene tubing.
	
	variable/G gDf   		= -6.8				// (cm2/sec)	LOG(diffusion coefficient for VOC diffusing in polymer wall)
	variable/G gDfG			= 0.067			// (cm2/sec)	Diffusion coefficient for VOC diffusing in gas phase while traveling through tube
	variable/G gGamma		= 0.02				// 				Availability parameter (percent of mass actually participating in partitioning)
	variable/G gSR			= 1.2				// 				Surface roughness factor (=8 for paint, =1 for smooth tube)
	variable/G gCstar		= 1e6				// (ug/m3)	Saturation mass concentration, C*, calculated using an assumed molar mass of 250 g/mol (Algrim, 2020)
	variable/G gFin  		= 1.86   			// (LPM) 		Flow rate  
	variable/G gDiam 		= 0.483 			// (cm)		Inner diameter of tube
	variable/G gLength 		= 117 				// (cm)		Length of tube
	variable/G gDensity 	= 1.21				// (g/cm3) Density of the absorbing material
	
	// The following parameters pertain to the 3-dimensional linear approximation this model uses, and require manual adjustment.
	// The model run time is very sensitive to gDeltaY, and gDeltaT
	// The code may not converge if 1 of these 3 parameters is too large or two small compared to the other 2. 
	// I recommend running a grid of values to make sure the model converges on an output time series properly.

	variable/G gDeltaX		= 2					// (cm)		Distance step for flow down the tube
	variable/G gDeltaY 		= 2					// 	(um)		Distance step for diffusion into the wall
	variable/G gDeltaT 		= 0.01				// (sec)		Time step for Euler method
	
	// The following parameters pertain to the timing of the simulated experiment, and require manual entry.
	
	variable/G gOG   		= 4071 			// (sec)		Time when concentration of VOC is changed from gCin to 0. This is useful for simulating a passivation period and then a depassivation period.
	variable/G gTotalSec	= 14000			// (sec)		Total time for model to simulate
	variable/G gCin  		= 1     			// 				Ambiguous conc. of VOC. (kept at 1)
   
   	// END OF INPUT PARAMETERS LIST
   // ============================================================================
	// ============================================================================

   // ============================================================================
	// ============================================================================
	// START OF OUTPUT PARAMETERS LIST
   
   Variable/G vCwTotal					// Final Cw value calculated from the modeled penetration depth
   Variable/G vPenetrationDepth		// Penetration depth of VOC into absorbing material
   Variable/G vThickness				// Total thickness that was modeled
   
   	// END OF OUTPUT PARAMETERS LIST
   // ============================================================================
	// ============================================================================
   
	// ============================================================================
	// ============================================================================
	// START OF LOCAL VARIABLE CREATION
   
	// Making local variables with standard units
   
	variable vDiam 		= gDiam / 100 		// (m)
	variable vDeltaX 	= gDeltaX / 100		// (m)
	variable vDeltaY	= gDeltaY / 1e6		// (m)
	variable vDensity 	= gDensity * 1e6 * (100)^3 // (ug/m3) 
	Variable vCin 		= gCin
   
	// Time variables
	variable starttime= datetime 						// helps calculate run time
	variable vTotalSec = gTotalSec					// changing to local variable
	variable vNumPtsAtTS 	= vTotalSec/gDeltaT	// number of points at time step
	variable vNumPtsAt1Hz	= vTotalSec				// number of points at 1 Hz
	variable iT1Hz											// loop variable for points at 1 Hz
	Variable vNumiTPerSec	= 1/gDeltaT					// time steps per second
	variable vCurrentTime1Hz							// used to calculate current time at 1 Hz from time step
	variable vPassivate = (round(gOG/gDeltaT))  // time step when emission ends
	
	// Checking alpha, which is a measure of model approximation stability
	variable vDf  = 10^gDf // (cm2/sec)
	variable valpha =  vDf * gDeltaT / (gDeltaY/1e4)^2	
	print "alpha "+num2str(valpha)
	If(valpha > 0.5)
		Abort "Aborting: Dimensionless difussion coefficient is greater than 0.5, model is unstable." 
	Endif 
		
	// Calculate required thickness to model based on diffusion coefficient
	vThickness = sqrt(2*vDf*vTotalSec)*4*1e4 // (um)
	
	// Loop variables
	variable nLayers = floor(vThickness/gDeltaY)	// Number of layers, rounded down to a whole number
	variable nBins 	= floor(gLength/gDeltaX)		// Number of bins, rounded down to whole number
	variable iT, iB, iL
	Variable vCond, vEvap, vinput, vexhaust
	
	// Calculating Cw per box (ug of wall per box / vol of air, m3)
	variable vCircumference = pi*vDiam // (m)
	variable vVolAir = pi*(vDiam/2)^2*vDeltaX // (m3) volume of air in bin i
	variable vCwL  = vCircumference*vDeltaX*vDeltaY*vDensity/vVolAir*gGamma  // (ug/m3) volume of box size * density of material / volume of air * availability parameter

	//===ka parameterization from McMurry & Stolzenburg 1987, with multiplier to account for 
	//===surface roughness of paint film, as explained in Algrim et al. 2019 text. 
	Variable vka  	= gSR*3*gDfG / ((gDiam)^2) 		// absorption rate constant 
	Variable vKgw 	= vCwL / gCStar 					// Gas wall equilibrium coef. for tube 
	variable vkd  	= vKa / vKgw				 			// desorption rate constant
	variable gV 		= pi*(gDiam/2)^2*gLength*0.001	// Volume of tubing in liters
	variable kflow 	= gFin/60/(gV/nBins)       		// Inverse of residence time (seconds^-1)
	
	// END OF LOCAL VARIABLE CREATION
   // ============================================================================
	// ============================================================================

	// ============================================================================
	// ============================================================================
	// START OF LOCAL WAVE CREATION

	// matrix of gas conc, row=Time, column=Bin, 1 Hz
	string sConcG = "wConcG_"+str
	make/o/n=(vNumPtsAt1Hz, nBins)	$sConcG =0
	wave wConcG = $sConcG  		
	
	// temp waves for Euler method, gas conc per Bin
	Make/o/n=(nBins) 	wTempConcG = 0
	Make/o/n=(nBins) 	wTempConcGprev = 0
	
	// temp waves for Euler method, wall conc per layer, bin
	string sConcW_Pass = "wConcW_Pass_"+str
	string sConcW_Depass = "wConcW_Depass_"+str
	make/o/n=(nLayers, nBins)  $sConcW_Pass = 0
	make/o/n=(nLayers, nBins)  $sConcW_Depass = 0
	wave wTempConcW_Pass = $sConcW_Pass // extra matrix that I will store the passivated wConcW matrix in
	wave wTempConcW = $sConcW_Depass	// the matrix that updates after every time step - will finish at the depassivated state
	make/o/n=(nLayers, nBins)  wTempConcWprev = 0
	
	// time wave 1 Hz
	Make/o/d/n=(vNumPtsAt1Hz) wTime = p

	
	// gas conc in last bin at all times - plot this for time series
	string sTubeEndConc = "wTubeEndConc_" + str
	make/o/n=(vNumPtsAt1Hz) $sTubeEndConc =0
	wave wTubeEndConc = $sTubeEndConc
	
	// END OF LOCAL WAVE CREATION
   // ============================================================================
	// ============================================================================

   // ============================================================================
	// ============================================================================
	// START MODEL
		
	//==================Loop Time===================================
	
	For(iT = 0; iT <vNumPtsAtTS; iT+=1)
	
		//================Loop Bins==================================
	
		For(iB = 0; iB < nBins; iB+=1)
	
			If(iB ==0  && iT < vPassivate) 		//at first bin , passivating 
				vCin = gCin 
			elseif (iB ==0 && iT >= vPassivate) 	//at first bin, depassivating 
				vCin = 0
			else							    		   	//All other bins, Cin is based on previous Bin concentration
				vCin  = wTempConcG[iB-1]
			endif
			  
			If(iT ==0)			// initialize conditions if time is zero
				vCond =0 
				vEvap = 0
				vExhaust = 0 
			
			else 					// calculate conditions if time is not zero
				vinput    =    		        vCin   	* kflow  * gDeltaT    
				vExhaust  = 		wTempConcGprev[iB]	* kflow  * gDeltaT
				vCond     = 		wTempConcGprev[iB]	* vKa    * gDeltaT    
				vEvap	    = wTempConcWprev[0][iB]	* vKd    * gDeltaT   		
			endif 
			
			wTempConcG[iB] = wTempConcGprev[iB] + vInput - vExhaust - vCond + vEvap
			
			
			//===============Loop Layers================================


			for(iL =0 ;iL<(nLayers); iL += 1)	
				If( iL == 0)						//Surface layer calculations
					wTempConcW[iL][iB] = wTempConcWprev[iL][iB] + vCond - vEvap - (vAlpha*(wTempConcWprev[iL][iB] - wTempConcWprev[iL+1][iB]))
				Elseif (iL == nLayers -1)  	//Last layer calculations
					wTempConcW[iL][iB] = wTempConcWprev[iL][iB] + (vAlpha*(wTempConcWprev[iL-1][iB] - wTempConcWprev[iL][iB]))
				else 								//All other layers calculations							
					wTempConcW[iL][iB] = wTempConcWprev[iL][iB] + (vAlpha*(wTempConcWprev[iL-1][iB] -(2*wTempConcWprev[iL][iB]) + wTempConcWprev[iL+1][iB]) )
				endif 
			endfor
			
			//===============End Loop Layers================================
		
		endfor
		
		//===============End Loop bins================================
		
		wTempConcGprev = wTempConcG
		wTempConcWprev = wTempConcW
		
		// Record data at 1 Hz
		if(mod(iT,vNumiTPerSec)==0)  // if remainder=0 for iT/vNumiTPerSec, record data for all bins
			
			vCurrentTime1Hz = iT/vNumiTPerSec  // calculate current time at 1 Hz
			
			// Loop through all bins and copy wTempConcG to wConcG
			for( iB = 0 ; iB < nBins ; iB += 1)
				wConcG[vCurrentTime1Hz][iB] = wTempConcG[iB]
			endfor						
		
		endif
				
		// Store wConcW_Pass at the last passivation time period
		if(iT == vPassivate)
			wTempConcW_Pass = wTempConcW
		endif
		
	endfor
	
	//===============End Loop Time================================
	
	// END MODEL
	// ============================================================================
	// ============================================================================
	
	// ============================================================================
	// ============================================================================
	// START CALCULATIONS FROM MODEL
	
	// Populate the tube end time series (gas conc from last bin at all times, 1 Hz)
	For(iT1Hz=0; iT1Hz < vNumPtsAt1Hz; iT1Hz +=1)
		wTubeEndConc[iT1Hz] = wConcG[iT1Hz][nBins-1]
	Endfor
	
	//	Display time series
		Display wTubeEndConc vs wTime
	
	// Print run time to command history
	Print (DateTime - StartTime ), "seconds to finish"
	
	saveexperiment
	
	// Calculate total Cw from penetration depth (10% of max concentration)
	wavestats/Q wTempConcW
	For(iL=0; iL < (nLayers-1); iL +=1)
		if(wTempConcW[nLayers-1-iL][0] >= v_Max*0.1)
			vCwTotal = (nLayers-1-iL) * nBins * vCwL // summed CwL to get Cwtotal
			vPenetrationDepth = nLayers-1-iL
			//			print "Penetration Depth = " + num2str(nLayers-1-iL) + " (um)"
			//			print "Cw = " + num2str(vCwTotal) + " (ug/m3)"
			// return [vCwTotal, vPenetrationDepth, vThickness]		
		endif		
	Endfor	
	 
	// END CALCULATIONS FROM MODEL
	// ============================================================================
	// ============================================================================
	
end 