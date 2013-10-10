
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <math.h>

//////////////////////////////////////////////////////////////////////////////

void 
returnminmax(int ntop,
        double *value, 
        double *min, 
        double *max);

void 
R_meebls(SEXP R_NPTS, 
    SEXP R_nP, 
    SEXP R_Pmin, 
    SEXP R_Pmax, 
	SEXP R_qmi_fix, 
    SEXP R_qma_fix, 
    SEXP R_nb_fix, 
    SEXP R_fix_qrange,
	SEXP R_t_vec, 
    SEXP R_y_vec, 
    SEXP R_Power_BLS, 
	SEXP R_Period_BLS, 
    SEXP R_Depth_BLS, 
    SEXP R_InCare_BLS,
	SEXP R_PowerBright_BLS, 
    SEXP R_PeriodBright_BLS, 
	SEXP R_DepthBright_BLS, 
    SEXP R_InCareBright_BLS,
    SEXP R_nbmax,
	SEXP R_want_print);

void 
meebls(int NPTS, 
    double nP, 
    double Pmin, 
    double Pmax, 
	double qmi_fix, 
    double qma_fix, 
    int nb_fix, 
    int fix_qrange, 
	double * t_vec, 
    double * y_vec, 
    double * Power_BLS, 
	double * Period_BLS, 
    double * Depth_BLS, 
    double * InCare_BLS,
	double * PowerBright_BLS, 
    double * PeriodBright_BLS, 
	double * DepthBright_BLS, 
    double * InCareBright_BLS,
    int nbmax,
	int want_print);

//////////////////////////////////////////////////////////////////////////////

void 
R_meebls(SEXP R_NPTS, 
    SEXP R_nP, 
    SEXP R_Pmin, 
    SEXP R_Pmax, 
	SEXP R_qmi_fix, 
    SEXP R_qma_fix, 
    SEXP R_nb_fix, 
    SEXP R_fix_qrange,
	SEXP R_t_vec, 
    SEXP R_y_vec, 
    SEXP R_Power_BLS, 
	SEXP R_Period_BLS, 
    SEXP R_Depth_BLS, 
    SEXP R_InCare_BLS,
	SEXP R_PowerBright_BLS, 
    SEXP R_PeriodBright_BLS, 
	SEXP R_DepthBright_BLS, 
    SEXP R_InCareBright_BLS,
    SEXP R_nbmax,
	SEXP R_want_print)
{
	int NPTS = INTEGER(R_NPTS)[0];
	double nP = REAL(R_nP)[0];
	double Pmin = REAL(R_Pmin)[0];
	double Pmax = REAL(R_Pmax)[0];
	double qmi_fix = REAL(R_qmi_fix)[0];
	double qma_fix = REAL(R_qma_fix)[0];
	int nb_fix = INTEGER(R_nb_fix)[0];
	int fix_qrange = INTEGER(R_fix_qrange)[0];
	double * t_vec = REAL(R_t_vec);
	double * y_vec = REAL(R_y_vec);
	double * Power_BLS = REAL(R_Power_BLS);
	double * Period_BLS = REAL(R_Period_BLS);
	double * Depth_BLS = REAL(R_Depth_BLS);
	double * InCare_BLS = REAL(R_InCare_BLS);
	double * PowerBright_BLS = REAL(R_PowerBright_BLS);
	double * PeriodBright_BLS = REAL(R_PeriodBright_BLS);
	double * DepthBright_BLS = REAL(R_DepthBright_BLS);
	double * InCareBright_BLS = REAL(R_InCareBright_BLS);
    int nbmax = INTEGER(R_nbmax)[0];
	int want_print = INTEGER(R_want_print)[0];

	Rprintf("Inside 'R_meebls'...\n");
	Rprintf("want_print = %d\n",want_print);

	if (want_print){
		Rprintf("Calling meebls with the following arguments:\n");
		Rprintf("NPTS = %d\n",NPTS);
		Rprintf("nP = %g\n",nP);
		Rprintf("Pmin = %g\n",Pmin);
		Rprintf("Pmax = %g\n",Pmax);
	}

	meebls(NPTS,nP,Pmin,Pmax, 
	qmi_fix, qma_fix, nb_fix, fix_qrange, 
	t_vec,y_vec,Power_BLS,Period_BLS,Depth_BLS,InCare_BLS,
	PowerBright_BLS,PeriodBright_BLS,DepthBright_BLS,InCareBright_BLS,
    nbmax,want_print);

	return;
}

void 
meebls(int NPTS, 
    double nP, 
    double Pmin, 
    double Pmax, 
	double qmi_fix, 
    double qma_fix, 
    int nb_fix, 
    int fix_qrange, 
	double * t_vec, 
    double * y_vec, 
    double * Power_BLS, 
	double * Period_BLS, 
    double * Depth_BLS, 
    double * InCare_BLS,
	double * PowerBright_BLS, 
    double * PeriodBright_BLS, 
	double * DepthBright_BLS, 
    double * InCareBright_BLS,
    int nbmax,
	int want_print)
{

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Inputs:
//		NPTS				= number of data points
//		nP					= number of period points in which the spectrum is computed
//		Pmin				= minimum period investigated by the EEBLS routine 
//		Pmax				= maximum period investigated by the EEBLS routine 
//		t_vec				= array t_vec[i], containing the time values of the time series
//		y_vec				= array y_vec[i], containing the data values (in magnitudes) of the time series
//		want_print			= whether one wants the EEBLS routine to be verbose (want_print=1), or concise (want_print=0)
//	Outputs
//		Power_BLS			= array Power_BLS[jf] recording the Signal Residue returned by the EEBLS algorithm for dimming signals
//		Period_BLS			= array Period_BLS[jf] recording the period of the putative transit
//      Depth_BLS			= array Depth_BLS[jf] recording the depth of the putative transit
//      InCare_BLS			= array InCare_BLS[jf] recording the phase of the putative transit,
//      PowerBright_BLS		= array PowerBright_BLS[jf] recording the Signal Residue returned by the EEBLS algorithm for the putative brightening anti-transit
//      PeriodBright_BLS	= array Period_BLS[jf] recording the period of the putative brightening anti-transit
//      DepthBright_BLS		= array Depth_BLS[jf] recording the depth of the putative brightening anti-transit 
//      InCareBright_BLS	= array InCare_BLS[jf] recording the phase of the putative brightening anti-transit,
////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (want_print){
		Rprintf("Setting up...\n");
	}

	// Declare workspace arrays:
    // Note: these need to be larger because of the wrapping
	double UU_BLS[NPTS+NPTS];
	double VV_BLS[NPTS+NPTS];

    // Binning arrays:
	double ibi[nbmax];
	double y[nbmax];

	double qma,qmi;
	int good = 1;
	int minbin = 5;
	int nb;

	// Fixed transit fraction ranges?
	if (fix_qrange){
		qmi = qmi_fix;
		qma = qma_fix;
		nb = nb_fix;
	} 
	
	//return the minimum and maximum time of the data-set
	double jdmini,jdmaxi;
	returnminmax(NPTS,t_vec,&jdmini,&jdmaxi);
	double tot = jdmaxi-jdmini;
	
	//determine the appropriate logarithmic period spacing
	if (want_print){
		Rprintf("Determining the appropriate logarithmic period spacing...\n");
	}
	
    double eta = pow(Pmax/Pmin,(double)1.0/((double)nP-1.0)) - 1.0;
	
    if (want_print==1){
		Rprintf("eta: %lf\n",eta);
	}

	double rn = (double)NPTS;
	double bpow = 0;
	double bpow_bright = 0;	
	int jnb;

	//=================================
	// Pre-process time series
	//=================================

    // s = storage for mean(y_vec):
	double s = 0.0;
    
    // t1 == initial time point
	double t1 = t_vec[0];
	
    // UU_BLS = (t_vec - t1) # i.e., zero-origin time vector
    for (int i=0; i<NPTS; i++){
     	UU_BLS[i] = t_vec[i]-t1;
	    s += y_vec[i];
	}

	// s = sum(y_vec), so make s = mean(y_vec):
	s = s/rn;

    // VV_BLS = y_vec - mean(y_vec) # i.e., zero-mean observation vector
	for (int i=0; i<NPTS; i++){
		VV_BLS[i] = y_vec[i]-s;
	}

	//******************************
	//     Start period search     *
	//******************************
    
    // Useful storage variables...
	double p0,f0; 
	double ph;
	double power,pow_pow;
	double power_bright,pow_pow_bright;
	int jn1,jn2;
		
    // Regular versions for dimming:
	int stoper=0;
	double s3,k,kk,nb2,rn1,rn3;
	double in1,in2,qtran,depth,bper;
	double phase_care;

	// These are the versions for brightening
	int jn1_bright,jn2_bright;
	double rn3_bright,s3_bright,phase_care_bright;
	double in1_bright,in2_bright,qtran_bright,depth_bright,bper_bright;

    // How often to print to screen:
	double jfprint_double = 1000.0;

	if (want_print){
		Rprintf("Beginning loop over periods...\n");
	}

	for (int jf = 0; jf<nP; jf++){

		// Logarithmic period space as recommend in Burke et al. (2006)
		p0 = Pmin*pow( (double)(1.0+eta),(double)jf );
		f0 = 1.0/p0;
		if (!fix_qrange){
			error("'fix_qrange' must be TRUE for meebls");
		}

        // Check number of bins does not exceed the maximum:
		if (nb > nbmax){
			nb = nbmax;
		}

		if (jf/jfprint_double == (int)(jf/jfprint_double)){
			if (want_print==1){
				Rprintf("jf: %d/%d, nb: %d\n",jf,(int)nP,nb);
			}
		}
		//printf("fmin: %f,f0: %f, p0: %f\n",fmin,f0,p0);

        // Note: Full series is binned into nb bins.

        // Minimum number of bins for a transit:
		int kmi = (int)(qmi*(double)nb);
        // Check for degeneracy:
		if (kmi < 1){ 
			kmi = 1;
		}
        // Maximum number of bins for a transit
		int kma = (int)(qma*(double)nb)+1;
        // Total number of pts * minimum number of bins for a transit,
        // i.e., rough/average number of points in smallest transit bin:
		int kkmi = (int)(rn*qmi);
		if (kkmi < minbin){
			kkmi=minbin;
			Rprintf("warning setting to minbin\nwarning\nwarning\nwarning\n");
		}
        // Number of bins + 1:
		int nb1 = nb+1;
        // Number of bins + maximum number of bins for a transit:
		int nbkma = nb+kma;
        // Debug?
        if (want_print){
            Rprintf("(max,min) number of bins for a transit at period %g = (%d,%d)\n",p0,kmi,kma);
        }
		///////////////////////////////////////////

		//======================================================
		//     Compute folded time series with  *p0*  period
		//======================================================
        
        // Set y and ibi to zero 
        // These are the number and sum of points in each bin 
        // (Note: used to be 1-->nb, switched to 0-->nb-1):
		for (int j=0; j<nb; j++){
			y[j] = 0;
			ibi[j] = 0;
		}
		
        // Compute phase:
        // Recall UU_BLS = zero-origin tie, VV_BLS = zero-mean version of y_vec
		for (int i=0; i<NPTS; i++){
            // ph = fractional part of (t[i]/P_0)
            // i.e., how far into the period t[i] falls 
			ph     = UU_BLS[i]*f0;
			ph     = ph - (int)ph;
            // j = bin membership of data point i
            // j = (number of bins * fractional part of (t[i]/P_0))
			int j  = (int)(nb*ph);
            // Number of points in bin j:
			ibi[j] += 1;
            // Sum of points in bin j:
			y[j] += VV_BLS[i];
		}

        // Check number of bins is correct:
        if (j != nb){
            Rprintf("j=%d, nb=%d\n",j,nb);
            error("'j' and 'nb' do not match after loop");
        }

		// Extend the arrays ibi() and  y() beyond  nb  by  wrapping
        // Recall: nb1 = number of bins + 1
		for (int j=nb; j<nbkma; j++){
            //   j loop: nb,nb+1,...,nb+kma-1 
            // jnb loop: 0,1,2,...,kma-1 (i.e., maximum number of bins for a transit - 1)
			jnb = j-nb;
            // Mirror values e.g., y[0] into y[nb+0], y[1] into y[nb+1] etc.
			ibi[j] = ibi[jnb];
			y[j] = y[jnb];
		}

        // Now, y looks like:
        // TODO
        // Now ibi looks like:
        // TODO
        
		//===============================================
		// Compute BLS statistics for this period
		//===============================================
			
		power = 0.0;
		power_bright = 0.0;

        // TODO: Change to zero indexing

        // Loop over number of bins:
		for (int i=1; i<=nb; i++){//for i from 1 to nbins
			s     = 0;
			k     = 0;
			kk    = 0;
			nb2   = i+kma;//the bin we are looking at plus the maximum length of transit
			//if(nb2 > nb){ //this line is likely incorrect
			//printf("nb2: %lf/%d\n",nb2,nb);
			//	nb2=nb;
			//}
			for (int j=i; j<=nb2; j++){
				k     = k+1;
				kk    = kk+ibi[j];
				s     = s+y[j];
				//printf("k: %d, kmi: %d, kk: %d, kkmi: %d\n",k,kmi,kk,kkmi);
				if (k < kmi || kk < kkmi){
					continue;
				} else {
					rn1 = (double)kk;
						
					if( -s*rn/(rn1*(rn-rn1)) < 0){
                        // This is a dimming:
						pow_pow = s*s/(rn1*(rn-rn1));//HELP
						pow_pow_bright = 0;
					} else{
                        // This is a brightening:
						pow_pow_bright = s*s/(rn1*(rn-rn1));
						pow_pow = 0;
					}
						
					// Continue, if the SR is not above previous values
					if (pow_pow < power && pow_pow_bright < power_bright){
						continue;
					}
					// This is dimming:
					if (pow_pow >= power){
						power = pow_pow;
						//printf("%d,pow_pow: %le\n",j,pow_pow);
						jn1   = i;
						//printf("jn1: %d\n",jn1);
						jn2   = j;
						rn3   = rn1;
						s3    = s;
					}
					// This is BRIGHTENING:
					if (pow_pow_bright >= power_bright){
						power_bright = pow_pow_bright;
						//printf("%d,pow_pow: %le\n",j,pow_pow);
						jn1_bright   = i;
						//printf("jn1: %d\n",jn1);
						jn2_bright   = j;
						rn3_bright   = rn1;
						s3_bright    = s;
					} 

				}  // END if (k < kmi || kk < kkmi){ ... } else { ... }
			} // END for-loop over j
		} // END for-loop over bins

		//DIMMING
		power = sqrt(power);
		Power_BLS[jf] = power;
		Period_BLS[jf] = p0;
		Depth_BLS[jf] = -s3*rn/(rn3*(rn-rn3)); 
		phase_care = 0.5*(jn2/(double)nb + jn1/(double)nb);
		InCare_BLS[jf] = phase_care - (int)phase_care;//this makes sure the phase is between 0 and 1

		//BRIGHTENING
		power_bright = sqrt(power_bright);
		PowerBright_BLS[jf] = power_bright;
		PeriodBright_BLS[jf] = p0;
		DepthBright_BLS[jf] = -s3_bright*rn/(rn3_bright*(rn-rn3_bright)); 
		phase_care_bright = 0.5*(jn2_bright/(double)nb + jn1_bright/(double)nb);
		InCareBright_BLS[jf] = phase_care_bright - (int)phase_care_bright;//this makes sure the phase is between 0 and 1
	
		//
		if (power < bpow && power_bright < bpow_bright){
			continue;
		}

		if (power >= bpow){

			bpow  =  power;
			in1   =  jn1;
			in2   =  jn2;
			if (in2 > nb){ 
				//This is the edge effect 
				in2 = in2-nb;
			}    
			qtran =  rn3/rn;
			depth = -s3*rn/(rn3*(rn-rn3));
			bper  =  p0;
			if (want_print==1){
				Rprintf("%04d/%04d, DIMMING     period: %f,freq: %f,power: %f,depth: %f,phase: %lf,\n",
					jf,(int)nP,bper,f0,bpow,depth, InCare_BLS[jf] );
			}

		} // END if (power >= bpow){ ... } 

		if (power_bright >= bpow_bright){

			bpow_bright  =  power_bright;
			in1_bright   =  jn1_bright;
			in2_bright   =  jn2_bright;
			if (in2_bright > nb){
				//This is the edge effect
				in2_bright = in2_bright-nb;
			}     
			qtran_bright =  rn3_bright/rn;
			depth_bright = -s3_bright*rn/(rn3_bright*(rn-rn3_bright));
			bper_bright  =  p0;
			if (want_print==1){
				Rprintf("%04d/%04d, BRIGHTENING period: %f,freq: %f,power: %f,depth: %f,phase: %lf,\n",
					jf,(int)nP,bper_bright,f0,bpow_bright,depth_bright, InCareBright_BLS[jf] );
			}

		} // END if (power_bright >= bpow_bright){ ... } 
			
 	} // END loop over periods

	return;
}

void 
returnminmax(int ntop, 
        double *value, 
        double *min, 
        double *max)
{
	double minner = R_PosInf;
	double maxxer = R_NegInf;
	for (int j=0; j<ntop; j++){
		if (value[j] < minner){
			minner = value[j];
		}
		if (value[j] > maxxer){
			maxxer = value[j];
		}
	}
	*min = minner;
	*max = maxxer;
	return;
}

