
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <math.h>

#define MAX 100000

extern "C" 
{

//////////////////////////////////////////////////////////////////////////////

void returnminmax(int ntop,double value[MAX], double *min, double *max);
void return_transit_time_q(double period, double *qmi, double *qma, int *nb);

void 
R_eebls_croll(SEXP R_NPTS, SEXP R_nP, SEXP R_Pmin, SEXP R_Pmax, 
	SEXP R_qmi_fix, SEXP R_qma_fix, SEXP R_nb_fix, SEXP R_fix_qrange,
	SEXP R_JD_BLS, SEXP R_MAG_BLS, SEXP R_Power_BLS, 
	SEXP R_Period_BLS, SEXP R_Depth_BLS, SEXP R_InCare_BLS,
	SEXP R_PowerBright_BLS, SEXP R_PeriodBright_BLS, 
	SEXP R_DepthBright_BLS, SEXP R_InCareBright_BLS,
	SEXP R_want_print);

void 
eebls_croll_modified(int NPTS, double nP, double Pmin, double Pmax, 
	double qmi_fix, double qma_fix, int nb_fix, int fix_qrange, 
	double JD_BLS[MAX], double MAG_BLS[MAX], double Power_BLS[MAX], 
	double Period_BLS[MAX], double Depth_BLS[MAX], double InCare_BLS[MAX],
	double PowerBright_BLS[MAX], double PeriodBright_BLS[MAX], 
	double DepthBright_BLS[MAX], double InCareBright_BLS[MAX],
	int want_print);

//////////////////////////////////////////////////////////////////////////////

void 
R_eebls_croll(SEXP R_NPTS, SEXP R_nP, SEXP R_Pmin, SEXP R_Pmax, 
	SEXP R_qmi_fix, SEXP R_qma_fix, SEXP R_nb_fix, SEXP R_fix_qrange,
	SEXP R_JD_BLS, SEXP R_MAG_BLS, SEXP R_Power_BLS, 
	SEXP R_Period_BLS, SEXP R_Depth_BLS, SEXP R_InCare_BLS,
	SEXP R_PowerBright_BLS, SEXP R_PeriodBright_BLS, 
	SEXP R_DepthBright_BLS, SEXP R_InCareBright_BLS,
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
	double * JD_BLS = REAL(R_JD_BLS);
	double * MAG_BLS = REAL(R_MAG_BLS);
	double * Power_BLS = REAL(R_Power_BLS);
	double * Period_BLS = REAL(R_Period_BLS);
	double * Depth_BLS = REAL(R_Depth_BLS);
	double * InCare_BLS = REAL(R_InCare_BLS);
	double * PowerBright_BLS = REAL(R_PowerBright_BLS);
	double * PeriodBright_BLS = REAL(R_PeriodBright_BLS);
	double * DepthBright_BLS = REAL(R_DepthBright_BLS);
	double * InCareBright_BLS = REAL(R_InCareBright_BLS);
	int want_print = INTEGER(R_want_print)[0];

	Rprintf("Inside 'eebls_crool_R'...\n");
	Rprintf("want_print = %d\n",want_print);

	if (want_print){
		Rprintf("Calling eebls_croll_modified with the following arguments:\n");
		Rprintf("NPTS = %d\n",NPTS);
		Rprintf("nP = %g\n",nP);
		Rprintf("Pmin = %g\n",Pmin);
		Rprintf("Pmax = %g\n",Pmax);
	}

	eebls_croll_modified(NPTS,nP,Pmin,Pmax, 
	qmi_fix, qma_fix, nb_fix, fix_qrange, 
	JD_BLS,MAG_BLS,Power_BLS,Period_BLS,Depth_BLS,InCare_BLS,
	PowerBright_BLS,PeriodBright_BLS,DepthBright_BLS,InCareBright_BLS,want_print);

	return;
}

void 
eebls_croll_modified(int NPTS, double nP, double Pmin, double Pmax, 
	double qmi_fix, double qma_fix, int nb_fix, int fix_qrange, 
	double JD_BLS[MAX], double MAG_BLS[MAX], double Power_BLS[MAX], 
	double Period_BLS[MAX], double Depth_BLS[MAX], double InCare_BLS[MAX],
	double PowerBright_BLS[MAX], double PeriodBright_BLS[MAX], 
	double DepthBright_BLS[MAX], double InCareBright_BLS[MAX],
	int want_print)
{
//          This routine is a modification of the EEBLS routine. For a description
//          of the modifications please see Croll et al. 2007, 658, 1328
//			A brief summary of the modifications made to the EEBLS routine are given below.
//           >>>>>>>>>>>> This routine computes the EEBLS spectrum <<<<<<<<<<<<<<
//     Please see Kovacs, Zucker & Mazeh 2002, A&A, Vol. 391, 369 for a description
//     of the original (BLS) implementation.
//     Also see http://www.konkoly.hu/staff/kovacs/ for a description of the BLS -> EEBLS modifications
//		For further details, please see the following papers:
//		Kovacs, Zucker & Mazeh 2002, A&A, Vol. 391, 369
//		Croll et al. 2007, 658, 1328
//      
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Summary of the modifications to the EEBLS routine
//
// A brief summary of the modifications made to the EEBLS routine to create eebls_croll_modified is presented here.
// For a full explanation of these details, please see Croll et al. (2007). 
// - Logarithmic period spacing is used, rather than even frequency spacing. This modification was suggested by Burke et al. (2006)
//   Logarithmic period spacing retains high sensitivity to transits with both short and long orbital periods
// - eebls_croll_modified records the Signal Residue (SR) of both dimming and brightening transits independently. This reduces
//   the noise, as it allows one to focus solely on astrophysically plausible transits. 
// - The minimum, and maximum fractional transit lengths were set to be variable, rather than constant as in the original formulation. 
//   This is because the expected fractional transit length, Qmp, is dependent on the orbital period. The small planet approximation 
//   of Mandel & Agol (2002) is used to generate the expected fractional transit length for each period of interest (this is not performed
//   in the code). In the code currently these values are set as Qmi = 0.75*Qmp, and Qma = 1.1*Qmp.
// - The number of bins that each phase plot was binned into was also set to be dependent on the orbital period. The rationale for this 
//   variable number of bins was again because the expected fractional transit length is dependent on the period. Thus it was decided that
//   the same number of bins was desired in the in-transit portion of a transit at each period. The number of bins was therefore
//   set as Nb = 20.0/Qmi.
// These extra steps that are involved in the eebls_croll_modified routine
// are of negligible computational duration, but in the opinion of B. Croll significantly increase
// the detection efficiency of the routine for continuous space-based photometry, when one is attempting
// to search for low signal-to-noise transits. 
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Inputs:
//		NPTS				= number of data points
//		nP					= number of period points in which the spectrum is computed
//		Pmin				= minimum period investigated by the EEBLS routine 
//		Pmax				= maximum period investigated by the EEBLS routine 
//		JD_BLS				= array JD_BLS[i], containing the time values of the time series
//		MAG_BLS				= array MAG_BLS[i], containing the data values (in magnitudes) of the time series
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

	if (want_print){
		Rprintf("Setting up...\n");
	}

	int nbmax = 2000;
	//DECLARE pointer arrays
	double * UU_BLS; UU_BLS = new double[NPTS+NPTS];//this need to be larger because of the wrapping
	double * VV_BLS; VV_BLS = new double[NPTS+NPTS];
	double * ibi; ibi = new double[nbmax];
	double * y; y = new double[nbmax];

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
	returnminmax(NPTS,JD_BLS,&jdmini,&jdmaxi);
	double tot = jdmaxi-jdmini;
	
	//determine the appropriate logarithmic period spacing
	if (want_print){
		Rprintf("Determining the appropriate logarithmic period spacing...\n");
	}
	double eta = pow( Pmax/Pmin, 1/double(nP-1) ) -1;
	if(want_print==1){
		Rprintf("eta: %lf\n",eta);
	}
	if(good){
		if(want_print==1){
			Rprintf("still good\n");
		}
		float rn = float(NPTS);
		
		//printf("kmi: %d, kma: %d, kkmi: %d\n",kmi,kma,kkmi);
		double bpow = 0;

		//BRIGHT
		double bpow_bright = 0;

		
		//c     The following variables are defined for the extension 
		//c     of arrays  ibi()  and  y()  [ see below ]
		
		int jnb;

		//=================================
		//     Set temporal time series
		//=================================
		double s = 0;
		double t1 = JD_BLS[0];
		//double t1 = JD_BLS[0];
		//double t1 = gJD_min;
		//printf("t1: %lf\n",t1);
		for(int i=0; i<NPTS; i++){
     		UU_BLS[i] = JD_BLS[i]-t1;
	      	s = s+MAG_BLS[i];
		}
		//s is just the mean
		s = s/rn;
		for (int i=0; i<NPTS; i++){
			VV_BLS[i] = MAG_BLS[i]-s;
			//printf("UU: %lf, VV: %f\n",UU_BLS[i],VV_BLS[i]);
		}

		//******************************
		//     Start period search     *
		//******************************
		double p0,f0; 
		double ph;
		double power,pow_pow;
		double power_bright,pow_pow_bright;//HELP_BRIGHT
		int jn1,jn2;
		
		int stoper=0;
		double s3,k,kk,nb2,rn1,rn3;
		double in1,in2,qtran,depth,bper;
		double phase_care;

		//These are the functions for brightening
		int jn1_bright,jn2_bright;
		double rn3_bright,s3_bright,phase_care_bright;
		double in1_bright,in2_bright,qtran_bright,depth_bright,bper_bright;
		double jfprint_double = 1000.0;

		if (want_print){
			Rprintf("Beggining loop over periods...\n");
		}

		for (int jf = 0; jf<nP; jf++){
			//Logarithmic period space as recommend in Burke et al. (2006)
			p0 = Pmin*pow( double(1+eta),double(jf) );
			f0 = 1/p0;
			if (!fix_qrange){
				return_transit_time_q(p0,&qmi,&qma,&nb);
				if (want_print){
					Rprintf("RTT :: (p0=%g, qmi=%g, qma=%g, nb=%d)\n",p0,qmi,qma,nb);
				}
			}
			if (nb > nbmax){
				nb = nbmax;
			}

			if (jf/jfprint_double == int(jf/jfprint_double)){
				if (want_print==1){
					Rprintf("jf: %d/%d, nb: %d\n",jf,int(nP),nb);
				}
			}
			//printf("fmin: %f,f0: %f, p0: %f\n",fmin,f0,p0);

			///////////////////////////
			//NOW PUT THESE AND DECLARE THEM EACH TIME
			
			int kmi = int(qmi*float(nb));//minimum # of bins for a transit
			if (kmi < 1){ 
				kmi=1;
			}
			int kma = int(qma*float(nb))+1;//maximum number of bins for a transit
			int kkmi = int(rn*qmi);//number of pts*minimum number of bins for a transit
			if (kkmi < minbin){
				kkmi=minbin;
				Rprintf("warning setting to minbin\nwarning\nwarning\nwarning\n");
			}
			int nb1 = nb+1;
			int nbkma = nb+kma;
			///////////////////////////////////////////


			//c======================================================
			//     Compute folded time series with  *p0*  period
			//======================================================
			for (int j=1; j<=nb; j++){
				y[j] = 0;
				ibi[j] = 0;
			}
			
			for (int i=0; i<NPTS; i++){
				ph     = UU_BLS[i]*f0;
				ph     = ph - int(ph);
				int j      = 1 + int(nb*ph);
				ibi[j] = ibi[j] + 1;
				y[j]  =   y[j] + VV_BLS[i];
			}

			//c-------------EEBLS-----------------------------
			//Extend the arrays  ibi()  and  y() beyond  
			//nb   by  wrapping
			for (int j=nb1; j<=nbkma; j++){//should this be less than or less than or equal
				jnb = j-nb;
				ibi[j] = ibi[jnb];
				y[j] =   y[jnb];
				//printf("ibi[%d]: %lf, y[%d]: %lf\n",j,ibi[j],j,y[j]);
			}
			//c-----------------------------------------------   


			//c===============================================
			//c     Compute BLS statistics for this period
			//c===============================================
			
			power = 0.0;
			power_bright = 0.0;//HELP_BRIGHT

			for (int i = 1; i<=nb; i++){//for i from 1 to nbins
				s     = 0;
				k     = 0;
				kk    = 0;
				nb2   = i+kma;//the bin we are looking at plus the maximum length of transit
				//if(nb2 > nb){ //this line is likely incorrect
				//printf("nb2: %lf/%d\n",nb2,nb);
				//	nb2=nb;
				//}
				for(int j=i; j<=nb2; j++){
					k     = k+1;
					kk    = kk+ibi[j];
					s     = s+y[j];
					//printf("k: %d, kmi: %d, kk: %d, kkmi: %d\n",k,kmi,kk,kkmi);
					if (k < kmi || kk < kkmi){
						continue;
					}
					else{
						rn1   = float(kk);
						
						if( -s*rn/(rn1*(rn-rn1)) < 0){//this is a dimming
							pow_pow = s*s/(rn1*(rn-rn1));//HELP
							pow_pow_bright = 0;
						}
						else{//THIS IS BRIGHTENING
							pow_pow_bright = s*s/(rn1*(rn-rn1));
							pow_pow = 0;
						}
						
						//Continue, if the SR is not above previous values
						if (pow_pow < power && pow_pow_bright < power_bright){
							continue;
						}
						//THIS IS DIMMING
						if(pow_pow >= power){
							power = pow_pow;
							//printf("%d,pow_pow: %le\n",j,pow_pow);
							jn1   = i;
							//printf("jn1: %d\n",jn1);
							jn2   = j;
							rn3   = rn1;
							s3    = s;
						}
						//This is BRIGHTENING
						if(pow_pow_bright >= power_bright){
							power_bright = pow_pow_bright;
							//printf("%d,pow_pow: %le\n",j,pow_pow);
							jn1_bright   = i;
							//printf("jn1: %d\n",jn1);
							jn2_bright   = j;
							rn3_bright   = rn1;
							s3_bright    = s;
						}
					}
				}
			}

			//DIMMING
			power = sqrt(power);
			Power_BLS[jf] = power;
			Period_BLS[jf] = p0;
			Depth_BLS[jf] = -s3*rn/(rn3*(rn-rn3)); 
			phase_care = 0.5*(jn2/float(nb)+jn1/float(nb));
			InCare_BLS[jf] = phase_care-int(phase_care);//this makes sure the phase is between 0 and 1

			//BRIGHTENING
			power_bright = sqrt(power_bright);
			PowerBright_BLS[jf] = power_bright;
			PeriodBright_BLS[jf] = p0;
			DepthBright_BLS[jf] = -s3_bright*rn/(rn3_bright*(rn-rn3_bright)); 
			phase_care_bright = 0.5*(jn2_bright/float(nb)+jn1_bright/float(nb));
			InCareBright_BLS[jf] = phase_care_bright-int(phase_care_bright);//this makes sure the phase is between 0 and 1

			
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
						jf,int(nP),bper,f0,bpow,depth, InCare_BLS[jf] );
				}
			}
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
						jf,int(nP),bper_bright,f0,bpow_bright,depth_bright, InCareBright_BLS[jf] );
				}
			}
			
		//end period search
 		}
	}
	else{
		Rprintf("Problem\n");
		Rprintf("waiting until you quit, this program needs to be shut-down\n");
	}
	delete UU_BLS; delete VV_BLS; delete ibi; delete y;
	return;
}

void return_transit_time_q(double period, double *qmi, double *qma, int *nb)
//return_transit_time_q uses an empirical fit to the length of transit for a specific
//exoplanetary system.
//From this empirical fit it calculates the recommended minimum, and maximum fractional
//transit length, and the recommended number of bins for the inputted period
{
	//Inputs
	//period - the period that is being investigated
	//Outputs
	//*qmi - the recommended minimum fractional transit length for the inputted period
	//*qma - the recommended maximum fractional transit length for the inputted period
	//*nb - the recommended number of bins for the inputted period
	//
	double Phase_Oversample = 20.0;
	//
	double a,b,c,d,e;
	//HD 209458
	a = 0.906;
	b = 0.0068;
	c = 0.015800;
	d = 0.176;
	e = 0.001580;
	
	double q_mid = exp(-a*period)*(b*period*period+e*period*period*period*period+d)+c;
	*qmi = q_mid*0.75;//0.75
	*qma = q_mid*1.1;//1.1
	
	//could add a small premium
	double premium = 0.7*exp( - period) + 1.0; // close to one
	*nb = int(premium*Phase_Oversample/(*qmi)); // roughly 20/qmi

	return;
}

void returnminmax(int ntop, double value[MAX], double *min, double *max){
	double minner = 9e99;
	double maxxer = -9e99;
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

} // END extern "C" {...}

