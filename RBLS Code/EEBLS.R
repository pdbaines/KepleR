
# // Inputs:
# //		NPTS				= number of data points
# //		nP					= number of period points in which the spectrum is computed
# //		Pmin				= minimum period investigated by the EEBLS routine 
# //		Pmax				= maximum period investigated by the EEBLS routine 
# //		JD_BLS				= array JD_BLS[i], containing the time values of the time series
# //		MAG_BLS				= array MAG_BLS[i], containing the data values (in magnitudes) of the time series
# //		want_print			= whether one wants the EEBLS routine to be verbose (want_print=1), or concise (want_print=0)
# //	Outputs
# //		Power_BLS			= array Power_BLS[jf] recording the Signal Residue returned by the EEBLS algorithm for dimming signals
# //		Period_BLS			= array Period_BLS[jf] recording the period of the putative transit
# //      Depth_BLS			= array Depth_BLS[jf] recording the depth of the putative transit
# //      InCare_BLS			= array InCare_BLS[jf] recording the phase of the putative transit,
# //      PowerBright_BLS		= array PowerBright_BLS[jf] recording the Signal Residue returned by the EEBLS algorithm for the putative brightening anti-transit
# //      PeriodBright_BLS	= array Period_BLS[jf] recording the period of the putative brightening anti-transit
# //      DepthBright_BLS		= array Depth_BLS[jf] recording the depth of the putative brightening anti-transit 
# //      InCareBright_BLS	= array InCare_BLS[jf] recording the phase of the putative brightening anti-transit,

"EEBLS_Croll" <- function(y,t,nP,Pmin,Pmax,qmi,qma,nb,fix_qrange=FALSE,verbose=FALSE)
{
	#############
	MAX <- 100000
	NBMAX <- 2000
	dyn.load("EEBLSmodified.so")
	#############

	n <- length(y)
	if (length(t)!=n){
		stop("'y' and 't' must have same length")
	}
	if (n>=MAX){
		stop("'MAX' not large enough for length of 'y' and 't'")
	}
	zero_mean <- TRUE
	if (zero_mean){
		y <- y-mean(y)
	} 
	storage.mode(y) <- storage.mode(t) <- "double"
	Power_BLS <- vector("double",nP)
	Period_BLS <- vector("double",nP)
	Depth_BLS <- vector("double",nP)
	InCare_BLS <- vector("double",nP)
	PowerBright_BLS <- vector("double",nP)
	PeriodBright_BLS <- vector("double",nP)
	DepthBright_BLS <- vector("double",nP)
	InCareBright_BLS <- vector("double",nP)
	want_print <- as.integer(verbose)
	fix_qrange <- as.integer(fix_qrange)
	if (fix_qrange){
		storage.mode(qmi) <- "double"
		storage.mode(qma) <- "double"
		storage.mode(nb) <- "integer"
		if (nb>NBMAX){
			stop("'nb' exceeds maximum allowable number of bins")
		}
	}

	.Call("R_eebls_croll",
	as.integer(n),as.numeric(nP),as.numeric(Pmin),as.numeric(Pmax),
	as.numeric(qmi),as.numeric(qma),as.integer(nb),as.integer(fix_qrange),t,y, 
	Power_BLS, Period_BLS, Depth_BLS, InCare_BLS,
	PowerBright_BLS, PeriodBright_BLS, DepthBright_BLS, InCareBright_BLS, 
	want_print)

	ret <- list("Power_BLS"=Power_BLS, 
		        "Period_BLS"=Period_BLS,
		        "Depth_BLS"=Depth_BLS,
		        "InCare_BLS"=InCare_BLS,
		        "PowerBright_BLS"=PowerBright_BLS,
		        "PeriodBright_BLS"=PeriodBright_BLS,
		        "DepthBright_BLS"=DepthBright_BLS, 
		        "InCareBright_BLS"=InCareBright_BLS)

	cat("\n")

	return(ret)	
}

