

# Simple-Kepler Code

"compute_q" <- function(t,parameters)
{
	# Everything out of transit is zero
	ret <- rep(0.0,length(t))
	# Compute time mod period:
	t_mod_P <- t %% parameters$P
	# Check if the time point is in the transit (vectorized):
	in_transit <- (t_mod_P>parameters$t_0) & (t_mod_P<(parameters$t_0 + parameters$t_d))
	# If it is, update the return vector:
	ret[in_transit] <- parameters$alpha
	return(ret)
}

"rinvchisq" <- function(n,nu,scale){
	return(nu*scale/rchisq(n=n,df=nu))
}

"rar1" <- function(n,rho,sigma2){
	# Simulate from an AR(1) in uber-slow style:
	if (n<1 || !is.numeric(n)){
		stop("'n' should be a positive number")
	}
	n <- as.integer(n)
	Y <- rep(NA,n)
	sigma <- sqrt(sigma2)
	Y[1] <- rnorm(n=1,mean=0,sd=sigma)
	if (n==1){
		return(Y)
	}
	for (i in 2:n){
		Y[i] = rnorm(n=1,mean=rho*Y[i-1],sd=sigma)
	}
	return(Y)
}

# PROBLEM :: Mean is one not zero!! #TODO
"dar1" <- function(y,rho,sigma2,log=FALSE)
{
	# Compute density of AR(1) process in uber-slow style:
	n <- length(y)
	sigma <- sqrt(sigma2)
	if (n==1){
		return(dnorm(y,sd=sigma,log=log))
	}
	omean <- c(0,rho*y[1:(n-1)])
	pieces <- dnorm(x=y,mean=omean,sd=sigma,log=log)
	if (log){
		ret <- sum(pieces)
	} else {
		ret <- prod(pieces)
	}
	return(ret)
}

"sim_from_prior" <- function(nsims=1,hyperparameters=list(a_alpha=0,b_alpha=1,
	a_logP=8000,b_logP=37000,a_t_d=15,b_t_d=30,
	a_rho=0.0,b_rho=1.0,nu_0=1.0,ss_0=1.0)){
	alpha <- runif(n=nsims,min=hyperparameters$a_alpha,max=hyperparameters$b_alpha)
	P <- exp(runif(n=nsims,min=hyperparameters$a_logP,max=hyperparameters$b_logP))
	t_d <- runif(n=nsims,min=hyperparameters$a_t_d,max=hyperparameters$b_t_d)
	t_0 <- runif(n=nsims,min=0,max=P-t_d)
	rho <- runif(n=nsims,min=hyperparameters$a_rho,max=hyperparameters$b_rho)
	sigma2 <- rinvchisq(n=nsims,nu=hyperparameters$nu_0,scale=hyperparameters$ss_0)
	parameters <- list("alpha"=alpha,"P"=P,"t_d"=t_d,"t_0"=t_0,"rho"=rho,"sigma2"=sigma2)
	return(parameters)
}

"dprior" <- function(priors,hyperparameters,log=TRUE){
	p.a <- dunif(priors$alpha,min=hyperparameters$a,max=hyperparameters$b,log=log)			
	return()
}

"sim_data" <- function(n,parameters){
	# Given parameters, simulate data
	t <- c(1:n)
	q_vec <- compute_q(t=t,parameters=parameters)
	z <- rar1(n=n,rho=parameters$rho,sigma2=parameters$sigma2)
	y <- (1-q_vec) + z
	return(list("y"=y,"z"=z,"q"=q_vec))
}


"full_sim" <- function(n,hyperparameters){
	parameters <- sim_from_prior(nsims=1,hyperparameters=hyperparameters)
	data <- sim_data(n=n,parameters=parameters)
	return(list("data"=data,"parameters"=parameters))
}

"mcmc_update" <- function(y,parameters,hyperparameters,tuning_parameters)
{
	# Preliminaries:
	n <- length(y)
	t <- c(1:n)
	q_vec <- compute_q(t=t,parameters=parameters)
	stats <- rep(0,length(parameters))
	names(stats) <- names(parameters)

	# Update sigma^{2}:
	nu_n <- hyperparameters$nu_0 + n
	mvec <- (y-(1-q_vec))
	mvec <- (mvec[2:n] - parameters$rho*mvec[1:(n-1)])
	scale_n <- (hyperparameters$nu_0*hyperparameters$ss_0 + sum(mvec^2))/nu_n
	parameters$sigma2 <- rinvchisq(n=1,nu=nu_n,scale=scale_n)
	stats["sigma2"] <- 1

	# Update \rho:
	rho_curr <- parameters$rho
	rho_prop <- rnorm(n=1,mean=rho_curr,sd=tuning_parameters$rho_prop_sd)
	if (rho_prop<hyperparameters$a_rho || rho_prop>hyperparameters$b_rho){
		# Reject
	} else {
		# Compute density (no prior term, since flat prior):
		mvec <- y-(1-q_vec)
		log_p_prop <- dar1(y=mvec,rho=rho_prop,sigma2=parameters$sigma2,log=TRUE)
		log_p_curr <- dar1(y=mvec,rho=rho_curr,sigma2=parameters$sigma2,log=TRUE)
		# Symmetric so no transition terms:
		log_MH_alpha <- log_p_prop - log_p_curr
		if (log(runif(1)) < log_MH_alpha){
			# Accept:
			parameters$rho <- rho_prop
			stats["rho"] <- 1
		}
	}

	# Update (t_0,t_d,P,\alpha):
	alpha_curr <- parameters$alpha
	alpha_prop <- runif(n=1,min=hyperparameters$a_alpha,max=hyperparameters$b_alpha)
	# rnorm(n=1,mean=alpha_curr,sd=tuning_parameters$alpha_prop_sd)
	if ((alpha_prop<hyperparameters$a_alpha) || (alpha_prop>hyperparameters$b_alpha)){
		# Reject
	} else {
		prop_parameters <- parameters
		prop_parameters$alpha <- alpha_prop
		q_curr <- compute_q(t=t,parameters=parameters)
		q_prop <- compute_q(t=t,parameters=prop_parameters)
		mvec_prop <- y-(1-q_prop)
		mvec_curr <- y-(1-q_curr)
		log_p_prop <- dar1(y=mvec_prop,rho=parameters$rho,sigma2=parameters$sigma2,log=TRUE)
		log_p_curr <- dar1(y=mvec_curr,rho=parameters$rho,sigma2=parameters$sigma2,log=TRUE)
		log_MH_alpha <- log_p_prop - log_p_curr
		if (FALSE){
			cat("============================================\n")
			cat("alpha update:\n")
			cat("Current state:\n") ; print(unlist(parameters))
			cat("Proposed state:\n") ; print(unlist(prop_parameters))
			cat(paste0("log_p_curr = ",log_p_curr,"\n"))
			cat(paste0("log_p_prop = ",log_p_prop,"\n"))
			cat(paste0("log_MH_alpha = ",log_MH_alpha,"\n"))
			cat("============================================\n\n")
		}
		if (log(runif(1))<log_MH_alpha){
			# Accept:
			parameters <- prop_parameters
			stats["alpha"] <- 1
		}
	}

	t_0_curr <- parameters$t_0
	t_0_prop <- rnorm(n=1,mean=t_0_curr,sd=tuning_parameters$t_0_prop_sd)

	t_d_curr <- parameters$t_d
	t_d_prop <- rnorm(n=1,mean=t_d_curr,sd=tuning_parameters$t_d_prop_sd)

	P_curr <- parameters$P
	P_prop <- rnorm(n=1,mean=P_curr,sd=tuning_parameters$P_prop_sd)


	return(list("parameters"=parameters,"stats"=stats))
}

tune_MH_RW_sd <- function(acc_rate,prop_sd)
{
	if (acc_rate<0.05){
		sf <- 4.0
	}
	if ((acc_rate>=0.05)&&(acc_rate<0.20)){
		sf <- 2.0
	}
	if ((acc_rate>=0.20)&&(acc_rate<0.50)){
		sf <- 1.0
	}
	if ((acc_rate>=0.50)&&(acc_rate<0.70)){
		sf <- 0.5
	}
	if ((acc_rate>=0.70)&&(acc_rate<0.90)){
		sf <- 0.35
	}
	if ((acc_rate>=0.90)&&(acc_rate<0.95)){
		sf <- 0.30
	}
	if (acc_rate>=0.95){
		sf <- 0.25
	}
	return(prop_sd*sf)
}

"kepler_fit" <- function(y,
	hyperparameters,
	burnin=1000,
	nsamples=10000,
	print_every=1000,
	tune_every=100,
	rho_prop_sd=0.02,
	P_prop_sd=1000,
	t_0_prop_sd=1000,
	t_d_prop_sd=4,
	alpha_prop_sd=0.01,
	debug=TRUE)
{
	####
	# Fit simple Kepler model using MCMC:
	###

	# Starting states... simulate from prior
	parameters <- sim_from_prior(nsims=1,hyperparameters=hyperparameters)
	
	# Error checks:
	if (burnin<0){
		stop("'burnin' must be positive")
	}
	nsamples <- as.integer(nsamples)
	if (nsamples<=0){
		stop("'nsamples' must be a positive integer")
	}
	niter <- burnin + nsamples
	npars <- length(parameters) # (t_0,t_d,P,alpha) + (rho,sigma2) = 6
	parnames <- names(parameters)

	# Create draws object:
	draws <- matrix(NA,nrow=nsamples,ncol=npars)
	colnames(draws) <- parnames
	draws <- mcmc(draws)

	# Store acceptance rates:
	acceptances <- rep(0,npars)
	names(acceptances) <- parnames

	# Tuning parameters:
	tuning_parameters <- 
		list("rho_prop_sd"=rho_prop_sd,
		     "P_prop_sd"=P_prop_sd,
		     "alpha_prop_sd"=alpha_prop_sd,
		     "t_0_prop_sd"=t_0_prop_sd,
		     "t_d_prop_sd"=t_d_prop_sd)

	# Begin the sampler:
	for (iter in 1:niter){

		# Update the parameters:
		update <- mcmc_update(y=y,
			parameters=parameters,
			hyperparameters=hyperparameters,
			tuning_parameters=tuning_parameters)

		parameters <- update$parameters
		acceptances <- acceptances + update$stats

		# Store the sample?
		if (iter>burnin){
			# Yup. Store the sample:
			draws[iter-burnin,] <- unlist(parameters)
		}

		# Tune the proposals?
		if ((iter%%tune_every == 0) && iter<burnin){
			# Acceptance rates since last check:
			acc_rates <- acceptances/tune_every
			# Tune according to some rules:
			tuning_parameters$rho_prop_sd <- 
				tune_MH_RW_sd(acc_rate=acc_rates["rho"],
					prop_sd=tuning_parameters$rho_prop_sd)
			# Period of transit:
			tuning_parameters$P_prop_sd <- 
				tune_MH_RW_sd(acc_rate=acc_rates["P"],
					prop_sd=tuning_parameters$P_prop_sd)
			# Transit depth:
			tuning_parameters$alpha_prop_sd <- 
				tune_MH_RW_sd(acc_rate=acc_rates["alpha"],
					prop_sd=tuning_parameters$alpha_prop_sd)
			# Initial transit time:
			tuning_parameters$t_0_prop_sd <- 
				tune_MH_RW_sd(acc_rate=acc_rates["t_0"],
					prop_sd=tuning_parameters$t_0_prop_sd)
			# Transit duration:
			tuning_parameters$t_d_prop_sd <- 
				tune_MH_RW_sd(acc_rate=acc_rates["t_d"],
					prop_sd=tuning_parameters$t_d_prop_sd)	
			# Reset acceptance rates until next time:
			acceptances[] <- 0
		}

		# Reset to get correct post-burnin acceptance rates:
		if (iter == burnin){
			acceptances[] <- 0
		}

		# Check status and print to user:
		if ((iter%%print_every) == 0){
			cat(paste0("Finished iteration ",iter,"...\n"))
		}
	}

	# Final acceptance rates:
	cat("Final acceptance rates:\n")
	print(round(acceptances/nsamples,4))

	if (debug){
		return(list("draws"=draws,
			"tuning"=tuning_parameters))
	}

	return(draws)
}



# n <- 1028
# rho <- 0.9
# sigma2 <- 10.0
# y <- rar1(n=n,rho=rho,sigma2=sigma2)
# plot(y,type="l")
# # lag-1 acf = rho
# acf(y,plot=F)
# print(rho)

# Number of data/time points:
n <- 48*365*3 # ~52k

# Prior for transit depth:
a_alpha <- 0.03 # 0.01
b_alpha <- 0.08 # 0.05

# Prior for period:
a_logP <- log(8000) # 8000
b_logP <- log(20000) # 37000

# Prior for transit duration:
a_t_d <- 50
b_t_d <- 80

# Prior for \rho:
a_rho <- 0.0
b_rho <- 1.0

# Prior for sigma^{2}:
nu_0 <- 1000.0
ss_0 <- 0.0001

# Random seed:
seed <- 5029879#2101323

hyperparameter_list <- 
	list(a_alpha=a_alpha,
		 b_alpha=b_alpha,
		 a_logP=a_logP,
		 b_logP=b_logP,
		 a_t_d=a_t_d,
		 b_t_d=b_t_d,
		 a_rho=a_rho,
		 b_rho=b_rho,
		 nu_0=nu_0,
		 ss_0=ss_0)

# Set the random seed:
set.seed(seed)

sim_time <- system.time({
	y <- full_sim(n=n,hyperparameters=hyperparameter_list)
})
cat(paste0("Simulation time (n=",n,"): ",round(sim_time["elapsed"],3)," seconds\n"))

plot(y$data$y)
lines(1-y$data$q,col="red",lwd=2.0)

###################
### MCMC Specs ####
nsamples <- 100#1000
burnin <- 10#100
tune_every <- 10
print_every <- 50
rho_prop_sd <- 0.05
P_prop_sd <- 1000
t_0_prop_sd <- 1000
t_d_prop_sd <- 4
alpha_prop_sd <- 0.01
debug <- TRUE
###################

"bls" <- function(t,y,nf,fmin,df,nb,qmi,qma)
{

# c     Input parameters:
# c     ~~~~~~~~~~~~~~~~~
# c
# c     n    = number of data points
# c     t    = array {t(i)}, containing the time values of the time series
# c     x    = array {x(i)}, containing the data values of the time series
# c     u    = temporal/work/dummy array, must be dimensioned in the 
# c            calling program in the same way as  {t(i)}
# c     v    = the same as  {u(i)}
# c     nf   = number of frequency points in which the spectrum is computed
# c     fmin = minimum frequency (MUST be > 0)
# c     df   = frequency step
# c     nb   = number of bins in the folded time series at any test period       
# c     qmi  = minimum fractional transit length to be tested
# c     qma  = maximum fractional transit length to be tested
# c
# c     Output parameters:
# c     ~~~~~~~~~~~~~~~~~~
# c
# c     p    = array {p(i)}, containing the values of the BLS spectrum
# c            at the i-th frequency value -- the frequency values are 
# c            computed as  f = fmin + (i-1)*df
# c     bper = period at the highest peak in the frequency spectrum
# c     bpow = value of {p(i)} at the highest peak
# c     depth= depth of the transit at   *bper*
# c     qtran= fractional transit length  [ T_transit/bper ]
# c     in1  = bin index at the start of the transit [ 0 < in1 < nb+1 ]
# c     in2  = bin index at the end   of the transit [ 0 < in2 < nb+1 ]

	n <- as.integer(length(y))
	if (length(t)!=n){
		stop("'t' and 'y' must have the same length")
	}
	mzy <- y - mean(y)
	# Workspace:
	u <- vector("double",n)
	v <- vector("double",n)
	# Outputs:
	p <- vector("double",nf)
	bper <- bpow <- depth <- qtran <- 0.0
	storage.mode(bper) <- storage.mode(bpow) <- storage.mode(depth) <- storage.mode(qtran) <- "double"
	in1 <- in2 <- 0
	storage.mode(in1) <- storage.mode(in2) <- "integer"
	.Fortran("eebls",
		as.integer(n),
		as.numeric(t),
		as.numeric(mzy),
		as.numeric(u),
		as.numeric(v),
		as.integer(nf),
		as.numeric(fmin),
		as.numeric(df),
		as.integer(nb),
		as.numeric(qmi),
		as.numeric(qma),
		as.numeric(p),
		as.numeric(bper),
		as.numeric(bpow),
		as.numeric(depth),
		as.numeric(qtran),
		as.integer(in1),
		as.integer(in2))

	ret <- list(p,bper,bpow,depth,qtran,in1,in2)
	
	return(ret)
}

min_period = exp(a_logP)
max_period = exp(b_logP)
freq_step = 1.1/n  #0.0001
min_duration = a_t_d
max_duration = b_t_d

df = freq_step
nf =  1000 # (1.0/min_period)/freq_step
nb = 1400 # < 2000
qmi = a_alpha
qma = b_alpha
fmin = 1.0/n # (max_period*1.1)

# Load BLS:
#dyn.load("eebls.so")
#out_bls <- bls(t=c(1:n),y=y$data$y,nf=nf,fmin=fmin,df=df,nb=nb,qmi=qmi,qma=qma)
#bail ## 

library(coda)

do_mcmc <- FALSE

if (do_mcmc){

mcmc_time <- system.time({
	ret <- kepler_fit(y=y$data$y,
		nsamples=nsamples,
		burnin=burnin,
		hyperparameters=hyperparameter_list,
		rho_prop_sd=rho_prop_sd,
		alpha_prop_sd=alpha_prop_sd,
		t_0_prop_sd=t_0_prop_sd,
		t_d_prop_sd=t_d_prop_sd,
		P_prop_sd=P_prop_sd,
		print_every=print_every,
		tune_every=tune_every,
		debug=debug)
})
cat(paste0("MCMC time (n=",n,"): ",round(mcmc_time["elapsed"],3)," seconds\n"))

if (debug){
	plot(ret$draws[,c("rho","sigma2")])
	print(summary(ret$draws))
} else {
	plot(ret[,c("rho","sigma2")])
	print(summary(ret))
}

} # END if (do_mcmc){...}

cat("Truth:\n")
print(unlist(y$parameters))

#plot(ret)

# Pseudo-BLS:
# -- Fold time series at each candidate period
# -- Fit step function to folded series
# -- For each (i1,i2)
# -- Compute mean(folded_y[i1:i2])
# -- Compute var(folded_y[i1:i2])
# -- Minimum variance for which (min_transit_depth < mean < max_transit_depth) is top candidate
# -- Repeat for all periods, i1, i2
# -- Rank candidates by variances

"pbls" <- function(y,
	P_lo,P_hi,P_inc,
	a_alpha,b_alpha,
	a_t_d,b_t_d
	)
{
	y <- y-mean(y)
	P_candidates <- seq(P_lo,P_hi,by=P_inc)
	ret <- matrix(NA,nrow=length(P_candidates),ncol=5)
	colnames(ret) <- c("P","t_0","t_d","alpha","var")
	ret[,1] <- P_candidates
	for (i in 1:length(P_candidates)){
		# Do stuff...
		P <- P_candidates[i]
		max_i1 <- as.integer(P-a_t_d)
		min_tvar <- Inf
		for (i1 in 1:max_i1){
			min_i2 <- i1+a_t_d
			max_i2 <- min(i1+b_t_d,n)
			for (i2 in min_i2:max_i2){
				# Do stuff.. ################# NO FOLDING!!!!! #############
				tbar <- mean(y[i1:i2])
				tvar <- var(y[i1:i2])
				if (tbar<a_alpha || tbar>b_alpha){
					# carry on
				} else if (tvar < min_tvar){
					ret[i,2:5] <- c(i1,(i2-i1),tbar,tvar)
				}
			}
		}
		cat(paste0("Finished P = ",P,"...\n"))
	}
	return(ret)
}

#crude_pbls <- pbls(y=y$data$y,
#	P_lo=exp(a_logP),P_hi=exp(b_logP),P_inc=6000,
#	P_lo=as.integer(y$parameters$P),P_hi=as.integer(y$parameters$P)+1,P_inc=1,
#	a_alpha=a_alpha,b_alpha=b_alpha,
#	a_t_d=a_t_d,b_t_d=b_t_d)


source("EEBLS.R")

cat("Calling EEBLS Croll...\n")
eebls_croll <- EEBLS_Croll(y=y$data$y,t=c(1:n),nP=10000,Pmin=exp(a_logP),Pmax=exp(b_logP),
	qmi=a_alpha,qma=b_alpha,fix_qrange=TRUE,nb=1000,verbose=TRUE)
cat("done.\n")

par(mfrow=c(3,1))
plot(y=y$data$y,x=(c(1:n)%%y$parameters$P))
plot(y=eebls_croll$Power_BLS,x=eebls_croll$Period_BLS) ; abline(v=y$parameters$P,col="blue",lwd=2.0)
plot(y=eebls_croll$Depth_BLS,x=eebls_croll$Period_BLS) ; abline(v=y$parameters$P,col="blue",lwd=2.0)


