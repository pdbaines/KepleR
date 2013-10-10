
setwd("~/Dropbox/GitHub/kepler/myscripts")
# Start anew:
rm(list=ls())

# Load Kepler functions:
source("kepler.R")

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

# Pseudo-BLS:
# -- Fold time series at each candidate period
# -- Fit step function to folded series
# -- For each (i1,i2)
# -- Compute mean(folded_y[i1:i2])
# -- Compute var(folded_y[i1:i2])
# -- Minimum variance for which (min_transit_depth < mean < max_transit_depth) is top candidate
# -- Repeat for all periods, i1, i2
# -- Rank candidates by variances

source("meebls.R")

cat("Calling MEEBLS...\n")
bls_out <- MEEBLS(y=y$data$y,t=c(1:n),nP=10000,Pmin=exp(a_logP),Pmax=exp(b_logP),
                  qmi=a_alpha,qma=b_alpha,fix_qrange=TRUE,nb=1000,nbmax=2000,verbose=TRUE)
cat("done.\n")

par(mfrow=c(3,1))
plot(y=y$data$y,x=(c(1:n)%%y$parameters$P))
plot(y=bls_out$Power_BLS,x=bls_out$Period_BLS) ; abline(v=y$parameters$P,col="blue",lwd=2.0)
plot(y=bls_out$Depth_BLS,x=bls_out$Period_BLS) ; abline(v=y$parameters$P,col="blue",lwd=2.0)

bail

# Load BLS:
#dyn.load("meebls.so")
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




