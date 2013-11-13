# Simple-Kepler Code
"compute_q" <- function(t,parameters) {
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

"dar1" <- function(y,rho,sigma2,log=FALSE) { # PROBLEM :: Mean is one not zero!! #TODO
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
                                                          a_rho=0.0,b_rho=1.0,nu_0=1.0,ss_0=1.0)) {
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
seed <- 5029879

gen_for_python = function(n, a_alpha, b_alpha, a_logP, b_logP, a_t_d, b_t_d, 
                          a_rho, b_rho, nu_0, ss_0, seed, plot_it=FALSE, save_it=TRUE) {
  
  # Set the random seed:
  set.seed(seed)
  y <- full_sim(n=n,hyperparameters=list(a_alpha=a_alpha, b_alpha=b_alpha, a_logP=a_logP,
                                         b_logP=b_logP, a_t_d=a_t_d, b_t_d=b_t_d, 
                                         a_rho=a_rho, b_rho=b_rho, nu_0=nu_0, ss_0=ss_0))
  
  #Filename is dependent on the seed. This is useful for writing multiple simulations to files. Just need to loop over
  #different values of the seed to generate unique files.
  if (plot_it) {
    jpeg(file=paste("Plots/Simulated_Data_", seed, ".jpg", sep=""), width=640, height=480)
    plot(y$data$y, ylab="Depth", xlab="Index", pch=1)
    lines(1-y$data$q,col="red",lwd=2.0)
    dev.off()
  }
  
  if (save_it) {
    write.table(y$data$y, file=paste("Data/y_", seed, ".csv", sep=""), sep=',', row.names=FALSE, col.names=FALSE)
    write.table(unlist(y$parameters), file=paste("Data/pars_", seed, ".csv", sep=""), sep=',', row.names=FALSE, col.names=FALSE)
  }
  ifelse(save_it, return(NULL), return(y)) #If saving to file, don't return anything. Else, return the object.
}

<<<<<<< HEAD
for (s in seed:(seed+999)) { #Generate 10 datasets using the given hyperpars.
  y = gen_for_python(n, a_alpha, b_alpha, a_logP, b_logP, a_t_d, b_t_d, a_rho, b_rho, nu_0, ss_0, seed=s, plot_it=FALSE, save_it=TRUE)
}
=======
for (s in seed:(seed+9)) { #Generate 10 datasets using the given hyperpars.
  gen_for_python(n, a_alpha, b_alpha, a_logP, b_logP, a_t_d, b_t_d, a_rho, b_rho, nu_0, ss_0, seed=s, save_it=TRUE)
}

>>>>>>> b47890c16c749a0376d053259408d186bb2caacd
