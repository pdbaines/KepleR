source("EEBLS.R")

###################
### Email Specs ####
nf = 5000 #From the email
Pmin=18000 #Upped the minimum period to get closer to the true data
Pmax=22000 #Might not be enough?

fmin=1/Pmax #From the email
fmax = 1/Pmin #From the email

df = (fmax - fmin)/nf #From the email
nb = 1500 
qmi = 0.0005
qma = 0.005
###################

cat("Calling EEBLS_Modified...\n")
bls_out <- EEBLS_Croll(y=y$data$y,t=c(1:n),nP=10000, Pmin, Pmax, qmi, qma, nb, fix_qrange=TRUE, verbose=TRUE)
cat("done.\n")
