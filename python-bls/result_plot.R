
# Read in true pars
seed_init = 5029879

true_alpha = sapply(seed_init:(seed_init+9), function(f) read.table(paste('../Data/pars_', f, '.csv', sep=''), header=FALSE, sep=',')[1,])
true_P = sapply(seed_init:(seed_init+9), function(f) read.table(paste('../Data/pars_', f, '.csv', sep=''), header=FALSE, sep=',')[2,])
true_td = sapply(seed_init:(seed_init+9), function(f) read.table(paste('../Data/pars_', f, '.csv', sep=''), header=FALSE, sep=',')[3,])
true_t0 = sapply(seed_init:(seed_init+9), function(f) read.table(paste('../Data/pars_', f, '.csv', sep=''), header=FALSE, sep=',')[4,])
true_rho = sapply(seed_init:(seed_init+9), function(f) read.table(paste('../Data/pars_', f, '.csv', sep=''), header=FALSE, sep=',')[5,])
true_sigma2 = sapply(seed_init:(seed_init+9), function(f) read.table(paste('../Data/pars_', f, '.csv', sep=''), header=FALSE, sep=',')[6,])
# Read in best fits
fake_P = sapply(seed_init:(seed_init+9), function(f) 
  scan(paste('../Results/out_', f, '.txt', sep=''), what=double()))

# Make plot
plot(true_P, fake_P, xlab="True Period", ylab="BLS Period Estimate")
dg <- par("usr") 
segments(dg[1],dg[3],dg[2],dg[4], col='red')