import bls
import numpy as np
from math import exp
from math import log

print "Reading data..."

data = np.genfromtxt('./obs.csv', delimiter=",",skiprows=1)

print "Setting up arrays..."

flux = np.empty(data.shape[0])
flux[:] = (data[:,1]).copy()
time = range(len(flux))
n = len(flux)
u = np.empty(n)
v = np.empty(n)

print "Forcing contiguity..."

flux = np.ascontiguousarray(flux,dtype=np.float64)
time = np.ascontiguousarray(time,dtype=np.float64)
u = np.ascontiguousarray(u,dtype=np.float64)
v = np.ascontiguousarray(v,dtype=np.float64)

print "Setting up BLS specs..."

# Prior for period:
a_logP = log(8000) # 8000
b_logP = log(20000) # 37000

min_period = exp(a_logP)
max_period = exp(b_logP)

###################
### Email Specs ####
nf = 10000
fmin = 1 / max_period
fmax = 1 / min_period
df = (fmax - fmin) / nf
nb = 1500
qmi = 0.0005
qma = 0.005
###################

print "Running EEBLS..."

results = bls.eebls(time, flux, u, v, nf, fmin, df, nb, qmi, qma)

print "Done. Results:"
print results

