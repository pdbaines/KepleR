import bls
import numpy as np
from math import exp
from math import log

flux = np.genfromtxt('./obs.csv', delimiter=",",skiprows=1)[:,1]
time = range(len(flux))
n = len(flux)
u = np.empty(shape=(n,1))
v = np.empty(shape=(n,1))

# Prior for period:
a_logP = log(8000) # 8000
b_logP = log(20000) # 37000

min_period = exp(a_logP)
max_period = exp(b_logP)

###################
### Email Specs ####
nf = 5000
fmin = 1 / max_period
fmax = 1 / min_period
df = (fmax - fmin) / nf
nb = 1500
qmi = 0.0005
qma = 0.005
###################

results = bls.eebls(time, flux, u, v, nf, fmin, df, nb, qmi, qma)
