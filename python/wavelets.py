
import matplotlib.pyplot as plt
from matplotlib import cm
import pywt
import scipy
import numpy as np
import time
import random
import math
import sys
import kplr
import wavelet_funcs


data = np.loadtxt("KIC_7874976_llc.dat")
data.shape # (60020,4)
# From Montet:
# Columns are time, raw-flux, detrended-flux, uncertainty

# Plot detrended light curve
timepts = data[:,0]
y_raw = data[:,1]
y_det = data[:,2]

# 102 NA's:
np.sum(np.isnan(y_det))

plot(y_raw)
plot(timepts,y_det)

# Wavelet specs:
wavelet = 'db2'
mode = 'zpd'

# Wavelet transform:
wvb = mwavelets(y,wavelet=wavelet,mode=mode)

# Look at distributions of wavelets at each scale:

# Look at autocorrelations of wavelets at each scale:

# Look at QQ-plots of wavelets at each scale (check for normality)






















assert 0

###################################
# TODO:
# (1) Handle non-power of 2 series
# (2) Handle nan's
###################################

# Set the seed:
rng_seed = 3123194149
random.seed(rng_seed)

######################
# Simulate:
N = pow(2,10)
g_mean = 0.0
h_width = 8
t_depth = 1.0
sigma = 0.2 # 0.2
#######################

# Flat mean...
m = np.empty(N)
m.fill(g_mean)

# Noisy mean...
sin_cycles = 2.0
sin_height = 0.4
m = np.empty(N)
for i in range(0,N):
    m[i] = sin_height*math.sin(2*math.pi*i*sin_cycles/N)

# Layer in a transit...
for i in range(0,2*h_width):
    ix = (N/2)+i
    m[ix] = m[ix]-t_depth

# Simulate:
y = np.empty(N)
for i in range(0,N):
    y[i] = np.random.normal(loc=m[i],scale=sigma)

#####################################

# Begin wavelet stuff...
# pywt.families()
# ['haar', 'db', 'sym', 'coif', 'bior', 'rbio', 'dmey']
# pywt.wavelist('haar')
# ['haar'] # other families have lots...
#
# y_dwt = pywt.dwt(y, wavelet='haar', mode='zpd')
# print y_dwt
#
#pywt.wavedec(y, wavelet='haar', level=1)
#
#full_y_dwt = pywt.wavedec(y, wavelet='haar')
#coef_pyramid_plot(full_y_dwt)
#
#fdwt_vec = np.empty(0)
#for i in range(0,len(fdwt)):
#    fdwt_vec = np.append(fdwt_vec,fdwt[i])
#
#plt.plot(fdwt_vec)
#plt.show()
#
# cm.coolwarm
# cm.afmhot

#f_ret = wvanalyze(y,wtype='father',wavelet='haar',mode='zpd',cmap=cm.afmhot)

m_ret = wvanalyze(y,wtype='mother',wavelet='haar',mode='zpd',cmap=cm.afmhot)

assert 0

# Now, REAL Kepler data... :)

import kplr
import pyfits

data = np.loadtxt("transit.txt")
len(data)

data = data - np.median(data)
data = data[0:512]

f_ret = wvanalyze(data,wtype='father',wavelet='haar',mode='zpd',cmap=cm.afmhot)



client = kplr.API()

planet = client.planet("62b")

do_fetch = False

if do_fetch:
    all_datasets = []
    ndatasets = 0
    for dataset in planet.data:
        print "Fetching Kepler dataset " + str(ndatasets) + "..."
        all_datasets.append(dataset.fetch())
        ndatasets = ndatasets + 1
else:
    print "TODO: Read local file..."

# Filename (unicode):
all_datasets[0].filename

data = kplr.Dataset(all_datasets[0].filename)

# data.saplux is lightcurve...

foo = data.sapflux.copy()
sum(foo)

f_ret = wvanalyze(data.sapflux,wtype='father',wavelet='haar',mode='zpd',cmap=cm.afmhot)


assert 0



plt.plot(fgrid)
plt.show()


fig = plt.figure(1)
# Top plot:
ax1 = plt.subplot(2,1,1)
p1 = ax1.plot(y)
# Bottom plot:
ax2 = plt.subplot(2,1,2)
p2 = ax2.plot(y)
plt.show()

scales = wave.autoscales(N=y.shape[0], dt=1, dj=0.25, wf='dog', p=2)
X = wave.cwt(x=x, dt=1, scales=scales, wf='dog', p=2)
p2 = ax2.imshow(np.abs(X), interpolation='nearest')
plt.show()


