
import matplotlib.pyplot as plt
from matplotlib import cm
import pywt
import scipy
import numpy as np
import pyfits
import time
import random
import math
import sys
import kplr
import os
import copy
import wavelet_funcs as wv

def get_kepler_data(file_name,field_name):
    fitsdata = pyfits.open(file_name)
    kplr_data = fitsdata[1].data.field(field_name)
    fitsdata.close()
    return kplr_data

def rename_fits_to_txt(file_name):
    if (file_name[-5:] == ".fits"):
        main_part = file_name[:-5]
        new_file_name = main_part + ".txt"
        return new_file_name
    else:
        return None

# Basic fits info if done manually:
#klc.info()
#print(klc[1].data)
#print(klc[1].columns)
#print(klc[1].columns.names)

# Set the seed:
rng_seed = 3123194149
random.seed(rng_seed)

# NOTE: pyfits doesn't seem to handle "~/..." etc
datadir = "/Users/pdbaines/Dropbox/GitHub/kepler/data/sandbox/"
cpath = "sandbox1/skygroup27brightvariable/8415109/kplr008415109-2010265121752_llc.fits"
#cpath = "sandbox1/skygroup27TCE/8013439/kplr008013439-2012032013838_slc.fits"
#cpath = "sandbox1/skygroup27examplesystematics/8414982/kplr008414982-2010355172524_llc.fits"
filepath = datadir + cpath

# Find all fits files in the data directory (recursively):
do_search = False
if do_search:
    full_fits_list = []
    for r,d,f in os.walk(datadir):
        for files in f:
            if files.endswith(".fits"):
                full_fits_list.append(os.path.join(r,files))
        print "Finished searching: " + str(r)

# Read in all files and re-save them as .txt:
# TODO

# Extract SAP_FLUX time series:

new_filename = rename_fits_to_txt(filepath)
if os.path.isfile(new_filename):
    print "Extracting pre-converted data from file: " + new_filename
    kepler_ts = np.loadtxt(new_filename)
else:
    print "Extracting data from file: " + filepath
    kepler_ts = get_kepler_data(file_name=filepath,field_name='SAP_FLUX')
    np.savetxt(fname=new_filename,X=kepler_ts,fmt='%10.6f')

kepler_ts = kepler_ts.astype('float') # VERY IMPORTANT!
original_n = len(kepler_ts)

np.sum(np.isnan(kepler_ts)) # contains NA's
kepler_ts = wv.fill_missing(kepler_ts) # also removes stretches of NA's
np.sum(np.isnan(kepler_ts)) # Not any more :)

kepler_ts = wv.truncate_to_power_of_two(kepler_ts)
truncated_n = len(kepler_ts)
if (truncated_n>0):
    print "Threw away " + str(original_n-truncated_n) + " observations for power-of-2 truncation..,"
    print "(corresponds to " + str(round(100.0*float(original_n-truncated_n)/float(original_n),4)) + "% of the series)"

plt.plot(kepler_ts)
plt.show()

# Select type of wavelet:
wavelet = 'haar'
mode = 'sym'

# Wavelet transform:
m_ret = wv.wvanalyze(kepler_ts,wtype='mother',wavelet=wavelet,mode=mode,cmap=cm.afmhot)

# Set certain levels to zero
# Only remove low-frequency patterns,
# medium-frequency (transit-like + other), and,
# high-frequency (outliers + high-frequency instrumental noise)
# remain. Denoised series has mean zero.
m_flatten = copy.deepcopy(m_ret) # must copy, can't do: m_flatten = m_ret!!!
for i in range(len(m_ret)-6,len(m_ret)):
    m_flatten[i][:] = 0.0

# Reconstruct:
kepler_processed = wv.wvreconstruct(m_flatten,wavelet=wavelet,mode=mode)
kepler_denoised = kepler_ts - kepler_processed

# Plot both:
plt.figure(3)
plt.plot(kepler_denoised)
plt.plot(kepler_processed-np.median(kepler_processed))
plt.show()

# Zoomed:
#plot_subset = range(4400,4600)
#plot_subset = range(1000,1500)
#plt.figure(4)
#plt.plot(kepler_denoised[plot_subset])
#plt.plot((kepler_processed-np.median(kepler_processed))[plot_subset])
#plt.show()



