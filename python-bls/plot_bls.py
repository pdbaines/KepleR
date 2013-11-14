import numpy as np
import matplotlib.pyplot as plt
from pylab import *

#To be run at a later time, after the saving.
true_periods = np.load("true_periods.npy")
bls_periods = np.load("bls_periods.npy")
signals = np.load("signals.npy")
noises = np.load("noises.npy")

diffs = np.subtract(bls_periods, true_periods)
SNR = np.divide(signals, np.sqrt(noises))

hi_SNR = SNR[SNR > 2]
hi_diff = diffs[SNR > 2]
lo_SNR = SNR[SNR < 2]
lo_diff = diffs[SNR < 2]

diffs_rel = np.divide(diffs, np.abs(true_periods))

plt.close('all')
subplot(221)
plt.scatter(SNR, diffs)
plt.ylabel("Estimated Period-True Period")
plt.xlabel("Signal to Noise Ratio")
plt.title("SNR v. Period Difference")

subplot(222)
plt.scatter(SNR, diffs_rel)
plt.ylabel("Est. Period - True Period / |True|")
plt.xlabel("Signal to Noise Ratio")

subplot(223)
plt.scatter(lo_SNR, lo_diff)
plt.ylabel("Estimated Period-True Period")
plt.xlabel("Signal to Noise Ratio < 2")

subplot(224)
plt.scatter(hi_SNR, hi_diff)
plt.ylabel("Estimated Period-True Period")
plt.xlabel("Signal to Noise Ratio > 2")

plt.tight_layout()

plt.savefig('SNR_vs_Diff.pdf')