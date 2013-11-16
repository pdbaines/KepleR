import numpy as np
import matplotlib.pyplot as plt
from pylab import *

#To be run at a later time, after the saving.
true_periods = np.load("../Results/true_periods.npy")
true_signals = np.load("../Results/true_signals.npy")
true_noises = np.load("../Results/true_noises.npy")
true_tds = np.load("../Results/true_tds.npy")
true_t0s = np.load("../Results/true_t0s.npy")
true_rhos = np.load("../Results/true_rhos.npy")
true_sigma2s = np.load("../Results/true_sigma2s.npy")

bls_periods = np.load("../Results/bls_periods.npy")


diffs = np.subtract(bls_periods, true_periods)
diffs_rel = 100 * np.divide(diffs, np.abs(true_periods))
SNR = np.divide(signals, np.sqrt(noises))

plt.close('all')
subplot(221)
plt.scatter(SNR, diffs)
plt.ylabel("Estimated Period-True Period")
plt.xlabel("Signal to Noise Ratio")

subplot(222)
plt.scatter(SNR[SNR > 2], diffs[SNR > 2])
plt.ylabel("Estimated Period-True Period")
plt.xlabel("Signal to Noise Ratio > 2")
plt.xlim([2.0, 8.0])

ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
plt.scatter(SNR[SNR < 2], diffs[SNR < 2])
plt.ylabel("Estimated Period-True Period")
plt.xlabel("Signal to Noise Ratio < 2")
plt.xlim([0.0, 2.0])
plt.tight_layout()
plt.savefig('../Results/SNR_vs_AbsDiff.pdf')


#Relative Differences now
plt.close('all')
subplot(221)
plt.scatter(SNR, diffs_rel)
plt.ylabel("Relative Error in Period (%)")
plt.xlabel("Signal to Noise Ratio")

subplot(222)
plt.scatter(SNR[SNR > 2], diffs_rel[SNR > 2])
plt.ylabel("Relative Error in Period (%)")
plt.xlabel("Signal to Noise Ratio > 2")
plt.xlim([2.0, 8.0])
#plt.ylim([-.1, .1])

ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
plt.scatter(SNR[SNR < 2], diffs_rel[SNR < 2])
plt.ylabel("Relative Error in Period (%)")
plt.xlabel("Signal to Noise Ratio < 2")
plt.xlim([0.0, 2.0])
plt.tight_layout()
plt.savefig('../Results/SNR_vs_RelDiff.pdf')