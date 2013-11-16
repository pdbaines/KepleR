import numpy as np
import matplotlib.pyplot as plt
from pylab import *

#To be run at a later time, after the saving.
true_period = np.load("../Results/true_period.npy")
true_signal = np.load("../Results/true_signal.npy")
true_noise = np.load("../Results/true_noise.npy")
true_td = np.load("../Results/true_td.npy")
true_t0 = np.load("../Results/true_t0.npy")
true_rho = np.load("../Results/true_rho.npy")
true_sigma2 = np.load("../Results/true_sigma2.npy")

bls_period = np.load("../Results/bls_period.npy")


diff = np.subtract(bls_period, true_period)
diff_rel = 100 * np.divide(diff, np.abs(true_period))
SNR = np.divide(true_signal, np.sqrt(true_noise))

plt.close('all')
subplot(221)
plt.scatter(SNR, diff)
plt.ylabel("Estimated Period-True Period")
plt.xlabel("Signal to Noise Ratio")

subplot(222)
plt.scatter(SNR[SNR > 2], diff[SNR > 2])
plt.ylabel("Estimated Period-True Period")
plt.xlabel("Signal to Noise Ratio > 2")
plt.xlim([2.0, 8.0])

ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
plt.scatter(SNR[SNR < 2], diff[SNR < 2])
plt.ylabel("Estimated Period-True Period")
plt.xlabel("Signal to Noise Ratio < 2")
plt.xlim([0.0, 2.0])
plt.tight_layout()
plt.savefig('../Results/SNR_vs_AbsDiff.pdf')


#Relative Differences now
plt.close('all')
subplot(221)
plt.scatter(SNR, diff_rel)
plt.ylabel("Relative Error in Period (%)")
plt.xlabel("Signal to Noise Ratio")

subplot(222)
plt.scatter(SNR[SNR > 2], diff_rel[SNR > 2])
plt.ylabel("Relative Error in Period (%)")
plt.xlabel("Signal to Noise Ratio > 2")
plt.xlim([2.0, 8.0])
#plt.ylim([-.1, .1])

ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
plt.scatter(SNR[SNR < 2], diff_rel[SNR < 2])
plt.ylabel("Relative Error in Period (%)")
plt.xlabel("Signal to Noise Ratio < 2")
plt.xlim([0.0, 2.0])
plt.tight_layout()
plt.savefig('../Results/SNR_vs_RelDiff.pdf')