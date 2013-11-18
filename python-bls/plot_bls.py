import numpy as np
import matplotlib.pyplot as plt
from pylab import *

#To be run at a later time, after the saving.
true_period = np.load("../Results/true_period.npy")
true_signal = np.load("../Results/true_signal.npy")
true_noise = np.load("../Results/true_noise.npy")
true_td = np.load("../Results/true_td.npy")
true_t0 = np.load("../Results/true_t0.npy")

true_sigma2 = np.load("../Results/true_sigma2.npy")
true_rho = np.load("../Results/true_t0.npy")

bls_period = np.load("../Results/bls_period.npy")
bls_power = np.load("../Results/bls_power.npy")
bls_depth = np.load("../Results/bls_depth.npy")
bls_q = np.load("../Results/bls_q.npy")
bls_in1 = np.load("../Results/bls_in1.npy")
bls_in2 = np.load("../Results/bls_in2.npy")

diff_period = np.subtract(bls_period, true_period)
diff_rel_period = 100 * np.divide(diff_period, np.abs(true_period))
#SNR = np.divide(true_signal, np.sqrt(true_noise))


def plot_var(variable, SNR=np.divide(true_signal, np.sqrt(true_noise))):
    diff = np.subtract(eval("bls_" + variable), eval("true_" + variable))
    rel_diff = 100 * np.divide(diff, np.abs(eval("true_" + variable)))

    plt.close('all')
    subplot(221)
    plt.scatter(SNR, diff)
    plt.ylabel("Discrepancy in " + variable.title())
    plt.xlabel("Signal to Noise Ratio")

    subplot(222)
    plt.scatter(SNR[SNR > 2], diff[SNR > 2])
    plt.ylabel("Discrepancy in " + variable.title())
    plt.xlabel("Signal to Noise Ratio > 2")
    plt.xlim([2.0, 8.0])

    ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
    plt.scatter(SNR[SNR < 2], diff[SNR < 2])
    plt.ylabel("Discrepancy in " + variable.title())
    plt.xlabel("Signal to Noise Ratio < 2")
    plt.xlim([0.0, 2.0])
    plt.tight_layout()
    plt.savefig('../Results/SNR_vs_AbsDiff' + variable.title() + '.pdf')

    #Relative Differences now
    plt.close('all')
    subplot(221)
    plt.scatter(SNR, rel_diff)
    plt.ylabel("Relative Error in " + variable.title() + " (%)")
    plt.xlabel("Signal to Noise Ratio")

    subplot(222)
    plt.scatter(SNR[SNR > 2], rel_diff[SNR > 2])
    plt.ylabel("Relative Error in " + variable.title() + " (%)")
    plt.xlabel("Signal to Noise Ratio > 2")
    plt.xlim([2.0, 8.0])

    ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
    plt.scatter(SNR[SNR < 2], rel_diff[SNR < 2])
    plt.ylabel("Relative Error in " + variable.title() + " (%)")
    plt.xlabel("Signal to Noise Ratio < 2")
    plt.xlim([0.0, 2.0])
    plt.tight_layout()
    plt.savefig('../Results/SNR_vs_RelDiff' + variable.title() + '.pdf')
    return 0

plot_var("period")