import numpy as np
import matplotlib.pyplot as plt
from pylab import *

#To be run at a later time, after the saving.
true_period = np.load("../Results/true_period.npy")
bls_period = np.load("../Results/bls_period.npy")

true_signal = np.load("../Results/true_signal.npy")
bls_signal = np.load("../Results/bls_depth.npy")

true_td = np.load("../Results/true_td.npy")
bls_q = np.load("../Results/bls_q.npy")
bls_td = np.multiply(bls_q, bls_period)

true_t0 = np.load("../Results/true_t0.npy")
bls_in1 = np.load("../Results/bls_in1.npy")
bls_t0 = bls_period * (bls_in1 / 1500)

true_snr = np.load("../Results/true_snr.npy")


def plot_var(variable, snr=true_snr):
    abs_diff = np.subtract(eval("bls_" + variable), eval("true_" + variable))
    rel_diff = 100 * np.divide(abs_diff, np.abs(eval("true_" + variable)))
    plt.close('all')
    subplot(221)
    plt.scatter(snr, abs_diff)
    plt.ylabel("Discrepancy in " + variable.title())
    plt.xlabel("Signal to Noise Ratio")

    subplot(222)
    plt.scatter(snr[snr > 2], abs_diff[snr > 2])
    plt.ylabel("Discrepancy in " + variable.title())
    plt.xlabel("Signal to Noise Ratio > 2")
    plt.xlim([2.0, 8.0])

    ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
    plt.scatter(snr[snr < 2], abs_diff[snr < 2])
    plt.ylabel("Discrepancy in " + variable.title())
    plt.xlabel("Signal to Noise Ratio < 2")
    plt.xlim([0.0, 2.0])
    plt.tight_layout()
    plt.savefig('../Results/Graphs/SNR_vs_AbsDiff' + variable.title() + '.pdf')

    #Relative Differences now
    plt.close('all')
    subplot(221)
    plt.scatter(snr, rel_diff)
    plt.ylabel("Relative Error in " + variable.title() + " (%)")
    plt.xlabel("Signal to Noise Ratio")

    subplot(222)
    plt.scatter(snr[snr > 3], rel_diff[snr > 3])
    plt.ylabel("Relative Error in " + variable.title() + " (%)")
    plt.xlabel("Signal to Noise Ratio > 2")
    plt.xlim([3.0, 8.0])

    ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
    plt.scatter(snr[snr < 3], rel_diff[snr < 3])
    plt.ylabel("Relative Error in " + variable.title() + " (%)")
    plt.xlabel("Signal to Noise Ratio < 3")
    plt.xlim([0.0, 3.0])
    plt.tight_layout()
    plt.savefig('../Results/Graphs/SNR_vs_RelDiff' + variable.title() + '.pdf')
    return 0

plot_var("period")
plot_var("signal")
plot_var("td")
plot_var("t0")