import numpy as np
import matplotlib.pyplot as plt
from pylab import *

#To be run at a later time, after the saving.
true_period = np.load("../Results/true_period.npy")
bls_period = np.load("../Results/bls_period.npy")

true_signal = np.load("../Results/true_signal.npy")
bls_signal = np.load("../Results/bls_depth.npy")

true_td = np.load("../Results/true_td.npy")
true_q = true_td/true_period

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

plt.close('all')
subplot(311)
plt.scatter(true_period[true_snr>2], bls_period[true_snr>2])
plt.ylabel("Estimated Period")
plt.xlabel("True Period (SNR > 2)")
plt.xlim([8000.0, 20000.0])
plt.ylim([8000.0, 20000.0])

subplot(312)
plt.scatter(true_signal[true_snr>2], bls_signal[true_snr>2])
plt.ylabel("Estimated Signal")
plt.xlabel("True Signal (SNR > 2)")
plt.xlim([0.025, 0.085])
plt.ylim([0.025, 0.085])


subplot(313)
plt.scatter(true_td[true_snr>2], bls_td[true_snr>2])
plt.ylabel("Estimated Transit Duration")
plt.xlabel("True Transit Duration (SNR > 2)")
plt.xlim([40.0, 100.0])
plt.ylim([40.0, 100.0])
plt.tight_layout()
#plt.show()
plt.savefig('../Results/Graphs/True_vs_Est1.pdf', dpi=gcf().dpi)

plt.close('all')
subplot(211)
plt.scatter(true_t0[true_snr>2], bls_t0[true_snr>2])
plt.ylabel("Estimated Time to First Transit")
plt.xlabel("True Time to First Transit (SNR > 2)")
plt.xlim([0.0, 20000.0])
plt.ylim([0.0, 20000.0])

subplot(212)
plt.scatter(true_q[true_snr>2], bls_q[true_snr>2])
plt.ylabel("Estimated Fractional Transit Duration")
plt.xlabel("True Fractional Transit Duration (SNR > 2)")
plt.xlim([0.002, 0.01])
plt.ylim([0.002, 0.01])

plt.tight_layout()
#plt.show()
plt.savefig('../Results/Graphs/True_vs_Est2.pdf', dpi=gcf().dpi)