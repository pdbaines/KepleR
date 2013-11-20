import bls
import math
import glob
import numpy as np

def run_bls(file, nf=10000, a_logP = math.log(8000), b_logP = math.log(20000),
            nb=1500, qmi=0.0005, qma=0.005, verbose=True):
    if verbose:
        print "Reading data..."

    flux = np.loadtxt(file, delimiter=",")

    if verbose:
        print "Setting up arrays..."

    time = range(len(flux))
    n = len(flux)
    u = np.empty(n)
    v = np.empty(n)

    if verbose:
        print "Setting up BLS specs..."

    #Most of this are set through the function call, so I commented this out. For use later, perhaps.
    # Prior for period:
    min_period = math.exp(a_logP)
    max_period = math.exp(b_logP)

    ### Email Specs ####
    fmin = 1 / max_period
    fmax = 1 / min_period
    df = (fmax - fmin) / nf

    if verbose:
        print "EEBLS parameters:"
        print "nf=" + str(nf)
        print "fmin=" + str(fmin)
        print "fmax=" + str(fmax)
        print "df=" + str(df)
        print "nb=" + str(nb)
        print "qmi=" + str(qmi)
        print "qma=" + str(qma)
        print "Running EEBLS..."

    results = bls.eebls(time, flux, u, v, nf, fmin, df, nb, qmi, qma)
    #f = fmin + (np.arange(len(results[0])))*df

    if verbose:
        print "Done. Results:"
        print results

    return results


def validate_bls(data_files=sorted(glob.glob("../Data/y_*.csv")), save=True, times=None):
    if times is None:
        k = len(data_files)  # If we have a times of nothing, process the whole dataset.
    else:
        k = times  # Otherwise, set it to whatever we wanted.
    done = 0  # Number of datasets analyzed so far.

    #Initialize estimated and true arrays.
    bls_period = np.zeros(k)
    bls_power = np.zeros(k)
    bls_depth = np.zeros(k)
    bls_q = np.zeros(k)
    bls_in1 = np.zeros(k)
    bls_in2 = np.zeros(k)
    true_period = np.zeros(k)
    true_signal = np.zeros(k)
    true_snr = np.zeros(k)
    true_td = np.zeros(k)
    true_t0 = np.zeros(k)
    true_sigma2 = np.zeros(k)
    true_rho = np.zeros(k)

    #Using an integer iterator with in-place numpy arrays drastically increases speed, at the cost of simplicity.
    for i, data in enumerate(data_files):
        datanum = data.split('_')[1].split('.')[0]  # Determine which dataset we're on
        best_period, best_power, depth, q, in1, in2 = run_bls(file=data, verbose=False)[1:]
        #Collect BLS estimates for later graphing.
        bls_period[i] = best_period
        bls_power[i] = best_power
        bls_depth[i] = depth
        bls_q[i] = q
        bls_in1[i] = in1
        bls_in2[i] = in2

        #Import and collect true values for later graphing.
        signal, period, td, t0, rho, sigma2 = np.loadtxt("../Data/pars_"+str(datanum)+".csv", delimiter=",")
        true_signal[i] = signal  # Append true signal to array
        true_period[i] = period  # Collect true periods in on array
        true_rho[i] = rho
        true_sigma2[i] = sigma2
        true_snr[i] = signal / (math.sqrt(sigma2/(1-math.pow(rho, 2))))
        true_td[i] = td
        true_t0[i] = t0

        #t0_guess = 8000 + (6*in1/5)
        #print "Index of First Transit: " + str(in1)
        #print "t0 guess: " + str(t0_guess)
        #print "True t0: " + str(t0)
        #print "Relative Error (%): " + str(100*(t0 - t0_guess)/t0)
        #print "dataset:" + datanum
        #print "SNR: " + str(true_snr[i])

        done += 1
        print "Finished dataset " + str(done)

        if done == times:
            break

    if save:
        #Save estimate and true arrays to avoid doing this computation every time.
        np.save("../Results/true_period.npy", true_period)
        np.save("../Results/true_signal.npy", true_signal)
        np.save("../Results/true_td.npy", true_td)
        np.save("../Results/true_t0.npy", true_t0)
        np.save("../Results/true_rho.npy", true_rho)
        np.save("../Results/true_sigma2.npy", true_sigma2)
        np.save("../Results/true_snr.npy", true_snr)

        np.save("../Results/bls_period.npy", bls_period)
        np.save("../Results/bls_power.npy", bls_power)
        np.save("../Results/bls_depth.npy", bls_depth)
        np.save("../Results/bls_q.npy", bls_q)
        np.save("../Results/bls_in1.npy", bls_in1)
        np.save("../Results/bls_in2.npy", bls_in2)

    return None

validate_bls(times=10)
