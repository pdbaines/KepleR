import bls
import math
import glob
import numpy as np
import matplotlib.pyplot as plt
<<<<<<< HEAD
=======
import csv
>>>>>>> b47890c16c749a0376d053259408d186bb2caacd

def run_bls(file='./obs.csv', nf=10000, a_logP = math.log(8000), b_logP = math.log(20000),
            nb=1500, qmi=0.0005, qma=0.005, verbose=True):
    if verbose:
        print "Reading data..."

    flux = np.genfromtxt(file, delimiter=",")

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
    #a_logP = math.log(8000) # 8000
    #b_logP = math.log(20000) # 37000
    min_period = math.exp(a_logP)
    max_period = math.exp(b_logP)

    ### Email Specs ####
    #nf = 10000
    fmin = 1 / max_period
    fmax = 1 / min_period
    df = (fmax - fmin) / nf

    #nb = 1500
    #qmi = 0.0005
    #qma = 0.005

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
    f = fmin + (np.arange(len(results[0])))*df
    
    if verbose:
        print "Done. Results:"
        print results

    return results, 1.0/f

#Get file list, take simulated observation CSVs
data_files = glob.glob("../Data/y_*.csv")
#print data_files
outfilestub = "../Results/out_"

true_periods = np.empty(1)
bls_periods = np.empty(1)
signals = np.empty(1)
noises = np.empty(1)

for data in data_files:  # This one takes a long time. Better to prototype on small sets.
#for data in ['../Data/y_5030037.csv', '../Data/y_5030038.csv', '../Data/y_5030039.csv', '../Data/y_5030040.csv']:
    bls_results = run_bls(file=data, verbose=False)
    datanum = data.split('_')[1].split('.')[0]
    print "Running dataset number " + str(datanum)
    if datanum == "5029888":
      qlk4heqlk
    outfile = outfilestub + datanum + ".txt"
    ofile = open(outfile, "wb")
    ofile.write(str(bls_results[1]))
    ofile.close()

<<<<<<< HEAD
    pars = np.genfromtxt("../Data/pars_"+str(datanum)+".csv", delimiter=",") # Get true parameters
    signal = pars[0]  # Get true signal
    noise = pars[5]/(1-math.pow(pars[4], 2))  # Get true noise, sigma^2 / (1-rho^2)
    true_period = pars[1]  # Get true periods

    bls_periods = np.append(bls_periods, bls_results[1])  # Collect all estimated periods into one array
    true_periods = np.append(true_periods, true_period)  # Collect those, too
    signals = np.append(signals, signal)
    noises = np.append(noises, noise)


true_periods = true_periods[1:]
bls_periods = bls_periods[1:]
signals = signals[1:]
noises = noises[1:]

diffs = true_periods - bls_periods
SNR = signals/noises
=======
plt.scatter(y=foo[0][0],x=foo[1])
plt.plot()
plt.show()
>>>>>>> b47890c16c749a0376d053259408d186bb2caacd

#plt.scatter(diffs, SNR)
#plt.xlabel("Difference in Estimated Period and True Period")
#plt.ylabel("Signal to Noise Ratio")
#plt.title("Signal to Noise Ratio versus Difference in Period")
#plt.savefig('SNR_vs_Diff.pdf')
