import bls
import math
import glob
import numpy as np
import matplotlib.pyplot as plt

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

# hack
todo = 100
done = 0

#for data in ['../Data/y_5030037.csv', '../Data/y_5030038.csv']:
for data in data_files:  # This one takes a long time. Better to prototype on small sets.
    bls_results = run_bls(file=data, verbose=False)
    datanum = data.split('_')[1].split('.')[0]
    outfile = outfilestub + datanum + ".txt"
    ofile = open(outfile, "wb")
    ofile.write(str(bls_results[0][1]))
    ofile.close()

    pars = np.genfromtxt("../Data/pars_"+str(datanum)+".csv", delimiter=",") # Get true parameters
    signal = pars[0]  # Get true signal
    noise = pars[5]/(1-math.pow(pars[4], 2))  # Get true noise, sigma^2 / (1-rho^2)
    #print "SD" + str(pars[5])
    #print "1-rho2" + str((1-math.pow(pars[4], 2)))
    #print noise
    true_period = pars[1]  # Get true periods

    bls_periods = np.append(bls_periods, bls_results[0][1])  # Collect all estimated periods into one array
    true_periods = np.append(true_periods, true_period)  # Collect those, too
    signals = np.append(signals, signal)
    noises = np.append(noises, noise)
    done = done+1
    print "Finished dataset " + str(done)
    if done == todo:
        break

true_periods = true_periods[1:]
bls_periods = bls_periods[1:]
signals = signals[1:]
noises = noises[1:]

diffs = np.subtract(true_periods, bls_periods)
SNR = np.divide(signals, np.sqrt(noises))

print diffs
print SNR

##Save files to avoid this computation later.
np.save("true_periods", true_periods)
np.save("bls_periods", bls_periods)
np.save("signals", signals)
np.save("noises", noises)

plt.scatter(SNR, diffs)
plt.ylabel("Difference in Estimated Period and True Period")
plt.xlabel("Signal to Noise Ratio")
plt.title("Signal to Noise Ratio versus Difference in Period")
plt.savefig('SNR_vs_Diff.pdf')

