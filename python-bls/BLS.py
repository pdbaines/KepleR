import bls, math, glob
import numpy as np
import csv

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

    if verbose:
        print "Done. Results:"
        print results

    return results

#Get file list, take only simulated observation CSVs
data_files = glob.glob("../Data/y_*.csv")
outfilestub = "../Results/out_"

for data in data_files:
    foo = run_bls(file=data, verbose=True)
    datanum = data.split('_')[1].split('.')[0]
    outfile = outfilestub + datanum + ".txt"
    ofile = open(outfile,"wb")
    ofile.write(str(foo[1]))
    ofile.close()



