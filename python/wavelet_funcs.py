import matplotlib.pyplot as plt
from matplotlib import cm
import pywt
import scipy
import numpy as np
import math # for isnan

def mwavelets(y, wavelet, mode):
    decomp = pywt.wavedec(y, wavelet=wavelet, mode=mode)
    # (Do not) discard first level (scaling) coefficients
    # need for full reconstruction!
    return decomp

def fwavelets_rev(y, wavelet, mode):
    ret = []
    # recursively compute...
    tmp_y = y
    nlevels = math.log(len(y),2)
    if math.fabs(nlevels-int(nlevels))>1e-10:
        print "Non power of two length!"
    nlevels = int(nlevels)
    for i in range(1,nlevels):
        tmp_decomp = pywt.wavedec(y, wavelet=wavelet, mode=mode, level=i)
        tmp_y = tmp_decomp[0]
        ret.append(tmp_y)
    return ret

def fwavelets(y, wavelet, mode):
    ret = []
    # recursively compute...
    tmp_y = y
    nlevels = math.log(len(y),2)
    if math.fabs(nlevels-int(nlevels))>1e-10:
        print "Non power of two length!"
    nlevels = int(nlevels)
    for i in range(1,nlevels):
        tmp_decomp = pywt.wavedec(y, wavelet=wavelet, mode=mode, level=i)
        tmp_y = tmp_decomp[0]
        ret.append(tmp_y)
    return ret


def wvreconstruct(y,wavelet,mode,wtype='mother'):
    if wtype=='mother':
        ret = pywt.waverec(y,wavelet,mode)
    elif wtype=='father':
        print "Reconstruction from father wavelets not implemented yet"
        ret = None
    else:
        print "Invalid 'wtype'"
        ret = None
    return ret


def mtodense(x):
    # First level is overall scale (remove)
    x = x[1:]
    # Now as usual:
    J = len(x)
    N = 2*len(x[J-1])
    out = np.empty([J,N])
    for i in range(0,J):
        pads = N/len(x[i])
        for j in range(0,len(x[i])):
            for k in range(0,pads):
                out[i,(j*pads)+k] = x[i][j]
    return out

def ftodense(x):
    J = len(x)
    N = 2*len(x[0])
    out = np.empty([J,N])
    for i in range(0,J):
        pads = N/len(x[i])
        for j in range(0,len(x[i])):
            for k in range(0,pads):
                out[i,(j*pads)+k] = x[i][j]
    return out

def wvplot(y,f,wtype='father',cmap=cm.coolwarm,origin='lower',interpolation='none'):
    if wtype=='father':
        fgrid = ftodense(f)
    elif wtype=='mother':
        fgrid = mtodense(f)
    else:
        print "Invalid input for 'wtype'"
        return None
    fig1 = plt.figure(1)
    asp = np.shape(fgrid)[1]/np.shape(fgrid)[0]
    plt.imshow(fgrid,aspect=asp,cmap=cmap,origin=origin,interpolation=interpolation)
    plt.colorbar(orientation='vertical')
    fig2 = plt.figure(2)
    plt.plot(y)
    plt.xlim([0,len(y)])
    fig1.show()
    fig2.show()

def wvanalyze(y,wtype,wavelet,mode,cmap=cm.coolwarm,origin='lower',interpolation='none'):
    if wtype=='father':
        wd = fwavelets(y=y,wavelet=wavelet,mode=mode)
    elif wtype=='mother':
        wd = mwavelets(y=y,wavelet=wavelet,mode=mode)
    else:
        print "Invalid input for 'wtype'"
        return None
    wvplot(y=y,f=wd,wtype=wtype,cmap=cmap,origin=origin,interpolation=interpolation)
    return wd

def fill_missing(y,boxwidth=3,min_non_na=2,warn=True):
    z = np.array(y) # copy to avoid long-stretches of NA's
    N = len(y)
    remove_idx = []
    for i in range(0,N):
        if math.isnan(y[i]):
            tmp_sum = 0
            non_na = 0
            for j in range(i-boxwidth,i+boxwidth+1):
                if (j<0) or (j>=N):
                    continue
                elif not(math.isnan(y[j])):
                    tmp_sum += y[j]
                    non_na += 1
            if (non_na>min_non_na):
                z[i] = (tmp_sum/non_na)
            else:
                remove_idx.append(i)
    # Remove any remaining stretches of NA's:
    z = np.delete(z,remove_idx)
    return z

def truncate_to_power_of_two(y):
    N = len(y)
    actual_b2len = math.log(N,2)
    new_N = int(pow(2.0,math.floor(actual_b2len)))
    ntrunc = N-new_N
    ix_list = []
    for i in range(0,ntrunc):
        if i%2:
            tix = (i/2) # 0=>0, 2=>1, 4=>2 etc.
            ix_list.append(tix) # cut from front
        else:
            tix = (N-1)-int(math.floor(i/2)) # 1=>N-1, 3=>N-2, 5=>N-3 etc.
            ix_list.append(tix)
    new_y = np.delete(y,ix_list)
    return new_y


