
import matplotlib.pyplot as plt
import wavelet_funcs as wv
import numpy as np
import math
import pywt

###
from scipy import sparse
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
###

test = True

##################################################
#
# Sampling multivariate normals:
#
# numpy.random.multivariate_normal(mean, cov[, size])
#
# mean : 1-D array_like, of length N
# cov  : 2-D array_like, of shape (N, N)
# size : tuple of ints, optional
#
# Given a shape of, for example, (m,n,k), m*n*k samples
# are generated, and packed in an m-by-n-by-k arrangement.
# Because each sample is N-dimensional, the output shape
# is (m,n,k,N). If no shape is specified, a single (N-D)
# sample is returned.
#
# Each entry out[i,j,...,:] is an N-dimensional value
# drawn from the distribution.
#

def rnorm(mu,sigma,size=1):
    return np.random.normal(loc=mu,scale=sigma,size=size)

def rmvnorm(mu,Sigma,size=1):
    return np.random.multivariate_normal(mean=mu,cov=Sigma,size=size)

##################################################

# (TO BE REMOVED LATER)
if test:
    size = [1000]
    mu = np.array([0.0,1.0,2.0])
    Sigma = np.array([[math.pow(1.0,0),0.0,0.0],[0.0,math.pow(2.0,2),0.0],[0,0,math.pow(3.0,2)]])
    foo = rmvnorm(mu=mu,Sigma=Sigma,size=size)
    print "Multivariate normal samples:"
    print foo
    # Check stats:
    print "Column means = " + str(np.mean(foo,axis=0)) # column means
    print "Column SDs   = " + str(np.std(foo,axis=0))  # column sd's

##################################################
#
# Evaluating gaussian densities:
#
#
####
#
# Sampling inverse-chi^{2}:
#
# R code:
#
#"rinvchisq" <- function(n=1, nu, scale)
#{
#  return(nu*scale/rchisq(n,df=nu))
#}
# Python:

def rinvchisq(nu, scale, size=[1]):
    return (nu*scale)/np.random.chisquare(df=nu, size=size)

##################################################

# (TO BE REMOVED LATER)
if test:
    nu = 12.0
    scale = 1.0
    size = [10000]
    foo = rinvchisq(nu=nu,scale=scale,size=size)
    print foo
    print "Mean   = " + str(np.mean(foo)) + " (should be roughly = " + str(nu*scale/(nu-2)) + ")"
    print "Median = " + str(np.median(foo)) + " (should be lower than " + str(nu*scale/(nu-2)) + ")"

##################################################


def error_check(val,name):
    if name=='psi':
        if np.sum(val<0.0):
            raise Exception, "'psi' must be positive"
    if name=='omega':
        if (val['t_0'] < 0.0):
            Exception, "'t_0' cannot be negative"
        if (val['t_d'] < 0.0):
            Exception, "'t_d' cannot be negative"
        if (val['P'] < 0.0):
            Exception, "'P' cannot be negative"
        if (val['alpha']<0.0) or (val['alpha']>=1.0):
            Exception, "'alpha' must be in [0.0,1.0)"
    if name=='Sigma_w':
        for j in range(0,len(val)):
            if (val[j]<0.0):
                Exception, "'Sigma_w' cannot contain negative entries"
    if name=='d':
        None # Nothing to check...


def eval_q(omega,omega_pars):
    # return q(t|psi)...
    n = omega_pars['n']
    q = np.zeros(n)
    for i in range(0,n):
        mod_time = i % omega['P']
        if (mod_time > omega['t_0']) and (mod_time < (omega['t_0']+omega['t_d'])):
            q[i] = omega['alpha']
    return q

def f_to_w(f,wavelet,mode):
    # wavelet transform from f to w:
    w = pywt.wavedec(f,wavelet=wavelet,mode=mode)
    return w

def w_to_vector(w):
    out = np.array([])
    for i in range(0,len(w)):
        out = np.concatenate([out,w[i]])
    return out

def w_to_array(w,scale_indices):
    out = []
    for i in range(0,len(scale_indices)):
        out.append(w[scale_indices[i]])
    return out

def w_to_f(w,wavelet,mode):
    # wavelet transform from w to f,
    # (requires array-like form for w):
    f = pywt.waverec(w,wavelet=wavelet,mode=mode)
    return f

def find_scale_indices(f,wavelet,mode):
    # Find scale indices...
    w = f_to_w(f=f,wavelet=wavelet,mode=mode)
    # J is defined to be number of wavelet levels:
    J = len(w)
    scale_sizes = np.empty(J).astype('int')
    scale_indices = []
    so_far = 0
    for i in range(0,J):
        n_i = len(w[i])
        scale_sizes[i] = n_i
        scale_indices.append(range(so_far,so_far+n_i))
        so_far += n_i
    return {'scale_indices': scale_indices, 'scale_sizes': scale_sizes}

def omega_to_vector(omega):
    ret = np.empty(len(omega))
    j = 0
    for i in omega:
        ret[j] = omega[i]
        j += 1
    return ret

def K_func(p,srp):
    return (srp * 26 * (p/(48 * 365.25)))

def check_and_reformat_breaks(x,n,name):
    # Error check (var,qtr)_changes:
    B = len(x)+1
    if (B>1):
        x = x.astype('int')
        if np.any(np.logical_or((x<=0),(x>(n-1)))):
            raise Exception, "Elements of '" + str(name) + "' must be between 0 and n-1"
        for i in range(0,len(x)-1):
            if (x[i+1] <= x[i]):
                raise Exception, "'" + str(name) + "' must be an increasing vector of indices"
        # Prepend zero as opening value for first period:
        x = np.append(0,x)
    else:
        x = np.array([0]).astype('int')
    # In all cases, append last value as ending period:
    x = np.append(x,n)
    return x


def sample(n, priors, \
        wavelet='haar', \
        mode='zpd', \
        fixed_pars=[], \
        fixed_par_values=[], \
        var_changes=[], \
        qtr_changes=[], \
        verbose=False):

    ################################################################
    #
    # Inputs:
    #
    # y      :: Light curve as 1D np.array
    # priors :: List of priors for each of psi, omega and Sigma_w
    #
    #           qtr_changes is a vector of indices denoting the start of
    #           each new quarter. This should not include zero, which
    #           is automatically added anyway.
    #
    #           var_changes is a vector of indices denoting variance
    #           changes in the noise process. psi is a vector of
    #           length Q=len(var_changes), for the moment with the same
    #           prior on each element. var_changes could be selected
    #           to match qtr_changes.
    #
    #################################################################
    #
    #           For psi (=sigma^{2}_{epsilon,1},...,sigma^{2}_{epsilon,Q}):
    #
    #           psi_prior = {'nu': (1.0,1.0,...,1.0),
    #                        'scale': (1.0,1.0,...,1.0)}
    #
    #################################################################
    #
    #           For omega (=[t_0,t_d,alpha,P]):
    #
    #           Let, srp = (stellar ratio )^{-1/3}
    #
    #           Priors are:
    #
    #             P   ~ Uniform[ p_{lo} , p_{hi} ]
    #           t_{d} ~ K(P,srp)*(sqrt(1-U[0,1]^{2}))
    #           t_{0} ~ Uniform[ 0.0, P-t_{d} ]
    #           alpha ~ ( Uniform[ 0.5 , 2.0 ]/110.0 )^{2}
    #
    #           so prior parameters are just p_{lo}, p_{hi} and srp.
    #           For example,
    #
    #           omega_prior = {
    #                'P': {'lo': 8000, 'hi': 37000}, \
    #                't_d': {'srp': 1.0}}
    #
    #################################################################
    #
    #           For Sigma_w (=[sigma^{2}_{1},...,sigma^{2}_{J}]):
    #           [[nu_{w,1,0},sigma^{2}_{w,1,0}],
    #            ...
    #            [nu_{w,J,0},sigma^{2}_{w,J,0}]]
    #
    #################################################################
    #
    #           For Sigma_d:
    #           [(1/sigma^{2}_{d,1}),...,(1/sigma^{2}_{d,J})]
    #           Note: precisions not variances to allow for zeros
    #
    #################################################################
    #
    # fixed_pars :: List of any parameters to fix (& values...???)
    #
    ################################################################

    # Useful preliminary stuff:
    n = int(n)
    # Check for power-of-two length (maybe removal later...)
    if math.fabs(math.log(n,2)-int(math.log(n,2)))>1e-10:
        raise Exception, "Non-power of two length for 'y': not implemented yet...!"

    # Indices for the different scaling levels:
    scale_stuff = find_scale_indices(np.zeros(n),wavelet,mode)
    scale_indices = scale_stuff['scale_indices']
    scale_sizes = scale_stuff['scale_sizes']
    J = len(scale_sizes)
    Q = len(var_changes)+1
    R = len(qtr_changes)+1

    # Error check qtr_changes and var_changes:
    var_changes = check_and_reformat_breaks(x=var_changes,n=n,name="var_changes")
    qtr_changes = check_and_reformat_breaks(x=qtr_changes,n=n,name="qtr_changes")

    # Check that var_changes is of length = Q+1:
    if len(var_changes) != (Q+1):
        raise Exception, "Internal error: 'var_changes' of incorrect length"

    # Check that qtr_changes is of length = Q+1:
    if len(qtr_changes) != (R+1):
        raise Exception, "Internal error: 'qtr_changes' of incorrect length"

    # Check that the prior on psi is correctly specified:
    if (len(priors['psi']['nu']) != Q):
        raise Exception, "Prior for psi (nu) of incorrect length"
    if (len(priors['psi']['scale']) != Q):
        raise Exception, "Prior for psi (scale) of incorrect length"

    if verbose and False:
        print "scale sizes:"
        print scale_sizes
        print "scale indices:"
        print scale_indices

    if verbose:
        print "priors:"
        print priors

    # Do omega (transit-related parameters):
    if 'omega' in fixed_pars:
        omega = fixed_par_values['omega']
        error_check(omega,'omega')
    else:
        # Sample...
        P = math.exp(np.random.uniform(low=math.log(priors['omega']['P']['lo']),high=math.log(priors['omega']['P']['hi'])))
        t_d = K_func(p=P,srp=priors['omega']['t_d']['srp'])*math.sqrt(1.0-math.pow(np.random.uniform(),2))
        t_0 = np.random.uniform(low=0.0,high=P-t_d)
        alpha = math.pow(np.random.uniform(low=0.5,high=2.0)/110.0,2.0)
        omega = {'t_0': t_0, 't_d': t_d, 'P': P, 'alpha': alpha}
    if verbose:
        print "Sampled omega: "
        print omega

    # Compute q (transit signature):
    q = eval_q(omega=omega,omega_pars={'n':n})
    if verbose:
        print "Computed q: "
        print q

    # Sample psi (the observation-level (co)variance parameters):
    if 'psi' in fixed_pars:
        psi = fixed_par_values['psi']
        error_check(psi,'psi')
    else:
        psi = rinvchisq(nu=priors['psi']['nu'],scale=priors['psi']['scale'],size=Q)
    if verbose:
        print "Sampled psi: "
        print psi

    # Check prior for d:
    n_d = len(priors['d']['mu'])
    if (priors['d']['Sigma'].shape != (n_d,n_d)):
        raise Exception, "Prior for 'd' of incorrect dimensions"

    # Sample d:
    d = rmvnorm(mu=priors['d']['mu'],Sigma=priors['d']['Sigma'])[0,:]
    if verbose:
        print "Sampled d: "
        print d

    # Sample sigma^{2}_{w,k}:
    if 'Sigma_w' in fixed_pars:
        sigma2_w = fixed_par_values['Sigma_w']
        error_check(sigma2_w,'Sigma_w')
    else:
        # Check prior for Sigma_w:
        if (len(priors['Sigma_w']['nu']) != n_d):
            raise Exception, "Prior for 'Sigma_w' (nu) of incorrect dimensions"
        if (len(priors['Sigma_w']['scale']) != n_d):
            raise Exception, "Prior for 'Sigma_w' (scale) of incorrect dimensions"
        sigma2_w = np.empty(n_d)
        for j in range(0,n_d):
            sigma2_w[j] = rinvchisq(nu=priors['Sigma_w']['nu'][j],scale=priors['Sigma_w']['scale'][j])

    if verbose:
        print "Sampled sigma^{2}_{w}: "
        print sigma2_w

    #
    #############################################################################################
    #

    # Sample w:
    w_vector = np.empty(n)
    so_far = 0
    for k in range(0,len(scale_sizes)):
        for m in range(0,scale_sizes[k]):
            w_vector[so_far] = rnorm(mu=d[k],sigma=math.sqrt(sigma2_w[k]))
            so_far = so_far+1

    if verbose:
        print "Sampled w: "
        print w_vector

    # Transform to f:
    w_array = w_to_array(w=w_vector,scale_indices=scale_indices)
    if verbose:
        print "Array-ified w: "
        print w_array

    f = w_to_f(w=w_array,wavelet=wavelet,mode=mode)
    if verbose:
        print "Converted value of f: "
        print f


    # Compute (1-q)*f:
    mu_y = f-q
    if verbose:
        print "Computed mu_y: "
        print mu_y

    # Sample y:
    y = np.empty(n)
    # Allow for variance changes (recall that len(var_changes)=Q+1):
    for j in range(0,Q):
        for i in range(var_changes[j],var_changes[j+1]):
            y[i] = rnorm(mu_y[i],math.sqrt(psi[j]))
    if verbose:
        print "Sampled y..."

    # Return values:
    theta = {'psi': psi, 'omega': omega, 'Sigma_w': sigma2_w}
    lv = {'d': d, 'f': f, 'w_vector': w_vector, 'w_array': w_array, 'q': q, 'mu_y': mu_y}
    ret = {'y': y, 'theta': theta, 'lv': lv}
    return ret


def posterior_sample(y,priors, \
        nsamples=10,burnin=0,starting_values=[], \
        wavelet='haar',mode='zpd', \
        fixed_pars=[],fixed_par_values=[], \
        store_lvs=False,print_every=10,verbose=False):

    ################################################################
    #
    # Inputs:
    #
    # y      :: Light curve as 1D np.array
    # priors :: List of priors for each of psi, omega and Sigma_w
    #
    #           For psi (=sigma^{2}_{epsilon}):
    #           np.array([nu_{0}, s_{0}^{2}])
    #
    #           For omega (=[t_0,t_d,alpha,P]):
    #           ???
    #
    #           For Sigma_w (=[sigma^{2}_{1},...,sigma^{2}_{J}]):
    #           [[nu_{w,1,0},sigma^{2}_{w,1,0}],
    #            ...
    #            [nu_{w,J,0},sigma^{2}_{w,J,0}]]
    #
    #           For Sigma_d:
    #           [(1/sigma^{2}_{d,1}),...,(1/sigma^{2}_{d,J})]
    #           #NOTE: Changed to variances again!
    #           Note: precisions not variances to allow for zeros
    #
    # fixed_pars :: List of any parameters to fix (& values...???)
    #
    # Parameters:
    #
    # psi     :: Covariance parameters of first level GP
    # omega   :: Parameters of transit
    #            omega['n'] = length of y (needed for computing q(.|omega))
    # Sigma_w :: Variance of wavelet coefficients
    #
    # theta = {'psi': ... , 'omega': ... , 'Sigma_w: ... }
    # lvs = {'f': ... , 'w': ... , 'd': ... , 'Y_mis': ... }
    #
    ################################################################

    # Starting values:

    # Useful preliminary stuff:
    n = len(y)
    # Check for power-of-two length (maybe removal later...)
    if math.fabs(math.log(n,2)-int(math.log(n,2)))>1e-10:
        raise Exception, "Non-power of two length for 'y': not implemented yet...!"

    # Indices for the different scaling levels:
    scale_stuff = find_scale_indices(y,wavelet,mode)
    scale_indices = scale_stuff['scale_indices']
    scale_sizes = scale_stuff['scale_sizes']
    J = len(scale_sizes)

    # Missing data stuff:
    any_missing = np.sum(np.isnan(y))
    if any_missing:
        mis_index = np.isnan(y)
        obs_index = not(mis_index)

    # List of posterior parameters (many are constant across iterations):
    posterior_psi = {'nu': priors['psi']['nu']+n, 'scale': np.nan} # i.e., sigma^{2}_{epsilon}
    posterior_omega = {'n': n}
    posterior_Sigma_w = {'nu': priors['Sigma_w']['nu'] + scale_sizes, 'scale': np.empty(J)}
    posterior_d = {'mu': np.zeros(J), 's2': np.ones(J)}
    posterior_Ymis = {'mu': [], 'Sigma': []}
    posterior_pars = {\
            'psi': posterior_psi, \
            'omega': posterior_omega, \
            'Sigma_w': posterior_Sigma_w, \
            'd': posterior_d, \
            'Ymis': posterior_Ymis}

    # Starting values and fixed quantities:
    if 'psi' in fixed_pars:
        psi_t = fixed_par_values['psi']
        error_check(psi_t,'psi')
    elif 'psi' in starting_values:
        error_check(starting_values['psi'],'psi')
        psi_t = starting_values['psi']
    else:
        # Start at prior median...
        psi_t = priors['psi']['nu']*priors['psi']['scale']/(priors['psi']['nu']+2.0)

    if 'omega' in fixed_pars:
        omega_t = fixed_par_values['omega']
        error_check(omega_t,'omega')
    elif 'omega' in starting_values:
        error_check(starting_values['omega'],'omega')
        omega_t = starting_values['omega']
    else:
        omega_t = {'t_0': 200, 't_d': 40, 'P': 10000.0, 'alpha': 0.01}
        Exception, "Need reasonable starting values for 'omega'"

    # 'transit' modification to signal (could be zero vector)
    q_t = eval_q(omega=omega_t,omega_pars=posterior_pars['omega'])

    if 'Sigma_w' in fixed_pars:
        Sigma_w_t = fixed_par_values['Sigma_w']
        error_check(Sigma_w_t,'Sigma_w')
    elif 'Sigma_w' in starting_values:
        error_check(starting_values['Sigma_w'],'Sigma_w')
        Sigma_w_t = starting_values['Sigma_w']
    else:
        # Start at prior medians...
        Sigma_w_t = np.empty(J)
        for j in range(0,J):
            Sigma_w_t[j] = priors['Sigma_w']['nu'][j]*priors['Sigma_w']['scale'][j]/(priors['Sigma_w']['nu'][j]+2.0)

    # w = W'f, so only one of f and w can be specified...
    if 'w' in fixed_pars:
        w_t = fixed_par_values['w']
        error_check(w_t,'w')
        w_array = w_to_array(w=w_t,scale_indices=scale_indices)
        f_t = w_to_f(w=w_array,wavelet=wavelet,mode=mode)
    elif 'w' in starting_values:
        w_t = starting_values['w']
        error_check(w_t,'w')
        w_array = w_to_array(w=w_t,scale_indices=scale_indices)
        f_t = w_to_f(w=w_array,wavelet=wavelet,mode=mode)
    else:
        # Starting values for f_{t}:
        f_t = y+q_t
        w_array = f_to_w(f=f_t,wavelet=wavelet,mode=mode)
        w_t = w_to_vector(w_array)

    # Compute DWT of (y+q):
    tilde_w_t = w_to_vector(f_to_w(f=y+q_t,wavelet=wavelet,mode=mode))

    if 'd' in fixed_pars:
        d_t = fixed_par_values['d']
        error_check(d_t,'d')
    elif 'd' in starting_values:
        error_check(starting_values['d'],'d')
        d_t = starting_values['d']
    else:
        # Start at mean of wavelet coefficients:
        d_t = np.empty(J)
        for j in range(0,J):
            d_t[j] = np.sum(w_t[scale_indices[j]])

    Ymis_t = []

    # Put it all together:
    theta_t = {'psi': psi_t, 'omega': omega_t, 'Sigma_w': Sigma_w_t}
    lv_t = {'f': f_t, 'w': w_t, 'd': d_t, 'q': q_t, 'tilde_w': tilde_w_t, 'Ymis': Ymis_t}

    if verbose:
        print "Starting values..."
        print "psi = " + str(psi_t)
        print "omega = "
        print omega_t
        print "Sigma_w = "
        print Sigma_w_t
        print "w = "
        print w_t
        print "tilde_w = "
        print tilde_w_t
        print "f = "
        print f_t
        print "q = "
        print q_t
        print "d = "
        print d_t

    # Setup storage for posterior samples:
    psi_draws = np.empty(nsamples)
    omega_draws = np.empty([nsamples,len(theta_t['omega'])])
    Sigma_w_draws = np.empty([nsamples,len(theta_t['Sigma_w'])])
    if store_lvs:
        d_draws = np.empty([nsamples,len(lv_t['d'])])
        w_draws = np.empty([nsamples,len(lv_t['w'])])
        f_draws = np.empty([nsamples,len(lv_t['f'])])
        q_draws = np.empty([nsamples,len(lv_t['f'])])

    # Sample posterior...
    for i in range(0,nsamples+burnin):

        # Sample Ymis_t (if needed):
        if verbose:
            print "Sampling Ymis^{(" + str(i) + ")}..."
        if any_missing:
            # Update posterior parameters:
            posterior_pars['Ymis']['mu'] = 0.0 # TODO
            posterior_pars['Ymis']['Sigma'] = 0.0 # TODO
            # Sample Y_{mis}:
            lv_t['Ymis'] = rmvnorm(\
                    mu=posterior_pars['Ymis']['mu'], \
                    Sigma=posterior_pars['Ymis']['Sigma'])

        # Sample psi_t (i.e., sigma^{2}):
        if verbose:
            print "Sampling psi^{(" + str(i) + ")}..."
        if not('psi' in fixed_pars):
            # Need (tilde_w - d)^{T}(I+\Sigma_{w})^{-1}(tilde_w - d)
            extra_bit = 0.0
            so_far = 0
            for k in range(0,len(scale_sizes)):
                for m in range(0,scale_sizes[k]):
                    extra_bit += math.pow(lv_t['tilde_w'][so_far]-lv_t['d'][k],2)/(1.0+theta_t['Sigma_w'][k])
                    so_far += 1
            posterior_pars['psi']['scale'] = \
                    (priors['psi']['nu']*priors['psi']['scale'] + extra_bit)/posterior_pars['psi']['nu']
            theta_t['psi'] = rinvchisq(nu=posterior_pars['psi']['nu'],scale=posterior_pars['psi']['scale'])
            if (verbose>1):
                print "Sampled psi^{(" + str(i) + ")} = " + str(theta_t['psi']) + " from " + \
                    "rinvchisq(nu=" + str(posterior_pars['psi']['nu']) + ",scale=" + \
                    str(posterior_pars['psi']['scale']) + ") [extra_bit = " + str(extra_bit) + "]"

        # Sample omega_t (i.e., [t_0,t_d,alpha,P]):
        if verbose:
            print "Sampling omega^{(" + str(i) + ")}..."
        if 'omega' not in fixed_pars:
            raise Exception, "Sampling for 'omega' not implemented yet!"

        # Update q_t:
        lv_t['q'] = eval_q(omega=omega_t,omega_pars=posterior_pars['omega'])

        # Compute DWT of (y-q):
        lv_t['tilde_w'] = w_to_vector(f_to_w(f=y+lv_t['q'],wavelet=wavelet,mode=mode))

        # Sample w_t:
        if verbose:
            print "Sampling w^{(" + str(i) + ")}..."
        if 'w' not in fixed_pars:
            so_far = 0
            for k in range(0,len(scale_sizes)):
                tmp_var = theta_t['Sigma_w'][k]/(1.0 + theta_t['Sigma_w'][k])
                tmp_sd = math.sqrt(theta_t['psi']*tmp_var)
                mean_part_one = lv_t['d'][k]/theta_t['Sigma_w'][k]
                for m in range(0,scale_sizes[k]):
                    tmp_mean = (mean_part_one + lv_t['tilde_w'][so_far])*tmp_var
                    lv_t['w'][so_far] = rnorm(tmp_mean,tmp_sd)
                    so_far += 1
                    #if verbose:
                    #    print "Sampling w_{" + str(k) + "," + str(m) + "}^{(" + str(i) + ")} from "
                    #    print "rnorm(mu=" + str(tmp_mean) + ",sd=" + str(tmp_sd) + ")"
            if (verbose>1):
                print "Sampled w^{(" + str(i) + ")} = "
                print lv_t['w']

        # Compute f_{t} for later sampling steps:
        w_array = w_to_array(w=lv_t['w'],scale_indices=scale_indices)
        lv_t['f'] = w_to_f(w=w_array,wavelet=wavelet,mode=mode)

        # Sample d_{t}:
        if verbose:
            print "Sampling d^{(" + str(i) + ")}..."
        if 'd' not in fixed_pars:
            # Compute posterior parameters for d_{t}:
            #TODO: Wrong now -- no longer precisions!
            assert 0
            for j in range(0,J):
                posterior_pars['d']['s2'][j] = 1.0/ \
                        (priors['Sigma_d'][j] + (scale_sizes[j]/(theta_t['Sigma_w'][j]*theta_t['psi']))) ###
                posterior_pars['d']['mu'][j] = posterior_pars['d']['s2'][j]* \
                        (priors['mu_d'][j] + (np.sum(lv_t['w'][scale_indices[j]])/(theta_t['Sigma_w'][j]*theta_t['psi'])))
                lv_t['d'][j] = rnorm(mu=posterior_pars['d']['mu'][j],sigma=math.sqrt(posterior_pars['d']['s2'][j]))
            if (verbose>1):
                for j in range(0,J):
                    print "Sampled d_{" + str(j) + "}^{(" + str(i) + ")} = " + \
                            str(lv_t['d'][j]) + " from " + \
                            "N(" + str(posterior_pars['d']['mu'][j]) + "," + \
                            str(math.sqrt(posterior_pars['d']['s2'][j])) + ")"


        # Sample Sigma_w_t (i.e., [sigma^{2}_{1},...,sigma^{2}_{J}])
        if verbose:
            print "Sampling Sigma_{w}^{(" + str(i) + ")}..."
        if 'Sigma_w' not in fixed_pars:
            for j in range(0,J):
                # d is just a constant within each level:
                w_minus_d = lv_t['w'][scale_indices[j]]-lv_t['d'][j]
                extra_bit = np.dot(w_minus_d,w_minus_d)/theta_t['psi'] ###
                posterior_pars['Sigma_w']['scale'][j] = \
                        (priors['Sigma_w']['nu'][j]*priors['Sigma_w']['scale'][j] + \
                                extra_bit)/posterior_pars['Sigma_w']['nu'][j]
                theta_t['Sigma_w'][j] = \
                        rinvchisq(nu=posterior_pars['Sigma_w']['nu'][j],scale=posterior_pars['Sigma_w']['scale'][j])
                if (verbose>1):
                    print "Sampled Sigma_{w," + str(j) + "}^{(" + str(i) + ")}"
                    print "from rinvchisq(nu=" + str(posterior_pars['Sigma_w']['nu'][j]) + \
                            ",scale=" + str(posterior_pars['Sigma_w']['scale'][j]) + ")"

        # Store the draws (as copies):
        if (i>=burnin):
            psi_draws[i-burnin] = theta_t['psi']
            omega_draws[i-burnin,:] = omega_to_vector(theta_t['omega'])
            Sigma_w_draws[i-burnin,:] = theta_t['Sigma_w']
            if store_lvs:
                d_draws[i-burnin,:] = lv_t['d']
                w_draws[i-burnin,:] = lv_t['w']
                f_draws[i-burnin,:] = lv_t['f']
                q_draws[i-burnin,:] = lv_t['q']

        if (i%print_every == 0):
            print "Finished iteration " + str(i) + "..."

    all_draws = {}
    all_draws['psi'] = psi_draws
    all_draws['omega'] = omega_draws
    all_draws['Sigma_w'] = Sigma_w_draws
    if store_lvs:
        all_draws['d'] = d_draws
        all_draws['w'] = w_draws
        all_draws['f'] = f_draws
        all_draws['q'] = q_draws

    return all_draws


######################################################################

# Example run...

# Limb darkening: 2-4 parameters depending on models
# May be possible to get priors for these for specific
# Prior stuff:
# n = 2^15 = 60,000

##########################
verbose = True
wavelet = 'haar'
mode = 'zpd'
fixed_pars = []
fixed_par_values = []
make_plots = True
##########################

J = 15
n = int(math.pow(2,J))

# Roughly 4320 points per quarter
ppq = 4320
R = int(n/ppq)+1
qtr_changes = np.empty(R-1)
for i in range(0,R-1):
    qtr_changes[i] = ppq*(i+1)

# Just assume variance changes are at quarter-changes too:
var_changes = qtr_changes

# Prior for psi:
psi_nu_prior = np.empty(R)
psi_scale_prior = np.empty(R)
for i in range(0,R):
    psi_nu_prior[i] = 20.0
    psi_scale_prior[i] = 0.01

psi_prior = {'nu': psi_nu_prior, 'scale': psi_scale_prior}

# t_d ~= 10hrs * P^(1/3)
# Period: 0.5yr - 2yrs, (8000 - 37000)
# Duration: ~10hrs, ~20 data points (1hr-30hrs, earth-like: 6.5-13hrs)
# Transit signal reduction: 80 * 10^(-6)

# srp = (stellar ratio)^(-1/3)
omega_prior = { \
        'P': {'lo': 8000, 'hi': 37000}, \
        't_d': {'srp': 1.0}}

nlevels = len(f_to_w(f=np.zeros(n),wavelet=wavelet,mode=mode))

#############################################################################
# TODO:

# w and d priors:

# Length of d:
n_d = nlevels

# Prior place holders:
nu_w = np.zeros(n_d)
scale_w = np.zeros(n_d)
mu_d_prior = np.zeros(n_d)
Sigma_d_prior = np.zeros([n_d,n_d])

# Overall scaling (not used, currently set to a constant for simulation):
nu_w[0] = 100.0
scale_w[0] = 100.0
mu_d_prior[0] = math.pow(2.0,0.5*J)
Sigma_d_prior[0,0] = 1.0

# Mother wavelet priors:
for j in range(1,n_d):
    nu_w[j] = 100.0
    scale_w[j] = 0.01
    Sigma_d_prior[j,j] = 1.0/10000.0

d_prior = {'mu': mu_d_prior, 'Sigma': Sigma_d_prior}
Sigma_w_prior = {'nu': nu_w, 'scale': scale_w}

# Full prior list:
prior_list = { \
        'psi': psi_prior, \
        'omega': omega_prior, \
        'Sigma_w': Sigma_w_prior, \
        'd': d_prior}

###################################################

sim_data = sample(n=n,priors=prior_list, \
        wavelet=wavelet,mode=mode, \
        fixed_pars=fixed_pars, \
        fixed_par_values=fixed_par_values, \
        var_changes = var_changes, \
        verbose=verbose)

###################################################

fig_ctr = 0

# Plots:
if make_plots:
    fig_ctr += 1
    plt.figure(fig_ctr)
    plt.subplot(2,1,1)
    plt.plot(sim_data['y'])
    plt.subplot(2,1,2)
    plt.plot(sim_data['lv']['q'])
    plt.show()
    fig_ctr += 1
    plt.figure(fig_ctr)
    plt.plot(sim_data['y'])
    plt.xlim((sim_data['theta']['omega']['t_0']-30,sim_data['theta']['omega']['t_0']+30))
    plt.show()
    fig_ctr += 1
    plt.figure(fig_ctr)
    plt.show()
    wv.wvplot(y=sim_data['y'],f=sim_data['lv']['w_array'],wtype='mother')

assert 0

# MCMC stuff:
burnin = 100
nsamples = 500
store_lvs = True
verbose = False # 2.0

fixed_par_values = {'omega': { \
        't_0': sim_data['theta']['omega']['t_0'], \
        't_d': sim_data['theta']['omega']['t_d'], \
        'P': sim_data['theta']['omega']['P'], \
        'alpha': sim_data['theta']['omega']['alpha']},
        'psi': sim_data['theta']['psi'],
        'd': sim_data['lv']['d'],
        'w': sim_data['lv']['w_vector'],
        'Sigma_w': sim_data['theta']['Sigma_w']}

#
# psi slightly over-estimated, maybe sigma_d not 1/sigma_d somewhere?
# check stuff more clearly...
#
#                    [d] | else ==> ok
#                  [psi] | else ==> ok (500 iterations, 100 burnin)
#====================================================================
#                    [w] | else ==> ok
#              [Sigma_w] | else ==> ok
#               [psi, w] | else ==> ok (500 iterations, 100 burnin)
#         [psi, Sigma_w] | else ==> ok
#           [w, Sigma_w] | else ==> ok (500 iterations, 100 burnin)
#      [psi, w, Sigma_w] | else ==> ok
#   [d, psi, w, Sigma_w] | else ==> ok! :)
#

fixed_pars = ['omega']
#fixed_pars = ['omega','Sigma_w','w','d']

###################################################

# Sample posterior:
out = posterior_sample(y=sim_data['y'], \
        priors=prior_list, \
        nsamples=nsamples,burnin=burnin, \
        fixed_pars=fixed_pars, \
        fixed_par_values=fixed_par_values, \
        store_lvs=store_lvs, \
        verbose=verbose)

###################################################

if 'psi' not in fixed_pars:
    fig_ctr += 1
    plt.figure(fig_ctr)
    plt.subplot(2,1,1)
    plt.hist(out['psi'][:],color='g')
    plt.axvline(sim_data['theta']['psi'], color='r', linestyle='dashed', linewidth=2)
    plt.subplot(2,1,2)
    plt.plot(out['psi'][:])
    plt.axhline(sim_data['theta']['psi'], color='r', linestyle='dashed', linewidth=2)
    plt.show()

if 'Sigma_w' not in fixed_pars:
    fig_ctr += 1
    plt.figure(fig_ctr)
    for i in range(0,out['Sigma_w'].shape[1]):
        plt.subplot(4,4,1+i)
        plt.hist(out['Sigma_w'][:,i],color='y')
        plt.axvline(sim_data['theta']['Sigma_w'][i], color='r', linestyle='dashed', linewidth=2)
    plt.show()

if store_lvs and ('w' not in fixed_pars):
    fig_ctr += 1
    plt.figure(fig_ctr)
    for i in range(0,16):
        plt.subplot(4,4,1+i)
        plt.hist(out['w'][:,i],color='b')
        plt.axvline(sim_data['lv']['w_vector'][i], color='r', linestyle='dashed', linewidth=2)
    plt.show()

if store_lvs and ('d' not in fixed_pars):
    fig_ctr += 1
    plt.figure(fig_ctr)
    for i in range(0,16):
        plt.subplot(4,4,1+i)
        plt.hist(out['d'][:,i],color='r')
        plt.axvline(sim_data['lv']['d'][i], color='b', linestyle='dashed', linewidth=2)
    plt.show()




