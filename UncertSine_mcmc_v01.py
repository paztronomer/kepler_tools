# Script to estimate uncertainties in sine fit to a Kepler light curve
# If use/modify/distribute, refer to: Francisco Paz-Chinchon,
#                                     francisco at dfte.ufrn.br
#                                     DFTE, UFRN, Brazil.


import emcee
from pandas import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq     #least squares optimization
import scipy.optimize as op  # optimization
from sys import exit
import time 
import os
import gc

#
#           F u n c t i  o n s
# 

#           Add prefix to kic 
def Kname(idk):
    idk = str(idk)
    if len(idk) < 9:
        control = True
        aux_id = idk
        while control:
            if len(aux_id) == 9: control = False
            else: aux_id = '0' + aux_id
        return 'kplr' + aux_id 
    else: 
        return 'kplr' + idk 


#           Read LC and estimate sigma_{LC-ppalSignal}
def Read(gname, vargs):
    prot, ampl, nsLev = vargs # expand 
    # read lc
    TS = read_csv(gname, sep=' ', names=['time','nflux'])
    tjd, pdc = TS['time'].values, TS['nflux'].values
    TS = None    
    # estimate sigma_{lc-ppalSignal}
    # 1st) optimize parameters
    guess_a, guess_b = 0.0, np.mean(pdc)
    optimize_func = lambda x: ampl*np.sin(tjd*2.0*np.pi*(1./prot)+x[0]) + x[1] - pdc
    est_a, est_b = leastsq(optimize_func, [guess_a, guess_b])[0]
    phi, off = est_a, est_b
    est_a, est_b, guess_a, guess_b = None, None, None, None
    # 2nd) estimate sigma
    dd = pdc - (ampl*np.sin(tjd*2.0*np.pi*(1./prot)+phi) + off); SS = np.std(dd)
    dd = None
    # 3rd) add the noise contribution
    SS = SS + nsLev
    return tjd, pdc, SS, phi, off


#.......................................................

#
#           C o r p u s
#

print '\n\tstarting!....'

global_t = time.time()

# Path to light curves
# path_lc = '/var/run/media/fj/Data/tesis.phd/2014/KOI-Dec13/pdc_s05_UNtreat/'
path_lc = '/Users/fj/Work/tesis.phd/2014/ConfirPlanet-May20/pdcM20/s01tos04_treat/'

# define number of walkers and number of iterations
Nit = 1000
Wks = 200 
THREADS = 1

#           L o a d  tables
#
#           1) Table of uncertainties
# 01_KIC,02_pdc_err,03_normpdc_err,04_time_corr
# this table has the error in flux for sample s04 + sample KOI Dec2013... so all stars that we
# possible need
e_tab = read_csv('s04d13.r05_ErrorFITS.csv')
#
e_tab['05_kic_number'] = map(lambda x: int(x[ x.find('kplr')+4: ]), e_tab['01_KIC'].values)


#           2) Table of period and other info
# 01_kic,02_period,03_amplitude,04_power,05_peak,06_noiseLev,07_snr
# info_tab = read_csv('final_KOInonPlan_all_v03.csv')
# 
# 01_kic,02_period,03_amplitude,04_power,05_peak,06_noiseLev,07_snr,08_Msun,09_Rsun,10_teff
info_tab = read_csv('final_Plan_fromKCPl_v02.csv')

#           List to harbor acceptance fraction
# Note: the acceptance fraction is like a diagnosis, which must be have values between
#          0.2 and 0.5
af = [ [],[],[] ]

#           List to harbor varios parameters
comp=[[],[],[],[],[], [],[],[],[],[], [],[],[],[]]

#           Walk through ids
for ind_kic,kic in enumerate( info_tab['01_kic'].values ):
    start_time = time.time()
    # free memory 
    gc.collect()
  
    # get values from tables for the actual KIC 
    period = info_tab['02_period'].values[ind_kic]
    ampl = info_tab['03_amplitude'].values[ind_kic]  
    nsLev = info_tab['06_noiseLev'].values[ind_kic]
    e_flux = e_tab.loc[ (e_tab['05_kic_number'] == kic), '03_normpdc_err'].values[0]
    corr_time = e_tab.loc[ (e_tab['05_kic_number'] == kic), '04_time_corr'].values[0]
  
    # Walk through files up to match light curve id
    for (path, dirs, files) in os.walk( path_lc ):
        for ind_FF, FF in enumerate(files):    
            kID = map(Kname,  [kic])[0]
            if ('.dat' in FF) and ( kID in FF ): 
                print '\n.........................'
                print '\t{0:^12}'.format(kID)
                
                # returns: time_arr, flux_arr, sigma_total, phase, offset
                tjd, flux, ss_total, ph, delta = Read( path_lc+FF, [period,ampl,nsLev])
                
                #           M a x i m u m  l i k e l i h o o d
                # making the chi-squared statistics as small as possible
                # optimize.minimize(function, x0, args=(extra_arguments_to_function), ...)
                # chi2 returns -2.*lnlike
                #           Likelihood function
                def lnlike(vargs, x, y, yerr):
                    A, freq, lnf = vargs
                    #phi, offset = rho
                    model = A * np.sin(2.0*np.pi*freq*x + ph) + delta
                    inv_sigma2 = 1.0/(yerr**2 + np.exp(2*lnf))
                    # note that exp(2 ln f)==f**2
                    return -0.5*(np.sum( (y-model)**2*inv_sigma2 - np.log(inv_sigma2) ))


                print 'fix it!!!! : np.log(ss_total)----------------------- it must be a ratio, not the total sigma'
                exit(0)

                chi2 = lambda *args: -2.0 * lnlike(*args)
                vaux = [ampl, 1./period, np.log(ss_total)]  # list 
                # ASI FUNCIONA BIEN MINIMIZE
                result = op.minimize( chi2, vaux, args=(tjd, flux, e_flux), \
                                          options={'disp': True}, method = 'Nelder-Mead')
                # Nelder-Mead: unconstrained minimization using Simplex algorithm (no 1st or 2nd derivate)
                # The outputs of maximum likelihood are passed as a list: 'OptimizeResult' with various
                # methods: x::shows results; sucess::whether or not optim sucess; message::descript of 
                # cause of termination...
                ampl_ml, freq_ml, lnf_ml = result.x
                print '\t\t::: elapsed time : {0} min'.format((time.time()-start_time)/60.)

                # Diagnosis plot: fit vs ML
                if False:
                    plt.close('all')
                    fig = plt.figure(); ax1 =  fig.add_subplot(211); ax2 =  fig.add_subplot(212) 
                    ax1.plot(tjd, flux,'b', color='0.7')
                    ax1.plot(tjd, ampl*np.sin(2.0*np.pi*(1./period)*tjd+ph)+delta, 'r-', label='Fit')
                    ax1.plot(tjd, ampl_ml*np.sin(2.0*np.pi*freq_ml*tjd+ph)+delta, 'g-', label='ML')
                    ax2.plot(tjd, (ampl*np.sin(2.0*np.pi*(1./period)*tjd+ph))-\
                                 (ampl_ml*np.sin(2.0*np.pi*freq_ml*tjd+ph)), 'k,')
                    plt.show()

                # END: maximum likelihood


                #           Prior-probability function
                def lnprior(vargs):
                    # lnprior=0 means prior=1.0
                    # note: np.log(0.00005)=-9.9;  np.log(1.0)=0.0
                    A, freq, lnf = vargs
                    if (0.5*ampl < A < 2.0*ampl) and (0.5*(1./period) < freq < 2.0*(1./period)) and ((np.log(0.5)+np.log(ss_total)) < lnf < (np.log(1.5)+np.log(ss_total))): #(phi-np.pi < phi < phi+np.pi)
                        return 0.0 
                    else:
                        return -np.inf


                #           Full-probability function
                def lnprob(vargs, x, y, yerr): 
                    lp = lnprior(vargs) # prior
                    # if prior is infinite, the lnprob (full log-probability) will 
                    # be infinite (negative), so probability zero.
                    if not np.isfinite(lp):
                        return -np.inf
                    # if the prior is finite, the full log-probability is the prior
                    # plus the likelihood. Remember is a log-prob.
                    else:
                        return lp + lnlike(vargs, x, y, yerr)


                #           M o n t e  C M C
                # 1) Initialize the walkers in a tiny Gaussian ball around the ML result
                ndim = 3 # number of variables to use: amplitude, period and ss_total
                perturb = e_flux+ss_total # perturbation added to gaussian random noise
                pos = [ result.x + perturb*np.random.randn(ndim) for i in range(Wks) ]
                # briefly...
                # result.x + e_flux*np.random.randn(ndim) : array of shape (ndim,)
                # pos is a list of arrays with length=ndim. pos=[ array[x,y,z], array[xx,yy,zz],...]
                
                # 2) Setup the modified Metropolis-Hasting sampler
                # with: N of walkers, ndim, full-log-prob, light-curve-data 
                sampler = emcee.EnsembleSampler( Wks, ndim, lnprob, \
                                                     args=(tjd, flux, e_flux), threads=THREADS )
                # Arguments: 
                # nwalkers (Goodman&Weare walkers), dim, lnpostfn (log-posterior probability), 
                # a=2.0 (proposal scale parameter), args=[] (list of extra positional arguments for
                # lnpostfn, it will called as lnpostfn(p,*args,**kwargs)), kwargs={} (list of extra 
                # arguments for lnpostfn, it will called as lnpostfn(p,*args,**kwargs)), postargs=None, 
                # threads=1 (number of threads used in parallel calls to lnpostfn), pool=None, 
                # live_dangerously=False, runtime_sortingfn=None
                #
                # Methods:
                # - acceptance_fraction (array length:Wks, of the fraction of steps accepted for each walker)
                # - acor (estimate of autocorr time for each parameter, length:ndim)
                # - chain (Markov chain array of shape(Wks, iterations, ndim))
                # - get_autocorr_time(window=50, fast=False) (estimate of autocorr time for each 
                #    parameter, legth:ndim. window--size of the windowing function, equivalent to 
                #    the maximum number of lags to use)
                # - get_lnprob(p) (return log-probability at the given position)
                # - run_mcmc(pos0,N,rstate0=None,lnprob0=None,**kwargs) (iterate sample() for N 
                # iterations and return result. pos0:initial posit vector, N:number of steps, 
                # lnprob0:the log-posterior proba at p0, if not, the initial value is called, ...)
                # - sample(p0, lnprob0=None, rstate0=None, blobs0=None, iterations=1, thin=1, 
                # storechain=True, mh_proposal=None) (advance the chain iterations steps as a generator)
                # ...and more....

                sampler.run_mcmc(pos, Nit) # tiny gaussian ball is the zero-position
                
                # Reshape samples array
                cut_point= 0  
                # remember: [ nwalkers, ndim, lnprob ]
                samples = sampler.chain[:, cut_point:, :].reshape((-1, ndim))

                print '\n--- Finished MCMC core.\n\t\t::: elapsed time (global):  {0:^15} min \n\t\t\t({1:^15} min)'.format((time.time()-start_time)/60., (time.time()-global_t)/60.)
                print 'Median/Mean of acceptance fraction: {0} / {1}'.format(np.median(sampler.acceptance_fraction), np.mean(sampler.acceptance_fraction))

                # Save results: 
                # - numpy array of period mcmc samples for each KIC
                # - median and mean of acceptance fraction
                # - few statistics
                
                # Save samples of rotational period and amplitude of variation
                fn_prot = 'npy_files/prot_'+kID+'_N'+str(int(Nit))+'Walk'+str(int(Wks))+'.npy'
                np.save(fn_prot, 1.0/samples[:,1])
                #
                fn_ampl = 'npy_files/ampl_'+kID+'_N'+str(int(Nit))+'Walk'+str(int(Wks))+'.npy'
                np.save(fn_ampl, samples[:,0])    
           
                # kic , median and mean of acceptance fraction
                af[0].append(kic)
                af[1].append(np.median(sampler.acceptance_fraction))
                af[2].append(np.mean(sampler.acceptance_fraction))

                # samples:: Ampl, Freq, lnf
                samples[:, 2] = np.exp(samples[:, 2])
                 # 1st and 3rd quartile statistics
                tmp = zip(*np.percentile(samples, [25,50,75], axis=0))
                A_mc, freq_mc, ss_mc =  map(lambda x: (x[1]-x[0], x[1], x[2]-x[1]), tmp) 
                 
                prot_tmp = 1.0/samples[:,1]
                aux_prot = np.percentile(prot_tmp, [25,50,75], axis=0)

                # stacking...
                comp[0].append(kic)
                comp[1].append(period)
                comp[2].append(aux_prot[1]) # median of mcmc period
                comp[3].append(aux_prot[0]) # 1st quartile of mcmc period
                comp[4].append(aux_prot[2]) # 3rd quartile of mcmc period
                comp[5].append(prot_tmp.std()) # stdev of mcmc period
                comp[6].append(ampl )
                comp[7].append(A_mc[1]) # median of mcmc amplitude
                comp[8].append(A_mc[0]) # 1st quartile of mcmc amplitude
                comp[9].append(A_mc[2]) # 3rd quartile of mcmc amplitude
                comp[10].append(samples[:,0].std()) #stdev of mcmc amplitude
                comp[11].append(ss_total)
                comp[12].append(ss_mc[1]) # median of mcmc LC noise
                comp[13].append(samples[:,2].std()) # stdev of mcmc LC noise

                print '\n\t....passing to the next LC\n_____________________'
        
    # CLOSE: walk through files
# CLOSE: walk through kics



# Write results to file

# wr_dict = {'01_kic':comp[0][:],'02_prot':comp[1][:],'03_prot_mc':comp[2][:],'04_prot_q1':comp[3][:],\
#                 '05_prot_q3':comp[4][:],'06_prot_std':comp[5][:], '07_ampl':comp[6][:],'08_ampl_mc':comp[7][:],\
#                 '09_ampl_q1':comp[8][:],'10_ampl_q3':comp[9][:],'11_ampl_std':comp[10][:],'12_noise':comp[11][:],\
#                 '13_noise_mc':comp[12][:],'14_noise_std':comp[13][:],'15_af_median':af[1][:],'16_af_mean':af[2][:]}
# DataFrame(wr_dict).to_csv('final_mcmcErr_KOI_v02_N'+str(int(Nit))+'Walk'+str(int(Wks))+'.csv', index=False, header=True)

wr_dict = {'01_kic':comp[0][:],'02_prot':comp[1][:],'03_prot_mc':comp[2][:],'04_prot_q1':comp[3][:],\
                '05_prot_q3':comp[4][:],'06_prot_std':comp[5][:], '07_ampl':comp[6][:],'08_ampl_mc':comp[7][:],\
                '09_ampl_q1':comp[8][:],'10_ampl_q3':comp[9][:],'11_ampl_std':comp[10][:],'12_noise':comp[11][:],\
                '13_noise_mc':comp[12][:],'14_noise_std':comp[13][:],'15_af_median':af[1][:],'16_af_mean':af[2][:]}
DataFrame(wr_dict).to_csv('final_mcmcErr_Plan_v02_N'+str(int(Nit))+'Walk'+str(int(Wks))+'.csv', index=False, header=True)

print '\n\n\t succesful finished! :-)'
