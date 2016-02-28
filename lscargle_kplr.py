'''Perform fast Lomb Scargle and Generalized Lomb Scargle 
on time series
'''

import fastLombScargle as fastLS
import os 
import glob          
import pandas as pd                
import numpy as np         
import time                  
import gc                       
from sys import exit           

# Phase plot
'''To Do: substitute pandas usage by numpy
'''
def PhaseD(df_lightC,Period1): 
    # function is called by a DataFrame and a float
    #           for time:          df_lc['time'].astype('float64').values
    #           for flux:           df_lc['nflux'].astype('float64').values
    #           for period:       P1
    Period1 = float(Period1)
    LONG_CAD = 29.4 # minutes, long cadence
    bin_days = (LONG_CAD / 60.) / 24. # days, 0.02042
    
    # Define the number windows to make the Phase plot
    N_window = abs( df_lightC['time'].values[-1] - \
                    df_lightC['time'].values[0] ) / Period1  
    
    # if NaN values were not drop:
    # N_window = abs( df_lightC['time'].astype('float64').values[-1] - \
    #                   df_lightC['time'].astype('float64').values[0] ) / Period1

    # Unique Time for the phase diagram. 
    
    # make a REAL axis and superpose the discrete (binned)
    #
    # make CONTINUOS axis and over this axis put the DISCRETE one
    # we have the real and the discrete: having the real, cut and paste the discrete
    
    # Dataframe to harbor the phase plot
    df_phase = None
    # Sort to avoid errors. Reset index.
    # df_lightC = df_lightC.sort(['time'], ascending=[True])
    idx_sort3 = df_lightC['time'].values[:].argsort()
    df_lightC['time'] = df_lightC['time'].values[idx_sort3]
    df_lightC['nflux'] = df_lightC['nflux'].values[idx_sort3]
    df_lightC['Q'] = df_lightC['Q'].values[idx_sort3]
    
    df_lightC = df_lightC.reset_index(drop=True)
  
    # Lets set zero the first element of time in DF
    df_lightC['time'] = df_lightC['time'] - df_lightC['time'].values[0]
  
    # Borders of the continuos axis:
    NN = 0

    # Counter of how many window phase plot are made
    counter = 0

    # Change the vision: iterate over the continuos axis up to reach the max of the 
    #                            time series array, and in the path assign the pieces to a phase plot
    
    while 1:
        Cont_up = (NN+1) * Period1
        Cont_dw = NN * Period1
        #print Cont_up, Cont_dw, Cont_up-Cont_dw
        if (Cont_up - df_lightC['time'].values[-1]) > Period1:
            break
        # value in LC of the nearest value to lower window border
        # this method returns the DF index
        D_dw_idx = (np.abs(df_lightC['time'].values - Cont_dw)).argmin()
        Discrete_dw = df_lightC['time'][ D_dw_idx ] 
        if (Discrete_dw >= Cont_dw) and (Discrete_dw < Cont_up):
            D_up_idx = (np.abs(df_lightC['time'].values - Cont_up)).argmin()
            Discrete_up = df_lightC['time'][ D_up_idx ]
            if Discrete_up > Cont_up:
                restar_ind = 1
                while 1:
                    aux_index = (np.abs(df_lightC['time'].values - Cont_up)).argmin() - restar_ind
                    aux_index -= 1
                if df_lightC['time'].values[aux_index] <= Cont_up:
                    D_up_idx = aux_index
                    Discrete_up = df_lightC['time'][D_up_idx]
                    break
                restar_ind += 1
            # Continue: offset and save in df
            # offset_up = Cont_up - df_lightC['time'].values[discrete_up]  # d2 <= c2
            offset_dw = df_lightC['time'].values[D_dw_idx] - Cont_dw # d1 >= c1
            df_tmp = df_lightC[D_dw_idx : D_up_idx+1]
            df_tmp['time'] = (df_tmp['time'] - df_tmp['time'].values[0]) + offset_dw
            df_phase = concat( [ df_phase, df_tmp ] )
            df_tmp = None
            counter += 1
        else:
            print 'no data in the range: {0}  --  {1} d'.format(Cont_dw,Cont_up)
        NN = NN + 1

    print 'Number of phase plots: {0}'.format(counter)

    # DataFrame sort() was replaced by numpy functions ONLY to
    # uniform script, because pandas.sort() worked well inside this function
    #df_phase = df_phase.sort(['time'], ascending=[True])
    #
    idx_sort2 = df_phase['time'].values[:].argsort()
    df_phase['time'] = df_phase['time'].values[idx_sort2]
    df_phase['nflux'] = df_phase['nflux'].values[idx_sort2]
    df_phase['Q'] = df_phase['Q'].values[idx_sort2]
    #
    df_phase = df_phase.reset_index(drop=True)

    # IF Normalize time:
    #time_arr = time_arr/Per1
  
    # Set DoPhasePlot to "True" in order to see phase plot
    #
    DoPhasePlot = False
    if DoPhasePlot:
        plt.plot(df_phase['time'].values, df_phase['nflux'].values, 'b,')
        plt.show()

    return df_phase, counter

def GLS(time,flux,flux_err):
    '''Statistics, Data Mining, and Machine Learning in Astronomy 
    Ivezic, Connolly, VanderPlas, and Gray
    Great job of Ivezic et al
    '''
    from astroML.time_series import lomb_scargle
    #significance levels 90%,99%,99.9%
    sig = np.array([0.1, 0.01, 0.001])
    #range of frequencies to performa analysis
    #from 0.5 to 100 days, or 0.01 to 2.0 d-1. Use 100,000 freqs
    Pomega = np.linspace(0.01, 2., 10000)
    Py,Z = lomb_scargle(time,flux,flux_err,Pomega,generalized=True,significance=sig)
    return Pomega,Py,Z 
    
    
if __name__=='__main__':
    print '\n\tStart at: {0}'.format(time.ctime())
    
    path_lcs = '/Users/fcoj/Research/16_methCompare_pazchinchon/LC2work'
    path_out = '/Users/fcoj/Research/16_methCompare_pazchinchon/ls_Me'
    grp = 'kcp_trt'
    print grp
    
    #fast LS parameters
    OVERSAMPL = 15
    HIFREQ = 2.0
    kicL,perL,freqL,powerL = np.empty([0]),np.empty([0]),np.empty([0]),np.empty([0])
    method,dT = np.empty([0]),np.empty([0])
    #iterate over files
    for fname1 in glob.glob(os.path.join(path_lcs, grp+'*')):
        kic_int = int(fname1[fname1.find('kplr')+4:fname1.find('.tsv')])
        lc = np.loadtxt(fname1)

        '''for both I'll save only 1st peak
        '''
        
        print '\tkic {0}  fastLS'.format(kic_int)
        time1 = time.time()
        '''Fast Lomb Scargle. Calculates NOT ANGULAR frequency
        px, py: periodogram axis
        nout: number of calculated frequecies
        jmax: array index corresponding to the ppal peak: argmax(py)
        prob: an estimate of the significance of that
        maximum against the hypothesis of random noise. A small value of prob
        indicates that a significant periodic signal is present.
        '''
        px,py,nout,jmax,prob = fastLS.fasper(lc[:,0],lc[:,1],OVERSAMPL,HIFREQ)
        f1 = px[jmax] #ppal peak frequency
        P1 = 1.0/f1 #ppal peak period
        elap1 = (time.time()-time1)/60. #elapsed time in minutes
        
        print '\t\t\tGLS'
        time2 = time.time()
        '''Generalized LS
        Get angular frequency
        '''
        Pgr_x,Pgr_y,Z = GLS(lc[:,0],lc[:,1],lc[:,2])
        Pgr_x = Pgr_x/(2.*np.pi) #convert angular to linear frequency
        elap2 = (time.time()-time2)/60.
        
        #Save fastLS and GLS info
        kicL = np.append(kicL,[kic_int],axis=0)
        perL = np.append(perL,[ 1.0/px[np.argsort(py)[-1]] ],axis=0)
        freqL = np.append(freqL,[ px[np.argsort(py)[-1]] ],axis=0)
        powerL = np.append(powerL,[ py[np.argsort(py)[-1]] ],axis=0)
        method = np.append(method,['fastLS'],axis=0)
        dT = np.append(dT,[elap1],axis=0)
        
        kicL = np.append(kicL,[kic_int],axis=0)
        perL = np.append(perL,[ 1.0/Pgr_x[np.argsort(Pgr_y)[-1]] ],axis=0)
        freqL = np.append(freqL,[ Pgr_x[np.argsort(Pgr_y)[-1]] ],axis=0)
        powerL = np.append(powerL,[ Pgr_y[np.argsort(Pgr_y)[-1]] ],axis=0)
        method = np.append(method,['GLS'],axis=0)
        dT = np.append(dT,[elap2],axis=0)

    kicL = kicL.astype(int)

    print 'saving data'
    fname_out = path_out+'/fLS_GLS_'+grp+'.csv'
    pd.DataFrame({'1_kic':kicL,'2_Prot':perL,'3_freq':freqL,'4_power':powerL,'5_method':method,'6_time':dT}).to_csv(fname_out,index=False, header=True)

    print '\n\tEnd at: {0}'.format(time.ctime())

else:
  print '\n\t\t NOTE: this script is was imported from another script/program'

  
'''PHASE PLOT      
Some thinhgs to remember, before call Phase plot: 
period of ppal peak:    P1
time serie w/o NaN:     TimeSer = df_lc['time_BKJD'].astype('float64').values
flux serie w/o NaN:      FluxSer = df_lc['norm_flux'].astype('float64').values
PhaseD() :                  returns a Dataframe with normalized flux and norm time
Phase, Nphases = PhaseD(df_lc, P1)
'''

