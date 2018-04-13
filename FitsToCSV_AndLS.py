'''
 Script: FitsToCSV_AndLS
 Version: v09
 -- This scripts takes FITS files (for the available quarters), remove 3.5sigma
 outliers, join quarters using linear fit, and calculate Fast Lomb Scargle periodogram.
 Display Phase diagram, by using the 1st period (peak 1 of Lomb Scargle periodogram)
         The scripts returns the LS for each light curve, up to 100th peak, and a resume
         table with periodogram info up to 4th peak, for all analyzed LCs.
         Also a list of KICs that were given as input but were not found in the target folder.

 This script is similar to FitsToCSV.py script, in this repo, but has an additional feature of Lomb Scargle
 Python version: 2.7n

 If use/modify, refer Francisco Paz-Chinchon, francisco at dfte.ufrn.br ,
 UFRN, Brazil. 2014.

 Have fun & play your game.
'''
#_______________________________________________________________________________
# further PYFITS resources see: http://pythonhosted.org/pyfits/users_guide/users_tutorial.html#opening-a-fits-file
#

import os                        # walk through directories
import pyfits as pyf             # handle fits
import pandas                    # dataframe use
from pandas import *
import numpy as np               # numpy
import time                      # to see time
import gc                           # to free memory
from sys import exit             # use exit(0) to exit programm

#...........
#......................
#.................................
#.......BEGIN OF LS ROUTINES.........

# """ Fast algorithm for spectral analysis of unevenly sampled data

# The Lomb-Scargle method performs spectral analysis on unevenly sampled
# data and is known to be a powerful way to find, and test the 
# significance of, weak periodic signals. The method has previously been
# thought to be 'slow', requiring of order 10(2)N(2) operations to analyze
# N data points. We show that Fast Fourier Transforms (FFTs) can be used
# in a novel way to make the computation of order 10(2)N log N. Despite
# its use of the FFT, the algorithm is in no way equivalent to 
# conventional FFT periodogram analysis.

# Keywords:
#   DATA SAMPLING, FAST FOURIER TRANSFORMATIONS, 
#   SPECTRUM ANALYSIS, SIGNAL  PROCESSING

# Example:
#   > import numpy
#   > import lomb
#   > x = numpy.arange(10)
#   > y = numpy.sin(x)
#   > fx,fy, nout, jmax, prob = lomb.fasper(x,y, 6., 6.)

# Reference: 
#   Press, W. H. & Rybicki, G. B. 1989
#   ApJ vol. 338, p. 277-280.
#   Fast algorithm for spectral analysis of unevenly sampled data
#   bib code: 1989ApJ...338..277P

# """
from numpy import *
from numpy.fft import *

def __spread__(y, yy, n, x, m):
  # """
  # Given an array yy(0:n-1), extirpolate (spread) a value y into
  # m actual array elements that best approximate the "fictional"
  # (i.e., possible noninteger) array element number x. The weights
  # used are coefficients of the Lagrange interpolating polynomial
  # Arguments:
  #   y : 
  #   yy : 
  #   n : 
  #   x : 
  #   m : 
  # Returns:
    
  # """
  nfac=[0,1,1,2,6,24,120,720,5040,40320,362880]
  if m > 10. :
    print 'factorial table too small in spread'
    return

  ix=long(x)
  if x == float(ix): 
    yy[ix]=yy[ix]+y
  else:
    ilo = long(x-0.5*float(m)+1.0)
    ilo = min( max( ilo , 1 ), n-m+1 ) 
    ihi = ilo+m-1
    nden = nfac[m]
    fac=x-ilo
    for j in range(ilo+1,ihi+1): fac = fac*(x-j)
    yy[ihi] = yy[ihi] + y*fac/(nden*(x-ihi))
    for j in range(ihi-1,ilo-1,-1):
      nden=(nden/(j+1-ilo))*(j-ihi)
      yy[j] = yy[j] + y*fac/(nden*(x-j))

def fasper(x,y,ofac,hifac, MACC=4):
  # """ function fasper
  #   Given abscissas x (which need not be equally spaced) and ordinates
  #   y, and given a desired oversampling factor ofac (a typical value
  #   being 4 or larger). this routine creates an array wk1 with a
  #   sequence of nout increasing frequencies (not angular frequencies)
  #   up to hifac times the "average" Nyquist frequency, and creates
  #   an array wk2 with the values of the Lomb normalized periodogram at
  #   those frequencies. The arrays x and y are not altered. This
  #   routine also returns jmax such that wk2(jmax) is the maximum
  #   element in wk2, and prob, an estimate of the significance of that
  #   maximum against the hypothesis of random noise. A small value of prob
  #   indicates that a significant periodic signal is present.
  
  # Reference: 
  #   Press, W. H. & Rybicki, G. B. 1989
  #   ApJ vol. 338, p. 277-280.
  #   Fast algorithm for spectral analysis of unevenly sampled data
  #   (1989ApJ...338..277P)

  # Arguments:
  #     X   : Abscissas array, (e.g. an array of times).
  #     Y   : Ordinates array, (e.g. corresponding counts).
  #     Ofac : Oversampling factor.
  #     Hifac : Hifac * "average" Nyquist frequency = highest frequency
  #          for which values of the Lomb normalized periodogram will
  #          be calculated.
      
  #  Returns:
  #     Wk1 : An array of Lomb periodogram frequencies.
  #     Wk2 : An array of corresponding values of the Lomb periodogram.
  #     Nout : Wk1 & Wk2 dimensions (number of calculated frequencies)
  #     Jmax : The array index corresponding to the MAX( Wk2 ).
  #     Prob : False Alarm Probability of the largest Periodogram value
  #     MACC : Number of interpolation points per 1/4 cycle
  #           of highest frequency

  # History:
  #   02/23/2009, v1.0, MF
  #     Translation of IDL code (orig. Numerical recipies)
  # """
  #Check dimensions of input arrays
  n = long(len(x))
  if n != len(y):
    print 'Incompatible arrays.'
    return

  nout  = 0.5*ofac*hifac*n
  nfreqt = long(ofac*hifac*n*MACC)   #Size the FFT as next power
  nfreq = 64L             # of 2 above nfreqt.

  while nfreq < nfreqt: 
    nfreq = 2*nfreq

  ndim = long(2*nfreq)
  
  #Compute the mean, variance
  ave = y.mean()
  ##sample variance because the divisor is N-1
  var = ((y-y.mean())**2).sum()/(len(y)-1) 
  # and range of the data.
  xmin = x.min()
  xmax = x.max()
  xdif = xmax-xmin

  #extirpolate the data into the workspaces
  wk1 = zeros(ndim, dtype='complex')
  wk2 = zeros(ndim, dtype='complex')

  fac  = ndim/(xdif*ofac)
  fndim = ndim
  ck  = ((x-xmin)*fac) % fndim
  ckk  = (2.0*ck) % fndim

  for j in range(0L, n):
    __spread__(y[j]-ave,wk1,ndim,ck[j],MACC)
    __spread__(1.0,wk2,ndim,ckk[j],MACC)

  #Take the Fast Fourier Transforms
  wk1 = ifft( wk1 )*len(wk1)
  wk2 = ifft( wk2 )*len(wk1)

  wk1 = wk1[1:nout+1]
  wk2 = wk2[1:nout+1]
  rwk1 = wk1.real
  iwk1 = wk1.imag
  rwk2 = wk2.real
  iwk2 = wk2.imag
  
  df  = 1.0/(xdif*ofac)
  
  #Compute the Lomb value for each frequency
  hypo2 = 2.0 * abs( wk2 )
  hc2wt = rwk2/hypo2
  hs2wt = iwk2/hypo2

  cwt  = sqrt(0.5+hc2wt)
  swt  = sign(hs2wt)*(sqrt(0.5-hc2wt))
  den  = 0.5*n+hc2wt*rwk2+hs2wt*iwk2
  cterm = (cwt*rwk1+swt*iwk1)**2./den
  sterm = (cwt*iwk1-swt*rwk1)**2./(n-den)

  wk1 = df*(arange(nout, dtype='float')+1.)
  wk2 = (cterm+sterm)/(2.0*var)
  pmax = wk2.max()
  jmax = wk2.argmax()


  #Significance estimation
  #expy = exp(-wk2)          
  #effm = 2.0*(nout)/ofac       
  #sig = effm*expy
  #ind = (sig > 0.01).nonzero()
  #sig[ind] = 1.0-(1.0-expy[ind])**effm

  #Estimate significance of largest peak value
  expy = exp(-pmax)          
  effm = 2.0*(nout)/ofac       
  prob = effm*expy

  if prob > 0.01: 
    prob = 1.0-(1.0-expy)**effm

  return wk1,wk2,nout,jmax,prob

def getSignificance(wk1, wk2, nout, ofac):
  # """ returns the peak false alarm probabilities
  # Hence the lower is the probability and the more significant is the peak
  # """
  expy = exp(-wk2)          
  effm = 2.0*(nout)/ofac       
  sig = effm*expy
  ind = (sig > 0.01).nonzero()
  sig[ind] = 1.0-(1.0-expy[ind])**effm
  return sig
#...........
#......................
#.................................
#............END OF LS ROUTINES........


#
#                                 FUNCTIONS
#


#                    Read LCs
#
def Read_lc(fileLC):
  hdu_lc = pyf.open(fileLC, memmap=True)  # hdu is a Python like list
  KIC = int( hdu_lc[0].header['KEPLERID'] )
  Q = int( hdu_lc[0].header['QUARTER']  )
  S = int( hdu_lc[0].header['SEASON']  )
  additional = [KIC, Q, S]
  timeJD = hdu_lc[1].data.field(0)
  pdc_initial = hdu_lc[1].data.field(7)
  # to remove NaN :
  timeJD = timeJD[~np.isnan(pdc_initial)]
  pdc = pdc_initial[~np.isnan(pdc_initial)]
  pdc = pdc / pdc.mean() # mean avoiding NaN, using numpy
  
  hdu_lc.close() # closes the hdu
  # time and pdc are numpy.arrays
  # I'm gonna construct a data frame for the signal, normalized by its
  # average

  # O u t l i e r s  removal
  # In each quarter, flux data out 3.5sigma will be erased
  # This is an arbitrary value (3.5). I highly recommend to use Iglewicz & Hoaglin
  # Z-score (see http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm)
  sigma_f = pdc.std()
  deg3_fit = np.poly1d( np.polyfit(timeJD, pdc, 3) )
  # evaluate in all the flux array. Abs value
  pdc_diff = np.abs( deg3_fit(timeJD) - pdc )
  pdc = pdc[pdc_diff < 3.5*sigma_f]
  timeJD = timeJD[pdc_diff < 3.5*sigma_f]

  # Add quarter as a third column, to made the linear fit
  Q_arr = np.linspace(int(Q),int(Q),len(pdc)) #number of elemnts is tha same as time series
  
  # To fix some error bigendian little-endian that would appear due to 
  # computer architecture
  df_tserie = DataFrame({'time':timeJD, 'nflux':pdc, 'Q':Q_arr})
  # note that the flux is already normalized, so
  # the quartes can be simply joined into one
    
  return df_tserie, additional


#                    Fit a 1st order polinome and adjust differences
#
def StepFit(df_timeser):
  # receives the temporal series as a dataframe with 3 columns: time, nflux, Q
  # NaN must be already droped
  # df must be already sorted by time
  
  # An array of quarters, non duplicate items.
  # nondup = df_timeser.drop_duplicates(cols=['Q'] )['Q'].values # only 'Q' column
  nondup = np.unique(df_timeser['Q'].values)
  print nondup
  #print nondup, len(nondup)
  
  # Go through time series, quarter by quarter
  if len(nondup) > 1:
    # Save the first quarter of the LC
    df_Fit = df_timeser[ df_timeser.Q == nondup[0] ]
      
    # Iterate up to the n-1 element
    for index,q_item in enumerate( nondup[:len(nondup)-1] ): # indexing is OK!!!
      df_tmp1 = df_timeser[ df_timeser.Q == q_item ]
      df_tmp2 = df_timeser[ df_timeser.Q == nondup[index+1] ]
      # fit the 2 linear fits using: np.polyfit, np.polyval, p = np.poly1d(z)
      p1 = np.poly1d( np.polyfit(df_tmp1['time'].values, df_tmp1['nflux'].values,1) )
      p2 = np.poly1d( np.polyfit(df_tmp2['time'].values, df_tmp2['nflux'].values,1) )
      # then evaluate the borders of each piece, in the corresponding fit
      # and determine the offset.
      Offset = p1(df_tmp1['time'].values[-1]) -  p2(df_tmp2['time'].values[-1])
      # now add the offset to the second piece
      df_tmp2['nflux'] += Offset
      # and concatenate the 2nd piece with the previous
      df_Fit = concat( [df_Fit, df_tmp2] )
  else:
    df_Fit = df_timeser
    print 'no fit made, only ONE quarter in LC'

  return df_Fit


#                      Phase plot
#
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
  # Ther is an unique Time for the phase diagram.
 
  # Dataframe to harbor the phase plot
  df_phase = None
  # Sort to avoid errors. Reset index.
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
  
  while 1:
    Cont_up = (NN+1) * Period1
    Cont_dw = NN * Period1
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
  # DoPhasePlot = False
  if DoPhasePlot:
    plt.plot(df_phase['time'].values, df_phase['nflux'].values, 'b,')
    plt.show()

  return df_phase, counter
  #
  # NOTE: time output isn't normalized [0,1] to reach desired 
  # precission in SIGMA calculation
  #



# <><><><><><><><><><><><><><><><><><><><>
#
#                C O R P U S
#
# <><><><><><><><><><><><><><><><><><><><>
if __name__ == "__main__":

  general_time = time.time()

  # Superior-Path to files
  path_fits = 'fits_folder'
  
  # Filename of the IDs list
  df_IDs = read_table('IDS_kplr.csv',sep=',') #list of type kplr00xxxxxxx

  print '\n\tworking...\n'

  # Walk around folders and subfolders, ID by ID
  info_aux = [[],[],[]]
  summary_ls = [] # summary list, wil show the LCs that aren't availables in cluster folders

  # list to save the KIC and the ppal period
  # periods_tab = [['KICname'],['Period_(d)'],['LinearFreq'],['Power'], ['Prob'],['Peak']]
  periods_tab = [[],[],[],[],[],[]]
  
  # See Paz-Chinchon et al 2015 for justification
  OVERSAMPL = 15
  HIFREQ = 2.0

  MAINSAMPL = 'sample01'
  RUN_i = 'run01'
  
  for index1,kplr_id in enumerate( df_IDs['kic_name'].values ):
    # free memory
    # gc.collect()
    counter_Qs = 0 # in each KIC it resets

    for (path, dirs, files) in os.walk(path_fits):
      for index2,FILE in enumerate(files):   #file is a string 
        if ("_llc.fits" in FILE) and ( str(kplr_id) in FILE ) :  
                                                                # control flux, if file contains 
                                                                # FITS in name. Only LONG cadence
          if counter_Qs == 0: # for the first
            print '\t first LC of \t {0}. time:\t {1} minutes'.format(kplr_id,(time.time()-general_time)/60.)
            df_lc = Read_lc(path_fits+'/'+FILE)[0]
            aux_1 = Read_lc(path_fits+'/'+FILE)[1]
            info_aux[0].append( aux_1[0] )
            info_aux[1].append( aux_1[1] )
            info_aux[2].append( aux_1[2] )
            counter_Qs += 1
          elif counter_Qs > 0:
            df_lc = concat( [ df_lc, Read_lc(path_fits+'/'+FILE)[0] ] )
            aux_2 = Read_lc(path_fits+'/'+FILE)[1]
            info_aux[0].append( aux_2[0] )
            info_aux[1].append( aux_2[1] )
            info_aux[2].append( aux_2[2] )
            counter_Qs += 1
            # each LC is concatenated with the previous
            # this way, a single LC is constructed for
            # each ID
            
            #### NOTE: partial LCs are not charged chronologically so I MUST SORT
            
            # Up to here, the dataframe of the entire Ligth Curve is in: df_lc
            # Remember it has 3 columns: flux, time, quarter
            #     
    # I must discriminate between LCs with one querter and LCs with more than one quarter,
    # in order to perform (or not):  sort, Fit, reset of indices
    
    # If KIC have no fits in folder
    if (counter_Qs == 0) or (len(df_lc['nflux'].values) == 0):
      summary_ls.append(kplr_id)

    # To know how time passed...
    if (index1+1)%100 == 0: 
      print '\t\t\t>>> LC number: {0} \t elapsed time: {1} h'.format(index1+1, (time.time()-general_time)/3600.)

    #
    #          CASE 1) O N L Y  O N E  QUARTER
    #
    if (counter_Qs == 1) and (len(df_lc['nflux'].values) > 4):
      print '\n\t\t>>> there is ONLY ONE quarter for {0}'.format(kplr_id)
      #                      FAST LOMB SCARGLE
      #                      ------------------------------
      #                      (we're inside the loop of kepler IDs...)
      #
      # calculates NOT ANGULAR frequency
      #    px, py: periodogram axis
      #    nout: number of calculated frequecies
      #    jmax: array index corresponding to the ppal peak: argmax(py)
      #    prob: an estimate of the significance of that
      #      maximum against the hypothesis of random noise. A small value of prob
      #      indicates that a significant periodic signal is present.
    
      oversampl = OVERSAMPL # 4 or more
      hifreq = HIFREQ        # 1 or more
      px,py,nout,jmax,prob = fasper(df_lc['time'].values, \
                                      df_lc['nflux'].values, oversampl, hifreq)
      f1 = px[jmax]  # PPAL PEAK FREQUENCY
      P1 = 1.0/f1      # PPAL PEAK PERIOD
      #
      #
      DO_info_LS = False
      if DO_info_LS:
        print '\n                Fast LombScargle'
        print 'Frequency of the ppal peak (not angular freq): {}'.format(px[jmax])
        print 'Freq power, peak A: {0}'.format(py[jmax])
        print 'Period of peak A: {0} days'.format(1/px[jmax])
        print 'Probability (lower, more persistent signal): {0}'.format(prob)
        print 'time elapsed up to LS... {0} sec'.format(time.time()-general_time)
        print 
    
      #           Save periodogram
      # I use pandas.sort() because contain the output np.array from LS
      df_LS = DataFrame({'1_Period(d)':1.0/px,'2_LinearFreq':px,'3_Power':py}) # all arrays must have the same length
      df_LS = df_LS.sort(['3_Power'], ascending=[False])    # sort by power of each period
      df_LS = df_LS.reset_index(drop=True)                # re-index data frame after sort
      N_freqs = nout
      print '\n\tnumber of calculated frequencies: {0}\n'.format(N_freqs)

      # To save up to 4th peak
      for pp in [0,1,2,3]:
        periods_tab[0].append(kplr_id[ kplr_id.find('kplr')+4: ]) # save the KIC
        periods_tab[1].append(df_LS['1_Period(d)'].values[pp])
        periods_tab[2].append(df_LS['2_LinearFreq'].values[pp])
        periods_tab[3].append(df_LS['3_Power'].values[pp]) # ppal Power, as sort ascending is the first
        periods_tab[4].append(prob) # v08
        periods_tab[5].append(pp+1) # v09
    
      #           XX.- erase variables
      px,py,nout,jmax,prob = None,None,None,None,None
      #...................................................................................................................

      #                        W R I T E  I N F O  T O  F I L E S (I)
      #                        ----------------------------------
      #
      # KIC by KIC:   To Save periodogram first 100 frequencies
      out_name = MAINSAMPL +'.' + RUN_i + '_LS_' + kplr_id + '.csv'
      df_LS[0:100].to_csv(out_name, sep=',', index=False, header=True)
      #...................................................................................................................


    #
    #          CASE 2) M O R E  THAN  O N E  QUARTER
    #
    elif (counter_Qs > 1) and (len(df_lc['nflux'].values) > 4): # by Kolmogorov
      #                       MAKE UP
      #                       ------------
      #           1.- Sort LC using time
      # Due to LC is a merge of LCs with NO CHRONOLOGICAL order, we must sort
      # using time.
      idx_sort = df_lc['time'].values[:].argsort()
      df_lc['time'] = df_lc['time'].values[idx_sort]
      df_lc['nflux'] = df_lc['nflux'].values[idx_sort]
      df_lc['Q'] = df_lc['Q'].values[idx_sort]
      # I used this way instead of pandas.sort to avoid the ugly
      # problems with big-endian / little-endian
      #           2.- Reset index
      # After sort, indices must be re-written
      # VERY IMPORTANT POINT!
      df_lc = df_lc.reset_index(drop=True)
      #           3.- Fit a line (a*x+b) Q by Q
      # A linear fit is performed quarter by quarter (in pairs), and the offset is applied 
      # to the 2nd
      df_lc = StepFit( df_lc )
      # Note that: 
      #              When use np.insfinite()...
      #              In order to recover data from dataframe, in format float, we must use:
      #
      #          Time = df_lcNaN['time_BKJD'].astype('float64').values    
      #          Flux = df_lcNaN['norm_flux'].astype('float64').values    
      #...................................................................................................................
      
      #                      FAST LOMB SCARGLE
      #                      ------------------------------
      #                      (we're inside the loop of kepler IDs...)
      #
      # calculates NOT ANGULAR frequency
      #    px, py: periodogram axis
      #    nout: number of calculated frequecies
      #    jmax: array index corresponding to the ppal peak: argmax(py)
      #    prob: an estimate of the significance of that
      #      maximum against the hypothesis of random noise. A small value of prob
      #      indicates that a significant periodic signal is present.
    
      oversampl = OVERSAMPL # 4 or more
      hifreq = HIFREQ        # 1
      px,py,nout,jmax,prob = fasper(df_lc['time'].values, \
                                      df_lc['nflux'].values, oversampl, hifreq)
      f1 = px[jmax]  # PPAL PEAK FREQUENCY
      P1 = 1.0/f1      # PPAL PEAK PERIOD
      #
      #
      DO_info_LS = False
      if DO_info_LS:
        print '\n                Fast LombScargle'
        print 'Frequency of the ppal peak (not angular freq): {}'.format(px[jmax])
        print 'Freq power, peak A: {0}'.format(py[jmax])
        print 'Period of peak A: {0} days'.format(1/px[jmax])
        print 'Probability (lower, more persistent signal): {0}'.format(prob)
        print 'time elapsed up to LS... {0} sec'.format(time.time()--general_time)
        print 
    
      #           Save periodogram
      # I use pandas.sort() because contain the output np.array from LS
      df_LS = DataFrame({'1_Period(d)':1.0/px,'2_LinearFreq':px,'3_Power':py}) # all arrays must have the same length
      df_LS = df_LS.sort(['3_Power'], ascending=[False])    # sort by power of each period
      df_LS = df_LS.reset_index(drop=True)                # re-index data frame after sort
      N_freqs = nout
      print '\n\tnumber of calculated frequencies: {0}\n'.format(N_freqs)

      # To save up to 4th peak
      for pp in [0,1,2,3]:
        periods_tab[0].append(kplr_id[ kplr_id.find('kplr')+4: ]) # save the KIC
        periods_tab[1].append(df_LS['1_Period(d)'].values[pp])
        periods_tab[2].append(df_LS['2_LinearFreq'].values[pp])
        periods_tab[3].append(df_LS['3_Power'].values[pp]) # ppal Power, as sort ascending is the first
        periods_tab[4].append(prob) # v08
        periods_tab[5].append(pp+1) # v09


      #           XX.- erase variables
      px,py,nout,jmax,prob = None,None,None,None,None
      #...................................................................................................................

      #                        W R I T E  I N F O  T O  F I L E S (I)
      #                        ----------------------------------
      #
      # KIC by KIC:   To Save periodogram first 100 frequencies
      out_name = 'tables/LS/'+ MAINSAMPL +'.' + RUN_i + '_LS_' + kplr_id + '.csv'
      df_LS[0:100].to_csv(out_name, sep=',', index=False, header=True)
      #...................................................................................................................
      #                      PHASE PLOT
      #                      -----------------
      #
      # Some thinhgs to remember, before call Phase plot: 
      #           period of ppal peak:    P1
      #           time serie w/o NaN:     TimeSer = df_lc['time_BKJD'].astype('float64').values
      #           flux serie w/o NaN:      FluxSer = df_lc['norm_flux'].astype('float64').values
      #           PhaseD() :                  returns a Dataframe with normalized flux and norm time
      Phase, Nphases = PhaseD(df_lc, P1)
      print 'elapsed time up to 1st phase plot: {0} sec'.format(abs(general_time - time.time()))
      #...................................................................................................................


  # closing the FOR...
  # OUT OF FOR THAT WALK THROUGH KICs...

  #                        W R I T E  I N F O  T O  F I L E S (II)
  #
  # GENERAL:    To save first peak info and KIC
  df_Ppal = DataFrame({'1_KIC':periods_tab[0][:], '2_Period(d)':periods_tab[1][:], '3_LinearFreq':periods_tab[2][:], \
                       '4_Power':periods_tab[3][:],'5_FAP_Prob':periods_tab[4][:], '6_Peak':periods_tab[5][:]})
  fname_Ppal = MAINSAMPL +'.' + RUN_i + '_4PeakInfo.csv'
  df_Ppal.to_csv(fname_Ppal, sep=',', index=False, header=True)
    
  # GENERAL:   To Save the additional info of each light curve
  df_AddInfo = DataFrame({'1_KIC':info_aux[0][:],'2_Quarter':info_aux[1][:],'3_Season':info_aux[2][:]})
  fname_AddInfo = MAINSAMPL +'.' + RUN_i + '_LCsInfo.csv'
  df_AddInfo.to_csv(fname_AddInfo, sep=',', index=False, header=True)
  
  # GENERAL: missed KIC
  if len(summary_ls) > 0:
    fn_miss = MAINSAMPL +'.' + RUN_i + '_noKIC.csv'
    DataFrame({'Missed_KIC':summary_ls}).to_csv(fn_miss,index=False, header=True)
  #...................................................................................................................
    
    
  print 'Thanks... The End!'
  print 'Total elapsed time: {0} hours'.format( (time.time()-general_time)/3600. )



else:
  print '\n\t\t NOTE: this script is was imported from another script/program'
  print '\t\t -------------------------------------------------------------------------------------'


