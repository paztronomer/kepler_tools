'''This script takes the Kepler fits Long cadence files (using funtion Read_lc() ), applies a function called StepFit() to fix the gaps in mean flux between quarters, and finally save LCs.
The script requires as input a list of all KIC to be processed (integer)
Outliers are removed using Iglewicz and Hoaglin criteria
Note:pandas.sort() and pandas.drop_duplicates() causes some troubles when run on some architectures (big-endian conflict with little-endian processor). So were replaced by numpy functions:
 pandas.sort() >>> idx_sort = df_lc['time'].values[:].argsort()
 pandas.drop_duplicates >>> numpy.unique(array)

'''
# Python version: 2.7n
# If use/modify, refer Francisco Paz-Chinchon 
# Have fun & play your game.
#
#_______________________________________________________________________________
# to indent in emacs: C-u C-x TAB and to select again C-x C-x
#

import os                        # walk through directories
import pandas                    # dataframe use
from pandas import *
import numpy as np               # numpy
import time                      # to see time
import gc                           # to free memory
from sys import exit             # use exit(0) to exit program
#depending of which is installed:
import pyfits as pyf             # handle fits
#import astropy.io.fits as pyf # handle fits
'''
import scipy                     # tools
from scipy.stats import nanmean  # package to calculate mean, ignoring nan
import matplotlib.pyplot as plt             # plotting tools
'''


def IHcriteria(flux_q):
    ''' Estimate Iglewicz and Hoaglin criteria for outlier 
    http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
    Formula:
    Z=0.6745(x_i - median(x)) / MAD
    if abs(z) > 3.5, x_i is a potential outlier
    Reference:
    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers", The ASQC Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    
    To have like an ace up the sleeve:
    from statsmodels.robust import scale #mad
    scale.mad(xarray, c=1, axis=0, center=np.median)
    '''
    from scipy.stats import norm

    # Percent point function (inverse of cdf -- percentiles) of a normal continous random variable
    cte = norm.ppf(0.75) #aprox 0.6745
    
    # MAD
    MAD = np.median(np.abs(flux_q-np.median(flux_q)))
    
    # Array, the Iglewicz and Hoaglin scorer.
    # If absolute value of Zscore_i is higher than 3.5 it is a possible outlier
    Zscore = np.abs( cte*(flux_q-np.median(flux_q))/MAD )
    
    if len(Zscore[Zscore<3.5])>0:
      #return a boolean array to be used for outlier removal
      return Zscore<3.5
    else:
      print '\n\tError in outlier removal'
      exit(0)
      return


def Read_lc(fileLC):
  hdu_lc = pyf.open(fileLC, memmap=True)  # hdu is a Python like list
  KIC = int( hdu_lc[0].header['KEPLERID'] )
  Q = int( hdu_lc[0].header['QUARTER']  )
  S = int( hdu_lc[0].header['SEASON']  )
  # takes the BJD of reference and adds the value to the time, and substract the delta (see Archive Manual)
  # time is in "kepler" BJD = BJD-2454833.0 
  BJD_ref = float( hdu_lc[1].header['BJDREFI']) + float(hdu_lc[1].header['BJDREFF'])-2454833.0 
  additional = [KIC, Q, S]
  timeJD = hdu_lc[1].data.field(0) + BJD_ref # !!!
  timeJD_corr = hdu_lc[1].data.field(1)
  pdc_initial = hdu_lc[1].data.field(7)
  pdc_err = hdu_lc[1].data.field(8)
    
  # to remove NaN :
  timeJD = timeJD[~np.isnan(pdc_initial)]
  timeJD_corr = timeJD_corr[~np.isnan(pdc_initial)]
  pdc = pdc_initial[~np.isnan(pdc_initial)]
  pdc_err = pdc_err[~np.isnan(pdc_initial)]
  toNorm = pdc.mean()
  
  pdc = pdc / toNorm # mean avoiding NaN, using numpy
  # In case of Nan values, can use:
  # pdc = pdc / nanmean(pdc)
  # nanmean calculates the mean of an array, without
  # consider nan values. This package is from scipy.stats
    
  # to deal with the error in the normalized flux:
  pdc_err = pdc_err / toNorm
  
  hdu_lc.close() # closes the hdu
  # time and pdc are numpy.arrays
  # I'm gonna construct a data frame for the signal, normalized by its
  # average
    
  # Outliers  removal
  # instead of use n.n*sigma, Hoaglin criteria will be employed
  # In each quarter,outliers will be removed, based on Iglewicz and Hoaglin criteria
  pdc_res = pdc[IHcriteria(pdc)]
  timeJD = timeJD[IHcriteria(pdc)]
  pdc_err = pdc_err[IHcriteria(pdc)]
  timeJD_corr = timeJD_corr[IHcriteria(pdc)]

  # Add quarter as column, to made the linear fit
  Q_arr = np.linspace(int(Q),int(Q),len(pdc_res)) #number of elemnts is tha same as time series

  df_tserie = DataFrame({'time':timeJD, 'time_corr':timeJD_corr, 'nflux':pdc_res, 'nflux_err':pdc_err, 'Q':Q_arr})
    
  return df_tserie, additional


# Fit a 1st order polinome and adjust differences
def StepFit(df_timeser):
  # receives the temporal series as a dataframe with 3 columns: time, nflux, Q
  # NaN must be already droped
  # df must be already sorted by time
  
  # An array of quarters, non duplicate items.
  nondup = np.unique(df_timeser['Q'].values)
  
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
  #return out of ELSE statement
  return df_Fit


def VerboseID(kic_int):
  kic_str = []
  kic_int = map(str,kic_int)
  for i in range(0,len(kic_int)):
    if len(kic_int[i]) == 5:
      kic_str.append('kplr0000' + kic_int[i])
    elif len(kic_int[i]) == 6:
      kic_str.append('kplr000' + kic_int[i])
    elif len(kic_int[i]) == 7:
      kic_str.append('kplr00' + kic_int[i])
    elif len(kic_int[i]) == 8:
      kic_str.append('kplr0' + kic_int[i])
    elif len(kic_int[i]) == 9:
      kic_str.append('kplr' + kic_int[i])
    else:
      print '\n\tDummy function encountered some error'
      exit(0)
  return kic_str


if __name__ == "__main__":

  general_time = time.time()
  kic_file = 'koi.txt'
  path_fits = '/dados/home/fcoj/Work/KOI/fits_Dec2013'
  path_out = '/dados/home/fcoj/Work/KOI/methComp_lcurves/'
  
  # Filename of the IDs list
  IDs = np.loadtxt(kic_file,dtype='int',comments='#')
  IDs = VerboseID(IDs)

  print '\n\tworking...\n'

  # Main sample to add to output files
  MAINSAMPL = 'd13'
  RUN_i = 'koi'
  print RUN_i
  
  # to harbor KICs thata are not in folder
  summary_ls = []

  # Walk around folders and subfolders, ID by ID
  info_aux = [[],[],[]]
  
  for index1,kplr_id in enumerate(IDs):
    # free memory
    gc.collect()
    counter_Qs = 0 # in each KIC it resets

    for (path, dirs, files) in os.walk(path_fits):
      for index2,FILE in enumerate(files):   #file is a string 
        if ("_llc.fits" in FILE) and ( str(kplr_id) in FILE ) :  
                                                                # control flux, if file contains 
                                                                # FITS in name. Only LONG cadence
          if counter_Qs == 0: # for the first
            print '\t first LC of \t {0}. time:\t {1} minutes'.format(kplr_id,(time.time()-general_time)/60.)
            df_lc = Read_lc(path_fits+'/'+FILE)[0] 
            info_aux[0].append( Read_lc(path_fits+'/'+FILE)[1][0] )
            info_aux[1].append( Read_lc(path_fits+'/'+FILE)[1][1] )
            info_aux[2].append( Read_lc(path_fits+'/'+FILE)[1][2] )
            counter_Qs += 1
          elif counter_Qs > 0:
            df_lc = concat( [ df_lc, Read_lc(path_fits+'/'+FILE)[0] ] )
            info_aux[0].append( Read_lc(path_fits+'/'+FILE)[1][0] )
            info_aux[1].append( Read_lc(path_fits+'/'+FILE)[1][1] )
            info_aux[2].append( Read_lc(path_fits+'/'+FILE)[1][2] )
            counter_Qs += 1
            # each LC is concatenated with the previous
            # this way, a single LC is constructed for
            # each ID
            
            # NOTE: partial LCs are not loaded chronologically so I MUST SORT
            
            # Up to here, the dataframe of the entire Ligth Curve is in: df_lc
            #     
    # I must discriminate between LCs with one quarter and LCs with more than one quarter,
    # in order to perform (or not):  sort, Fit, reset of indices

    # If KIC have no fits in folder
    if counter_Qs == 0:
      summary_ls.append(kplr_id)
      print 'No llc file found for {0}. Jumping to the next KIC'.format(kplr_id)
      continue
        
    # Use only quarters from 1, not including quarter 0
    df_lc = df_lc.loc[df_lc['Q'] > 0]
    
    # To know how time passed...
    if (index1+1)%100 == 0: 
      print '\t\t\t>>> LC number: {0} \t elapsed time: {1} h'.format(index1+1, (time.time()-general_time)/3600.)

    ''' CASE 1) O N L Y  O N E  QUARTER
    '''
    if counter_Qs == 1:
      print '\n\t\t>>> there is ONLY ONE quarter for {0}'.format(kplr_id)
      # Save CSV light curves
      # 3 types, the classical (comma-separated) and the full (comma-separated)
      fname = path_out + MAINSAMPL + '.' + RUN_i + '_' + kplr_id + '.csv'
      fname2 =  path_out + MAINSAMPL + '.' + RUN_i + '_' + kplr_id + '_full.csv'
      df_tmp = df_lc[['time','nflux','nflux_err']]
      df_tmp.to_csv(fname, sep=',', header=False, index=False) # without header!!!
      df_lc.to_csv(fname2, sep=',', header=True, index=False)
      
      '''CASE 2) MORE THAN ONE QUARTER
      '''
    elif counter_Qs > 1:
      # Sort LC using time
      # Due to LC is a merge of LCs with NO CHRONOLOGICAL order, we must sort
      # using time.
      idx_sort = df_lc['time'].values[:].argsort()
      df_lc['time'] = df_lc['time'].values[idx_sort]
      df_lc['nflux'] = df_lc['nflux'].values[idx_sort]
      df_lc['Q'] = df_lc['Q'].values[idx_sort]
      # I used this way instead of pandas.sort to avoid the ugly
      # problems with big-endian / little-endian
      # df_lc = df_lc.sort(['time'], ascending=[True])
        
      # Reset index
      # After sort, indices must be re-written
      # VERY IMPORTANT POINT!
      df_lc = df_lc.reset_index(drop=True)

      # Banjai
      # Fit a line (a*x+b) Q by Q
      # A linear fit is performed quarter by quarter (in pairs), and the offset is applied 
      # to the 2nd
      df_lc = StepFit( df_lc )

      # Save CSV light curves
      # 3 types, the classical (comma-separated) and the full (comma-separated)
      fname = path_out + MAINSAMPL + '.' + RUN_i + '_' + kplr_id + '.csv'
      fname2 =  path_out + MAINSAMPL + '.' + RUN_i + '_' + kplr_id + '_full.csv'
      df_tmp = df_lc[['time','nflux','nflux_err']]
      df_tmp.to_csv(fname, sep=',', header=False, index=False) # without header!!!
      df_lc.to_csv(fname2, sep=',', header=True, index=False)
      
  # Addditional info
  # To Save the additional info of each light curve
  df_AddInfo = DataFrame({'1_KIC':info_aux[0][:],'2_Quarter':info_aux[1][:],'3_Season':info_aux[2][:]})
  fn_addinfo = MAINSAMPL +'.' + RUN_i + '_pdc_LCsInfo.csv'
  df_AddInfo.to_csv(fn_addinfo, sep=',', index=False, header=True)
  # Missed KIC
  if len(summary_ls) > 0:
    fn_miss = MAINSAMPL + '.' + RUN_i + '_missed.csv'
    DataFrame({'Missed_KIC':summary_ls}).to_csv(fn_miss,index=False, header=True)    
    
  print '\t\t>>> fits2csv ended'
  print 'Total elapsed time: {0} hours'.format( (time.time()-general_time)/3600. )

else:
  print '\n\t\t NOTE: this script is was imported from another script/program'
  print '\t\t -------------------------------------------------------------------------------------'


