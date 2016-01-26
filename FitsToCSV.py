'''
 Script: FitsToCSV
 Version: 04
 This script takes the Kepler fits Long cadence files (using funtion Read_lc() ), applies a function 
 called StepFit() to fix the gaps in mean flux between quarters, and finally save 2 types of 
 light curves: 
  --- a complete version (time,time_corr,nflux,nflux_corr,quarter) and 
  --- a reduced version (time, nflux)

 The script requires as input a list of all KIC to be processed, in format 'kplrNNNNNNNNN'
 Also...
 1) remove 3.5sigma outliers (flux), quarter by quarter, in the Read_lc() function
 2) pandas.sort() and pandas.drop_duplicates() causes some troubles due to pc architecture 
 (big-endian conflict with little-endian processor).
 So were replaced by numpy functions:
            pandas.sort() >>> idx_sort = df_lc['time'].values[:].argsort()
            pandas.drop_duplicates >>> numpy.unique(array)
 3) save only from quarter 1 forward, exclude quarter 0

 Python version: 2.7n

 If use/modify, refer Francisco Paz-Chinchon, francisco at dfte.ufrn.br , 
 UFRN, Brazil. 2014.

 Have fun & play your game.
'''
#_______________________________________________________________________________
# further PYFITS resources see: http://pythonhosted.org/pyfits/users_guide/users_tutorial.html#opening-a-fits-file
# to indent in emacs: C-u C-x TAB and to select again C-x C-x
#

import os                        # walk through directories
import pyfits as pyf             # handle fits
import pandas                    # dataframe use
from pandas import *
import numpy as np               # numpy
import time                      # to see time
import gc                           # to free memory
from sys import exit             # use exit(0) to exit programm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                    F U N C T I O N S
#

#                    Read LCs
#
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
  
  pdc = pdc / pdc.mean() # mean avoiding NaN, using numpy
    
  # to deal with the error in the normalized flux:
  pdc_err = pdc_err / pdc.mean()
  
  hdu_lc.close() # closes the hdu
  # time and pdc are numpy.arrays
  # I'm gonna construct a data frame for the signal, normalized by its
  # average
    
  # O u t l i e r s  removal
  # In each quarter, flux data out 3.5sigma will be erased
  # This value is arbitrary. I suggest to use the Z score based on MAD
  sigma_f = pdc.std()
  deg3_fit = np.poly1d( np.polyfit(timeJD, pdc, 3) )
  # evaluate in all the flux array. Abs value
  pdc_diff = np.abs( deg3_fit(timeJD)-pdc )
  pdc = pdc[pdc_diff < 3.5*sigma_f]
  timeJD = timeJD[pdc_diff < 3.5*sigma_f]
  pdc_err = pdc_err[pdc_diff < 3.5*sigma_f]
  timeJD_corr = timeJD_corr[pdc_diff < 3.5*sigma_f]

  # Add quarter as a third column, to made the linear fit
  Q_arr = np.linspace(int(Q),int(Q),len(pdc)) #number of elemnts is tha same as time series
  
  # df_tserie = DataFrame({'1.time':timeJD, '2.time_corr':timeJD_corr, '3.nflux':pdc, '4.nflux_err':pdc_err, '5.Q':Q_arr})
  df_tserie = DataFrame({'time':timeJD, 'time_corr':timeJD_corr, 'nflux':pdc, 'nflux_err':pdc_err, 'Q':Q_arr})
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


# <><><><><><><><><><><><><><><><><><><><>
#
#                C O R P U S
#
# <><><><><><><><><><><><><><><><><><><><>
if __name__ == "__main__":

  general_time = time.time()

  # Superior-Path to files
  path_fits = 'somewhere_folder'
  
  # Filename of the IDs list
  df_IDs = read_table('IDS_list.csv',sep=',') #list of type kplr00xxxxxxx

  print '\n\tworking...\n'

  # Main sample to add to output files
  MAINSAMPL = 'sample01'
  RUN_i = 'run01'
  
  # to harbor KICs that are not found in the folder
  summary_ls = []

  # Walk around folders and subfolders, ID by ID
  info_aux = [[],[],[]]
  
  for index1,kplr_id in enumerate( df_IDs['kic_name'].values ):
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
            
            # NOTE: partial LCs are not charged chronologically so I MUST SORT THEM!
            
            # Up to here, the dataframe of the entire Ligth Curve is in: df_lc
            # Remember it has 3 columns: flux, time, quarter
            #     
    # I must discriminate between LCs with one querter and LCs with more than one quarter,
    # in order to perform (or not):  sort, Fit, reset of indices

    # Use only quarters from 1, not including quarter 0
    df_lc = df_lc.loc[df_lc['Q'] > 0]
    
    # If KIC have no fits in folder
    if counter_Qs == 0:
      summary_ls.append(kplr_id)
    
    # To know how time passed...
    if (index1+1)%100 == 0: 
      print '\t\t\t>>> LC number: {0} \t elapsed time: {1} h'.format(index1+1, (time.time()-general_time)/3600.)

    #
    #          CASE 1) O N L Y  O N E  QUARTER
    #
    if counter_Qs == 1:
      print '\n\t\t>>> there is ONLY ONE quarter for {0}'.format(kplr_id)
      #                       MAKE UP
      #                       ------------
      #           1.- Sort LC using time
      # Do not need to sort because we assume the time series is already sorted in the quarter
        
      #           2.- Reset index
      # As no sort is performed, we not need to reset indices

      #           3.- Fit a line (a*x+b) Q by Q
      # As only one quarter is present, no Fit can be done between quarters

      #           4.- Save CSV light curves
      # 2 types, the classical and the full
      fname =  MAINSAMPL + '.' + RUN_i + '_' + kplr_id + '.dat'
      fname2 = MAINSAMPL + '.' + RUN_i + '_' + kplr_id + '_full.csv'
      df_tmp = df_lc[['time', 'nflux']]
      df_tmp.to_csv(fname, sep=' ', header=False, index=False) # without header!!!
      df_lc.to_csv(fname2, sep=',', header=True, index=False)
    
      #...................................................................................................................



    #
    #          CASE 2) M O R E  THAN  O N E  QUARTER
    #
    elif counter_Qs > 1:
      #           1.- Sort LC using time
      # Due to LC is a merge of LCs with NO CHRONOLOGICAL order, we must sort
      # using time.
      idx_sort = df_lc['time'].values[:].argsort()
      df_lc['time'] = df_lc['time'].values[idx_sort]
      df_lc['nflux'] = df_lc['nflux'].values[idx_sort]
      df_lc['Q'] = df_lc['Q'].values[idx_sort]
      # I used this way instead of pandas.sort to avoid the ugly
      # problems with big-endian / little-endian
      # df_lc = df_lc.sort(['time'], ascending=[True])
        
      #           2.- Reset index
      # After sort, indices must be re-written
      # VERY IMPORTANT POINT!
      df_lc = df_lc.reset_index(drop=True)

      #           3.- Fit a line (a*x+b) Q by Q
      # A linear fit is performed quarter by quarter (in pairs), and the offset is applied 
      # to the 2nd
      df_lc = StepFit( df_lc )

      #          4.- Save CSV light curves
      # 2 types, the classical and the full
      fname = MAINSAMPL +  '.' + RUN_i + '_' + kplr_id + '.dat'
      fname2 = MAINSAMPL +  '.' + RUN_i + '_' + kplr_id + '_full.csv'
      df_tmp = df_lc[['time', 'nflux']]
      df_tmp.to_csv(fname, sep=' ', header=False, index=False) # without header!!!
      df_lc.to_csv(fname2, sep=',', header=True, index=False)
      
  # closing the FOR... 
  # OUT OF FOR THAT WALK THROUGH KICs...

  #                        W R I T E  I N F O  T O  F I L E S (II)
  #
  # GENERAL:   To Save the additional info of each light curve
  df_AddInfo = DataFrame({'1_KIC':info_aux[0][:],'2_Quarter':info_aux[1][:],'3_Season':info_aux[2][:]})
  fn_addinfo = MAINSAMPL +'.' + RUN_i + '_pdc_LCsInfo.csv'
  df_AddInfo.to_csv(fn_addinfo, sep=',', index=False, header=True)
  
  # GENERAL: missed KIC
  if len(summary_ls) > 0:
    fn_miss = MAINSAMPL + '.' + RUN_i + '_FitsToCSV_noKIC.csv'
    DataFrame({'Missed_KIC':summary_ls}).to_csv(fn_miss,index=False, header=True)
  #...................................................................................................................
    
    
  print '\t\t>>>FITS to CSV ended!!!'
  print 'Total elapsed time: {0} hours'.format( (time.time()-general_time)/3600. )



else:
  print '\n\t\t NOTE: this script is was imported from another script/program'
  print '\t\t -------------------------------------------------------------------------------------'


