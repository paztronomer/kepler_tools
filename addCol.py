'''Adds error column flux to treated lcs
Also save all (treated and untreated) as space-separated values
Search pairs of tables and match them
'''

import numpy as np
import os
import sys
import glob
import matplotlib.pyplot as plt


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


def Pairs(path1,path2,kic_str):
    pathList = [[],[]]
    for i in range(0,len(kic_str)):
        for fname1 in glob.glob(os.path.join(path1, 'treat_*')):
            if kic_str[i] in fname1:
                pathList[0].append(fname1)  
        for fname2 in glob.glob(os.path.join(path2, 'd13.kcp*')):
            if kic_str[i] in fname2:
                pathList[1].append(fname2)
    return pathList


def NearestPos(arr1,value2):
    return np.argmin(np.abs(arr1-value2))

  
#function matches elements from both lists, and create updated data
def Match2(path_out,tabs2pair):
    for j in range(0,len(tabs2pair[0][:])):
        #treated cols: time|flux    .dat
        trt = np.loadtxt(tabs2pair[0][j],delimiter=' ')
        aux_fn1 = tabs2pair[0][j][ tabs2pair[0][j].find('kplr'):tabs2pair[0][j].find('.dat') ]

        #with errors cols: time|flux|flux_err    .csv
        werr = np.loadtxt(tabs2pair[1][j],delimiter=',')
        aux_fn2 = tabs2pair[1][j][ tabs2pair[1][j].find('kplr'):tabs2pair[1][j].find('.csv') ]

        print '\n\tworking on:   {0}'.format(aux_fn1)
        
        time,flux,flux_err = np.empty([0]),np.empty([0]),np.empty([0])
        for p in xrange(0,trt.shape[0]):
            time = np.append(time,[trt[p,0]],axis=0)
            flux = np.append(flux,[trt[p,1]],axis=0)
            flux_err = np.append(flux_err, [ werr[NearestPos( werr[:,0],trt[p,0] ),2] ] )

        '''After rotate array is ok, but cols must be inserted last-to-first to appear
        firs- to-last
        '''
        out1 = path_out+'kcp_trt_'+aux_fn1+'.tsv'
        nrot = 3
        np.savetxt(out1,np.rot90(np.vstack([flux_err,flux,time]),nrot),delimiter='   ')

        out2 = path_out+'kcp_raw_'+aux_fn2+'.tsv'
        np.savetxt(out2,werr,delimiter='   ')
          
    return True

  
if __name__=='__main__':
    path_treat = 's01tos04_treat'
    path_werr = 'kcp_lcs'
    path_tables = 'LC2work/'
    
    #generate list of paths, to match lists
    list2 = Pairs(path_treat,path_werr,VerboseID(np.loadtxt('kcp.txt',dtype='int')))

    #match tables
    transf = Match2(path_tables,list2)

    if transf:
      print 'All worked fine'
      
else:
    print '\n\tcalled from another script\n'
