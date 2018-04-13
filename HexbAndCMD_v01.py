# Trying to integrate hexbin with CMD plots
#
# If use/reference, cite: Francisco Paz-Chinchon
#                         francisco at dfte.ufrn.br
#                         DFTE-UFRN, Natal, Brazil
#

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpacks
import rpy2.rinterface as rinterface
# to handle garbage
import gc
#
import numpy as np
from sys import exit
# charge libraries to PlotCMD
import matplotlib as mpl
import matplotlib.pyplot as plt
# charge pandas to SaveTab
from pandas import *
# to walk around files
import os


#____________
#       TIPS:
# (1) to know from where packages come: print rpacks.wherefrom('lm')
# (2) to call a script from the interpreter
#
# This solution works well on R:
# aux <- sapply(c("count", "xcm", "ycm", "cell"), function(x) slot(hbin,x))
# write.table(aux, file="NGC6584.0092mag.sub", row.names=FALSE,col.names=TRUE )


#==================================================================================
#
#          Defining names
#
base = importr('base')
utils = importr('utils')
hexbin = importr('hexbin')
help_where = utils.help_search("help")
#print (robjects.r)
#print help_where


#==================================================================================
#
#          Functions
#
def ReadTab(fname,Ndiv):
  table = utils.read_table(fname,comment_char='#')#,colClasses="numeric"
  print '\nstructure type containing table: {0}'.format(type(table))
  # is a rpy2.robjects.vectors.DataFrame
  # each dimension is a column of the original table
  
  # 
  # Run hexbin
  v_i, v, i = table[0], table[1], table[2]
  #
  #
  # import pylab as plab
  # plab.plot(v_i, v, 'ro')
  # plab.show()
  #
  #
  hbin = hexbin.hexbin(v_i, v, xbins=Ndiv)
  #
  # Traduce result to rpy2 vectors:
  count = robjects.r.slot(hbin, 'count')
  xcm = robjects.r.slot(hbin, 'xcm')
  ycm = robjects.r.slot(hbin, 'ycm')
  cell = robjects.r.slot(hbin, 'cell')
  print '\nclass of the created R vectors: {0}'.format(count.rclass)
  print 'length of the R-vectors: {0}, {1}, {2}, {3}'.format(len(count),
						    len(xcm), len(ycm),
						    len(cell))
  #
  # free memory from R
  gc.collect()
  #
  # Fill a python array with these resulting rpy2 vectors:
  hx_arr_tmp = np.empty((len(count),4))
  hx_arr_tmp[:,:] = np.nan
  hx_arr_tmp[:,0], hx_arr_tmp[:,1], hx_arr_tmp[:,2], hx_arr_tmp[:,3] = xcm, ycm, count, cell
  #
  return hx_arr_tmp
  #____________

def PlotCMD(cmd_arr):
  # sort array and set the vectors using DataFrame
  df_cmd = DataFrame(cmd_arr, columns=['color','mag','count'] )#{'color':cmd_arr[:,0],'mag':cmd_arr[:,1],'count':cmd_arr[:,2]})  
  df_cmd = df_cmd.sort(['count'], ascending=[True])  # fixed on 28, March, 2014
  df_cmd = df_cmd.reset_index(drop=True) # fixed on 28, March, 2014
  color, mag, count = df_cmd['color'].values, df_cmd['mag'].values, df_cmd['count'].values
    
  # close prevoious plots
  plt.close('all')
  # setup the plot
  fig, ax = plt.subplots(1,1, figsize=(10,10))
  
  # define the colormap http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
  cmap = plt.cm.hsv
  # extract all colors from the .jet map
  cmaplist = [cmap(i) for i in range(cmap.N)]
  # force the first color entry to be grey
  cmaplist[0] = (.5,.5,.5,1.0)
  # create the new map
  cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
  
  # define the bins and normalize
  bounds = np.linspace(1,6,5+1)
  norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
  
  # define range of the CMD plot
  xmin, ymin = -0.4, min(mag)
  plt.axis([xmin, 1.6, ymin, 24])

  # make the scatter
  scat = ax.scatter(color,mag,c=count,s=5, marker= 'o', cmap=cmap, norm=norm, \
    edgecolors='gray')#, marker='o',s=np.random.randint(100,500,20)
  
  # create a second axes for the colorbar
  ax2 = fig.add_axes([0.9, 0.1, 0.01, 0.8]) # to set where the axis of colormap starts/ends
  cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', \
    ticks=bounds, boundaries=bounds, format='%1i') 

  #           Configure plot
  #
  ax.set_title('CMD')
  ax2.set_ylabel('Very custom cbar [-]', size=10)
  ax2.set_yticklabels(['1', '2', '3', '4', '5', '> 5']) # to set the names ox the ticks
  # Invert Y--axis
  ax.invert_yaxis()
  # Experimental tool to adjust plot to window
  # plt.tight_layout()
  #
  # modify space between subplots
  # plt.subplot_tool()
  plt.subplots_adjust(right=0.85)
  
  plt.show()
  #
  return
  #______

def SaveTab(aux_hx_arr, aux_N_div, aux_f_name):
  df_hbin = DataFrame( aux_hx_arr[:,:2+1], columns=['Color','Mag','N'] )
  df_bin = df_bin.sort(['N'], ascending=[True]) # fixed on 28, March, 2014
  df_bin = df_bin.reset_index(drop=True) # fixed on 28, March, 2014
  df_hbin.info()
  #
  # print the result on a table
  divisions = ( max(aux_hx_arr[:,0]) - min(aux_hx_arr[:,1]) ) / aux_N_div
  divisions = '{:.4f}'.format(divisions)
  aux_nm = aux_f_name[aux_f_name.find('.sub')-7:aux_f_name.find('.sub')]+'_'+divisions+'mag.hxb'
  df_hbin.to_csv(aux_nm, sep='\t', index=False)
  print '\n\nsaving the file... {0}'.format(aux_nm)
  #
  return
  #______


#==================================================================================
#
#          Corpus
#
if __name__ == "__main__":
  #
  # Call the hexbin
  thefolder = '/var/run/media/fj/Data/clusters-aldo/DatosHST/2nd_try/tables/' 
  #'/Users/fj/Work/clusters aldo/DatosHST/2nd_try/tables/sub_WorkingArea/'
  #f_name =  + 'NGC' + '0104' + '.sub'
  
  for (path, dirs, files) in os.walk(thefolder): #read ALL subfolders
    for each_file in files:   #walk over files  
      if ('.sub' in each_file) and ('NGC' in each_file):
	#
	# print the name of the file
	print '\n\t--->\tFILE: {0}\n'.format(each_file[each_file.find('.sub')-7:])
	#
	# set the complete path to file:
	f_name = path + each_file
	#
	# To control the while statement
	Control_continue = True
        #
	# first value to start
        N_div = 1000

	while Control_continue:
          hx_arr = ReadTab(f_name,N_div)
	  #
	  # Plot the resulted CMD
	  print '............shape of array to call Plot: {0}'.format(hx_arr[:,:2+1].shape)
	  PlotCMD(hx_arr[:,:2+1])#PlotCMD(camd_arr) arr-->color,mag,count
	  
	  #
	  # If plot looks good, then save it
	  SAVE_tbl = int( raw_input('\nis the hexbinning ok? (N={0}) (1/0): '.format(N_div)) ) 
	  #
	  #     when answer is '1'
	  if SAVE_tbl:
	    SaveTab(hx_arr, N_div, f_name)
	    #
	    Control_continue = False
	    #
	  # when answer is '0', ask for new N_div
	  if not SAVE_tbl:
	    N_div = int( raw_input('\n\tinsert the number of divisions for hxb (actual={0}) : '\
	    	.format(N_div)) )
	    print '\n...doing again\n'
 

else:
  print "main is being called form other script"


