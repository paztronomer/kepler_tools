#		FIT & INTERESTING POINTS
#
# Programm to:
#             fit cmd and find interesting points
# If use/modify, please reference: F. Paz-Chinchon
#                                                francisco@dfte.ufrn.br
#                                                DFTE, UFRN, Brazil.

import os
from pandas import *
import numpy as np
import pylab as plab
import scipy
import itertools                                       # to double iterate 
from sys import exit                                   # use of exit(0)
#import rpy2.robjects as R

# _____________________________________________________________________
# _____________________________________________________________________

#                  F U N C T I O N S  &  P R O C
#                      various tasks, easy life

#
# ::::::::::::::::::: Find Derivates and Roots in some interval
#
# df_dx = f.deriv()
# roots_g = np.roots(g)
#
def Inflexion(polinome,x_i,x_f): # Call: polynomial,x_init,x_final
    dx_P = polinome.deriv() #object, polyn class
    shot = 0
    #print 'initial: {0}  final: {1}'.format(x_i, x_f)
    for r in np.roots(dx_P): #array of roots
        if (r >= x_i) and (r <= x_f):
            root = r
            shot += 1
    if shot != 1: 
        print 'ERROR: no single root in the range'
        exit(0)
    polinome= dx_P= None #cleanning variables
    return root

#
# ::::::::::::::::::: Find x value, given y value in a polinomial
#
# I need to inverse evaluate
# P(x) to know x_i, given a P(x_i)
#
def xInPol(poly,y_value,x_init,x_fin):
    preciss, x_value, steps = 0.001, None, 1000
    print preciss
    while not x_value:
        print steps
        for it in np.linspace(x_init,x_fin,num=steps):
            if abs(poly(it) - y_value) <= preciss:
                x_value = it
        steps = 2*steps
    return x_value

# _____________________________________________________________________
# _____________________________________________________________________


#                       Calls &
#                  Variables to run

if __name__ == "__main__":

    path_hxb = '.' #walk through path
    path_split = 'split_7099.txt' # table containing critical points
    
    # Columns: NGC  cut.mag_TO  cut.colr_TO  split.colr_TO.sRGB  cut.mag_limit
    df_spl = read_table(path_split, sep='\s*', names=['ngc','toff_mag','toff_colr','split_colr','cut_mag'], skiprows=1)
    
    # the names of hexbin tables are of the form: NGC6441_0.0142mag.hxb
    for (path,dirs,files) in os.walk(path_hxb):
        for file in files:
            if ('NGC' in file) and ('mag.hxb' in file): 
                
                # I will charge the tables, using pandas, in order to avoid
                # self-functions. Vectorization seems to be more efficient
                # in this scenario.
                #
                #          C H A R G E   T A B L E S
                #
                fname = path + '/' + file
                # hexbin tables format '*.hxb': | Color | Mag | N |
                df_hxb = read_table(fname, sep='\s*')
                df_hxb_original = df_hxb
                
                # index of DF of the split table, when matchs the NGC of the hxb file
                index_ngc = df_spl['ngc'] == int( file[file.find('C')+1:file.find('C')+5] )
                #df_spl[ind1]['cut_mag']
                
                # sort df_hxb by N
                df_hxb = df_hxb.sort(['N'], ascending=[0])
                
                # cut by N, leaving only bins with some value or higher
                df_hxb = df_hxb[ df_hxb['N'] > 1 ]
                # df_hxb.info()
                
                # define main sequence and cut magnitude about 1 mag below TOff, 
                # using values of df_spl[]
                df_hxbMS = df_hxb[ df_hxb['Mag'] > float(df_spl[index_ngc]['toff_mag']) ]
                df_hxbMS = df_hxbMS[ df_hxbMS['Mag'] <= float(df_spl[index_ngc]['cut_mag']) ]
                # df_hxbMS.info()
                #
                df_hxb = df_hxb[ df_hxb['Mag'] <= float(df_spl[index_ngc]['toff_mag']) ]
                # df_hxb.info()

                # cut color near TOff
                df_hxb = df_hxb[ df_hxb['Color'] >= float(df_spl[index_ngc]['toff_colr']) ]
                # df_hxb.info()

                # separe upper CMD, by color
                df_hxbLeft = df_hxb[ df_hxb['Color'] <= float(df_spl[index_ngc]['split_colr']) ]
                df_hxbRight = df_hxb[ df_hxb['Color'] > float(df_spl[index_ngc]['split_colr']) ]
                # and sort by magnitude
                df_hxbLeft = df_hxbLeft.sort(['Mag'], ascending=[0])
                df_hxbRight = df_hxbRight.sort(['Mag'], ascending=[0])
                # df_hxbLeft.info()
                # df_hxbRight.info()
                
                # inside the loop I will perform the fit of the polinomials
                # the fit is divided in 3 sections: TOff, MS, sRGB
                #
                #          F I T   C M D   B Y   P A R T S
                #
                # http://docs.scipy.org/doc/scipy/reference/interpolate.html
                # Works fine with sRGB and TOff
                # Simplest polynomial
                # http://docs.scipy.org/doc/numpy/reference/generated/numpy.poly1d.html

                # 1) TOff
                degree = 3
                Fy_toff = np.poly1d( np.polyfit(df_hxbLeft['Mag'].values, df_hxbLeft['Color'].values, \
                                                    deg=degree, w=df_hxbLeft['N'].values) , variable='y') # x=f(y)
                # 2) sRGB
                #         Note: 
                #         Degree 5 fit well the desired subRGB region.
                Fx_rgb = np.poly1d( np.polyfit(df_hxbRight['Color'].values, df_hxbRight['Mag'].values, \
                                                   deg=degree+2, w=df_hxbRight['N'].values), variable='x') # y=f(x)
                # 3) MS
                Fx_ms = np.poly1d( np.polyfit(df_hxbMS['Color'].values, df_hxbMS['Mag'].values, \
                                                  deg=3,w=df_hxbMS['N'].values) , variable='x') # y=f(x)
                
                # now, evaluate the polynomials in order to get array-like values
                Ev_toff = Fy_toff( df_hxbLeft['Mag'].values )
                Ev_rgb = Fx_rgb( df_hxbRight['Color'].values )
                Ev_ms = Fx_ms( df_hxbMS['Color'].values )
                
                # plot the CMD with the fitted polynomials
                DO_plot = False
                if DO_plot:
                    plab.close('all')
                    # TOff + sRGB
                    plab.plot( df_hxbLeft['Color'].values, df_hxbLeft['Mag'].values, 'bo', \
                                  df_hxbRight['Color'].values, df_hxbRight['Mag'].values, 'go', alpha=0.75)
                    plab.plot( Ev_toff, df_hxbLeft['Mag'].values, 'y-', df_hxbRight['Color'].values, Ev_rgb, 'r-')
                    plab.gca().invert_yaxis()
                    plab.show()
                    #
                    # MS
                    plab.plot( df_hxbMS['Color'].values, df_hxbMS['Mag'].values, 'yo', alpha=0.75)
                    plab.plot( df_hxbMS['Color'].values, Ev_ms, 'o')
                    plab.gca().invert_yaxis()
                    plab.show()
                    
                # now, using the polynomial fits, evaluate critical points.
                # the most important point is the TOff, because it will be 
                # used as start point.
                #
                #          D E F I N E   C R I T I C A L   P O I N T S
                #

                # 1) Turn Off (color:c, magnitude:m)
                #     calls the Inflexion function, which determines when the derivate becomes zero
                m_1 = Inflexion( Fy_toff, min(df_hxbLeft['Mag'].values), \
                                     max(df_hxbLeft['Mag'].values) ) #TOff is well determinated! (inflexion point)
                #     then, evaluate the magnitude in the polinome
                c_1 = Fy_toff(m_1) 

                # 2) sub RGB (color:c, magnitude:m)
                #     first, define the offset by hand, because up to now no analytical method worked
                DeltaMag = -1
                m_2 = m_1 + DeltaMag #might justify it with some references
                #     calls the function xInPol, which get 'x', given 'f(x)'
                #     precission=0.001
                c_2 = xInPol( Fx_rgb, m_2, min(df_hxbRight['Color'].values), \
                                  max(df_hxbRight['Color'].values)) 

                # 3) Main Sequence
                #     finally, to set the critical points in the MS, I used the same color than in sRGB,
                #     applied on the MS polynomial fit
                c_0 = c_2
                m_0 = Fx_ms(c_0)

                # plot the points and the cmd
                DO_plot2 = True
                if DO_plot2:
                    plab.close('all')
                    # TOff + sRGB
                    plab.plot( df_hxbLeft['Color'].values, df_hxbLeft['Mag'].values, 'bo', \
                                   df_hxbRight['Color'].values, df_hxbRight['Mag'].values, 'go', alpha=0.75, markersize=2)
                    plab.plot( Ev_toff, df_hxbLeft['Mag'].values, 'y-', df_hxbRight['Color'].values, Ev_rgb, 'r-', markersize=3)
                    # MS
                    plab.plot( df_hxbMS['Color'].values, df_hxbMS['Mag'].values, 'yo', alpha=0.75, markersize=2)
                    plab.plot( df_hxbMS['Color'].values, Ev_ms, 'b.', markersize=3)
                    plab.title( file[ : file.find('C')+5] )
                    # critical points
                    plab.plot( c_0,m_0,'ws', markersize=5)
                    plab.plot( c_1,m_1, 'ws', markersize=5)
                    plab.plot( c_2,m_2, 'ws', markersize=5)
                    plab.gca().invert_yaxis()
                    plab.show()

                # 4) Transform into relative coordinates
                #     using the TOff as origin (0,0), I set an relative-to-TOff framework
                c_0, m_0 = c_0-c_1, m_0-m_1 # >0, positive
                c_2, m_2 = c_2-c_1, m_2-m_1 # <0, negative
                c_1, m_1 = 0.0, 0.0                # (0,0), origin


else:
  print 'this module has been called by other script'


# Pending!!!
# 1) Save the critical points into a DF and write out
# 2)  Make a table with input critical points











#__________________

# ____________________________________________________________________________________
# ____________________________________________________________________________________
# ____________________________________________________________________________________


# TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#                                     TESTING AREA
# _LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
# 1st: determine inflexion points and save it --- OK for LEFT side!
#                                             --- right side have no inflexion point. Set BY HAND. --- DONE!
#                                             --- MS
# 2nd: do that for all CMD
# 3rd: compare and match theorical tracks with CMD
# 4th: erase unused variables and run for all






#dxF_r = Fx_r.deriv()
#print np.roots(dxF_r)
#print Inflexion(Fx_r,min(x_right),max(x_right)) 

#print np.random.normal(scale=0.1)

# TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
# TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#        _____________________
#
#            N  O  T  E  S
#
# R.r.median(R.IntVector([1,2,3,4]))[0]
#
# add noise
#x+= random.normal(scale=0.1, size=x.shape)
#y+= random.normal(scale=0.1, size=y.shape)
#plab.scatter(x,y)
#plab.show()


# Methods to open file:
#with open('somefile.txt', 'r') as FILE:
#    for i in FILE:


#import fileinput
#for i in fileinput.input('somefile.txt'):
#    print "hola"
# Easy to plot
#plab.xlim(min(time), max(time)) #set axis extremes
#            plab.ylim(min(flux)-np.std(flux), max(flux)+np.std(flux)) #set axis extremes, using StdDev
            # scatter plot
#            plab.scatter(time, flux, s=5, color="green", linewidth=0.5, label= file[:file.find(".tsv")])
#            plab.legend(loc='upper left') #locate label in figure
#save figure using 72 dots per inch
#            fname_plot = 'pl.'+file[file.find("Q"):file.find("i")]+'png'
#            plab.savefig(fname_plot, dpi=72)
# showing plot on 
#plab.show()
