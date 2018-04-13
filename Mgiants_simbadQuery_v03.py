# Script to query over SIMBAD, using a table from Hipparcos
# If use/modify, please refer to:
#                                                 Francisco Paz-Chinchon
#                                                 UFRN, Natal, Brazil
#                                                 francisco at dfte.ufrn.br

from time import time
from pandas import *
import numpy as np
from astroquery.simbad import Simbad
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table
from sys import exit
import sys, traceback
import gc

#
#           Function to adopt the format "HIP NNNNNN"
#
def Identifier(hip):
    hip = str(int(hip))
    control = True
    while control:
        if len(hip) < 6:
            hip = '0' + hip
        elif len(hip) == 6:
            control = False
        elif len(hip) > 6:
            print '\n ERROR!'; exit(0)
    hip = 'HIP ' + hip
    return hip

#
#           Absolute magnitude
#
def AbsMag(flux, parallax):
    print 'ORGANIZE!!!!!!!!!'
    plx = np.array(s_tab['PLX_VALUE'])
    # to ap/abs magnitude
    v_mag, k_mag = -2.5*np.log10(v_flux), -2.5*np.log10(k_flux)  
    print v_flux-k_flux
    V_mag, K_mag = v_mag+5.*np.log10(10.*plx), k_mag+5.*np.log10(10.*plx)
    return

# ............................................

# Start time
t_init = time()


#
#           Load table of M giants
#
M = read_csv('hipp_Mgiants.csv')
# _RAJ2000,_DEJ2000,HIP,Vmag,RA(ICRS),DE(ICRS),Plx,pmRA,pmDE,e_Plx,e_pmRA,e_pmDE,B-V,e_B-V,Period,rho,e_rho,Notes,SpType,r_SpType

M.info

# extract and set hipparcos identifiers in the Simbad way
aux_hip = map( Identifier, M['HIP'].values)

# set an astroquery-object to search in
myS = Simbad()
# configure it
myS.reset_votable_fields()
myS.add_votable_fields('id(HIP)', 'coo(GAL)', 'sptype', 'flux(B)', 'flux(V)', 'flux(K)', 'parallax', 'pm', 'pmra', 'pmdec')
myS.TIMEOUT = 120 
myS.ROW_LIMIT = 10000 #None 
# Output columns (from console):
#                            dtype = [('MAIN_ID', 'O'), ('RA', 'S13'), ('DEC', 'S13'), ('RA_PREC', '<i2'), ('DEC_PREC', '<i2'), ('COO_ERR_MAJA', '<f4'), ('COO_ERR_MINA', '<f4'), ('COO_ERR_ANGLE', '<i2'), ('COO_QUAL', 'S1'), ('COO_WAVELENGTH', 'S1'), ('COO_BIBCODE', 'O'), ('ID_HIP', 'O'), ('RA_GAL', 'S13'), ('DEC_GAL', 'S13'), ('SP_TYPE', 'O'), ('SP_QUAL', 'S1'), ('SP_NATURE', 'S1'), ('SP_BIBCODE', 'O'), ('FLUX_B', '<f4'), ('FLUX_V', '<f4'), ('FLUX_K', '<f4'), ('PLX_VALUE', '<f8'), ('PLX_PREC', '<i2'), ('PLX_ERROR', '<f4'), ('PLX_QUAL', 'S1'), ('PLX_BIBCODE', 'O'), ('PMRA', '<f8'), ('PMDEC', '<f8'), ('PMRA_2', '<f8'), ('PMDEC_2', '<f8')]

# We now must to search in the table for each M giant companions if there is a star that fit
# into the binary criteria
#
sel_all =  [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
sel_hp, sel_comp = [[],[],[],[],[],[],[],[],[],[],[],[],[]], [[],[],[],[],[],[],[],[],[],[],[],[],[]]
sel_except = []
for ind_hp, hp in enumerate(aux_hip):
    gc.collect()
    try:
        print "INDEX: {0}".format(ind_hp)
        # Do the search in a radius of 20 arcmin, for a single hipparcos M giant
        s_tab = myS.query_region(hp, radius=20. * u.arcmin)
        # Constraints (updated):
        #

        idx_vk = [n for n,m in enumerate( np.array(s_tab['ID_HIP']) ) if 'HIP '+str(int(hp[hp.find('HIP')+3:])) in m]
        # a priori, the M III must be the first element in the output table. This is only a caution
        # procedure.
        # control there is at least one coincidence
        if len(idx_vk) == 1:
            idx_vk = idx_vk[0]
            # as 'idx_vk' is a list, I'll use the first element, otherwise, the condition do not work:
            # [1.85] > 3.5 >>> FALSE
            # 1.82 > 3.5   >>> TRUE
            # When inside a list, comparison do not work.
        else:
            print '\t>>>{0} and its HIP identifier do not match'.format(hp)
            continue
        s_tab_Giant = s_tab[idx_vk]
        #print s_tab_Giant

        # 1) M giant must have V - K >= 3.5
        # Flux:
        V, K = s_tab_Giant['FLUX_V'], s_tab_Giant['FLUX_K'] #np.array(s_tab['FLUX_V'])[idx_vk], np.array(s_tab['FLUX_K'])[idx_vk]
        # if V-K >= 3.5, the loop continues, otherwise, breaks here and continue the next iteration
        if V-K < 3.5:
            print '\t>>>{0} has a value of V-K lower than 3.5'
            continue

        # 1.2) Possible companion must have < 3.5
        # flux of the companion
        v_flux, k_flux = np.array(s_tab['FLUX_V']), np.array(s_tab['FLUX_K'])
        idx_vk_comp = [n for n,m in enumerate(v_flux) if ( m - k_flux[n] ) < 3.5]
        # control 
        if len(idx_vk_comp) > 0:
            # replace the table by its sub version
            s_tab = s_tab[idx_vk_comp]
        else:
            print '\n\t>>>{0} companion has no (V-K) < 3.5'.format(hp)
            continue

        # 2) FGK
        # take only the spectral type column and transform to array. This way the masked values     
        # of the astropy-table are transformed to 'nan' in the array
        sptype = np.array(s_tab['SP_TYPE'])
        # search the index of the stars having spectral type F, G or K. Due to that,
        # there is no reason to exclude of the iteration the target hipparcos M giant, because
        # it do not obey the criteria of spectral type
        idx_sptype = [j for j,k in enumerate(sptype) if ( ('F' in k) or ('G' in k) or ('K' in k) )]
        # control there is at least one coincidence
        if len(idx_sptype) > 0:
            # replace the table by its sub version
            s_tab = s_tab[idx_sptype]
        else:
            print '\n\t>>>{0} has no near FGK'.format(hp)
            continue

        # 3) Parallax within 2 times the higher error
        ind_plx_sel = []
        error_M, plx_M =  s_tab_Giant['PLX_ERROR'], s_tab_Giant['PLX_VALUE']
        if np.isfinite( plx_M ):
            # list of indices where parallax is not missing, except the index of the M giant
            indx_plx = [i_cc for i_cc,cc in enumerate(s_tab['PLX_VALUE']) if np.isfinite( cc )]
            if len(indx_plx) > 0:
                TMP = s_tab[indx_plx]
                for j,k in enumerate( np.array(TMP['PLX_ERROR']) ):
                    two_err = np.array([ k, error_M ])
                    # drop NaN
                    two_err = two_err[ np.logical_not( np.isnan(two_err) ) ]
                    if (len(two_err) > 0) and (np.abs( plx_M-np.array(TMP['PLX_VALUE'])[j]  ) < 2.0*two_err.max()) and (np.abs( np.float(plx_M)-np.float(k) ) != 0.0):
                        # select indices
                        #print ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;', np.abs( np.float(plx_M)-np.array(TMP['PLX_VALUE'])[j] ), 2.0*two_err.max(), np.float(plx_M), np.array(TMP['PLX_VALUE'])[j]
                        ind_plx_sel.append(j)
        else:
            print '\n\t>>>{0} has no data of parallax for the M III'.format(hp)
            continue

        if len(ind_plx_sel) > 0:
            # print s_tab[ind_plx_sel]
            # replace the table by its sub version
            s_tab = s_tab[ind_plx_sel]
        else:
            # print error_M
            print '\n\t>>>{0} has no stars with similar parallax'.format(hp)
            continue 

        # Up to here, the selection has been made, if all criteria were matched.
        print "\n\t\t>>> {0} has {1} possible companions".format(hp, len(s_tab))
        #
        # Store in an array!
        # print s_tab_Giant['MAIN_ID'], len(s_tab_Giant['MAIN_ID']), type(s_tab_Giant['MAIN_ID'])
        # Resume all
        sel_all[0].append( hp ); sel_all[1].append( M['SpType'].values[ind_hp] )
        sel_all[2].append( s_tab_Giant['MAIN_ID'] ); sel_all[3].append( s_tab_Giant['SP_TYPE'] )
        sel_all[4].append( s_tab_Giant['FLUX_B'] ); sel_all[5].append( s_tab_Giant['FLUX_V'] )
        sel_all[6].append( s_tab_Giant['FLUX_K'] ); sel_all[7].append( 'A' )
        sel_all[8].append( s_tab_Giant['PLX_VALUE'] ); sel_all[9].append( s_tab_Giant['PMRA'] )
        sel_all[10].append( s_tab_Giant['PMDEC'] )
        sel_all[11].append( s_tab_Giant['RA'] ); sel_all[12].append( s_tab_Giant['DEC'] )
        sel_all[13].append( s_tab_Giant['RA_GAL'] ); sel_all[14].append( s_tab_Giant['DEC_GAL'] )
        sel_all[15].append( s_tab_Giant['ID_HIP'] )

        # M III Giant
        sel_hp[0].append( s_tab_Giant['MAIN_ID'] ); sel_hp[1].append( hp )
        sel_hp[2].append( s_tab_Giant['RA'] ); sel_hp[3].append( s_tab_Giant['DEC'] )
        sel_hp[4].append( s_tab_Giant['SP_TYPE'] ); sel_hp[5].append( s_tab_Giant['FLUX_B'] )
        sel_hp[6].append( s_tab_Giant['FLUX_V'] ); sel_hp[7].append( s_tab_Giant['FLUX_K'] )
        sel_hp[8].append( s_tab_Giant['PLX_VALUE'] ); sel_hp[9].append( s_tab_Giant['PMRA'] )
        sel_hp[10].append( s_tab_Giant['PMDEC'] )
        sel_hp[11].append( s_tab_Giant['PLX_ERROR'] ); sel_hp[12].append( s_tab_Giant['PLX_QUAL'] )
        #
        for ii,cc in enumerate(s_tab['MAIN_ID']):
            sel_comp[0].append(cc); sel_comp[1].append( s_tab['ID_HIP'][ii] )
            sel_comp[2].append( s_tab['RA'][ii] ); sel_comp[3].append( s_tab['DEC'][ii] )
            sel_comp[4].append( s_tab['SP_TYPE'][ii] ); sel_comp[5].append( s_tab['FLUX_B'][ii] )
            sel_comp[6].append( s_tab['FLUX_V'][ii] ); sel_comp[7].append( s_tab['FLUX_K'][ii] )
            sel_comp[8].append( s_tab['PLX_VALUE'][ii] ); sel_comp[9].append( s_tab['PMRA'][ii] )
            sel_comp[10].append( s_tab['PMDEC'][ii] )
            sel_comp[11].append( s_tab['PLX_ERROR'][ii] ); sel_comp[12].append( s_tab['PLX_QUAL'][ii] )
            #
            sel_all[0].append( hp ); sel_all[1].append( M['SpType'].values[ind_hp] )
            sel_all[2].append( s_tab['MAIN_ID'][ii] ); sel_all[3].append( s_tab['SP_TYPE'][ii] )
            sel_all[4].append( s_tab['FLUX_B'][ii] ); sel_all[5].append( s_tab['FLUX_V'][ii] )
            sel_all[6].append( s_tab['FLUX_K'][ii] ); sel_all[7].append( 'B' )
            sel_all[8].append( s_tab['PLX_VALUE'][ii] ); sel_all[9].append( s_tab['PMRA'][ii] )
            sel_all[10].append( s_tab['PMDEC'][ii] )
            sel_all[11].append( s_tab['RA'][ii] ); sel_all[12].append( s_tab['DEC'][ii] )
            sel_all[13].append( s_tab['RA_GAL'][ii] ); sel_all[14].append( s_tab['DEC_GAL'][ii] )
            sel_all[15].append( s_tab['ID_HIP'][ii] )
            #
    except:
        sel_except.append(hp)
        traceback.print_exc(file=sys.stdout)
        print '___________________INDEX OF CRASH: {0}'.format(ind_hp)
        pass


# Write results to a file

# exceptions
DataFrame({'HIP_main':sel_except}).to_csv('MIII_exceptions_v03.csv', index=False, header=True)

# companion selection
dictionary = {'01_id':sel_comp[0][:],'02_hip':sel_comp[1][:],'03_ra':sel_comp[2][:],'04_dec':sel_comp[3][:],'05_sptype':sel_comp[4][:],'06_B_flux':sel_comp[5][:],'07_V_flux':sel_comp[6][:],'08_K_flux':sel_comp[7][:],'09_plx':sel_comp[8][:],'10_pmra':sel_comp[9][:],'11_pmdec':sel_comp[10][:],'12_err_plx':sel_comp[11][:],'13_quality_plx':sel_comp[12][:]}
out_comp = DataFrame(dictionary)
#out_comp = out_comp.replace({'--', np.nan})
out_comp.to_csv('MIII_comp_v03.csv', index=False, header=True)

# M III selection
dictionary_2 = {'01_id':sel_hp[0][:],'02_hip':sel_hp[1][:],'03_ra':sel_hp[2][:],'04_dec':sel_hp[3][:],'05_sptype':sel_hp[4][:],'06_B_flux':sel_hp[5][:],'07_V_flux':sel_hp[6][:],'08_K_flux':sel_hp[7][:],'09_plx':sel_hp[8][:],'10_pmra':sel_hp[9][:],'11_pmdec':sel_hp[10][:],'12_err_plx':sel_hp[11][:],'13_quality_plx':sel_hp[12][:]}
out_hp = DataFrame(dictionary_2)
#out_hp = out_hp.replace({'--', np.nan})
out_hp.to_csv('MIII_giant_v03.csv', index=False, header=True)

# The resume of all
dictionary_3 = {'01_src_ID':sel_all[0][:],'02_src_sptype':sel_all[1][:],'03_xm_ID':sel_all[2][:],'04_xm_sptype':sel_all[3][:],'05_b_flux':sel_all[4][:],'06_v_flux':sel_all[5][:],'07_k_flux':sel_all[6][:],'08_flag':sel_all[7][:],'09_plx':sel_all[8][:],'10_pmra':sel_all[9][:],'11_pmdec':sel_all[10][:],'12_ra':sel_all[11][:],'13_dec':sel_all[12][:],'14_ra_gal':sel_all[13][:],'15_dec_gal':sel_all[14][:],'16_xm_hip_ID':sel_all[15][:]}
out_all = DataFrame(dictionary_3)
out_all.to_csv('MIII_both_v03.csv', index=False, header=True)

print '\n\n\t\t>.>.>. Total elpased time: {0:.4f} minutes'.format((time()-t_init)/60.)






exit(0)
s_tab = myS.query_region(aux_hip[0], radius=20. * u.arcmin) #<class 'astropy.table.table.Table'>

print myS.get_votable_fields()
print myS.get_field_description('mk')
print s_tab['PMRA']
print s_tab.colnames
s_tab = s_tab[np.array([0,1,2])]
print len(s_tab), s_tab
exit(0)
a = s_tab['PMRA']
print type(a)
print a[:20]

exit(0)






#           S I M B A D  search
#           Examples of search:
# result_table = Simbad.query_region("m81", radius=0.1 * u.arcmin)
# result_table = Simbad.query_region(coord.SkyCoord(31.0087, 14.0627,\
#                                        unit=(u.deg, u.deg), frame='galactic'),\
#                                        radius='0d0m2s')
# More info: http://astroquery.readthedocs.org/en/latest/simbad/simbad.html
#
# Some functions:
#
# Simbad.list_votable_fields() 
# customSimbad = Simbad()
# customSimbad.get_field_description('mk')
# customSimbad.get_votable_fields()
# customSimbad.add_votable_fields('mk', 'rot', 'bibcodelist(1800-2014)')
# customSimbad.remove_votable_fields('mk', 'coordinates')
# customSimbad.reset_votable_fields()
# Simbad.TIMEOUT = 60 # sets the timeout to 60s
# Simbad.ROW_LIMIT = 20 # now any query fetches at most 15 rows
