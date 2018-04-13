''' 
Perform aperture and PSF photometry on Kepler Target Pixel files

'''

from __future__ import division
from scipy import stats #multivariate_normal
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import time
import os                        # walk through directories
import glob
import astropy.io.fits as fits
from progressbar import ProgressBar



def Read_DataCube(path_targ):
    ''' open FITS file using astropy/pyfits
    mmap=True helps when reading heavy files
    '''
    hdu=fits.open(path_targ, mmap=True)
    # hdu.info() # print info 

    ''' FITS file are divided in 3-main bodies:
    0: PRIMARY     PrimaryHDU
    1: TARGETTABLES  BinTableHDU
    2: APERTURE    ImageHDU
    the 2nd body, BinTables stores the aperture images
    '''
    
    # Store BinTable in a numpy array
    hdu1 = hdu[1].data
    # print type(hdu1)

    '''
    Header Info
    '''
    KIC = str( hdu[0].header['KEPLERID'] )
    skygroup = str( hdu[0].header['SKYGROUP'] )
    Y_mask = str( hdu[2].header['CRVAL2P'] )
    X_mask = str( hdu[2].header['CRVAL1P'] )
    cadence = hdu[0].header['OBSMODE']
    # reference coordinates of object
    ref_x = hdu[1].header['1CRPX4']
    ref_y = hdu[1].header['2CRPX4']
    # Observation time in BJD-BJDREF
    OBS_START = hdu[1].header['TSTART']
    OBS_END = hdu[1].header['TSTOP']
    # Midpoint of first cadence in MJD
    LC_START = hdu[1].header['TSTART']
    LC_END = hdu[1].header['TSTART']
    # Frame times
    T_INIT = hdu[1].header['INT_TIME'] # photon accumul. time per frame; [s] 
    T_READ = hdu[1].header['READTIME'] # readout time per frame; [s]
    T_FRM = hdu[1].header['FRAMETIM'] # frame time = T_INIT + T_READ
    # Gain and spatial separation
    GAIN = hdu[1].header['GAIN'] # ccd gain factor; [e- / ADU]
    DT_BIN = hdu[1].header['TIMEDEL'] # time resolution of data; [d]
    # Cadence: long=270, short=9
    N_FRM = hdu[1].header['NUM_FRM'] # number of frames per time stamp

    # Close FITS instance
    hdu.close()
    
    ''' Get the datacube
    hdu[***][3] --> raw counts
    hdu[***][4] --> (calibrated) flux counts
    get also the flux (with the corrections applied) matrix
    '''
    # Determine number of layers of the datacube:
    # Nlayers = len(hdu1[:]) or more elegant
    Nlayers = hdu1.shape[0]
    
    # Show header
    # time | timecorr | cadenceno | raw_cnts | flux | flux_err | flux_bkg |
    # | flux_bkg_err | cosmic_rays | quality | pos_corr1 | pos_corr2 
    # print hdu[1].header

    '''Stack all layers of the raw and corrected photometry into a 3D matrix.
    I'll save the 1st element to define the list-matrix
    '''
    # Htime will harbor time, while Hraw and Hcalib will harbor raw and calibrated
    # data, respectively.
    Htime = []
    Htime.append(hdu1[0][0])
    Hraw = hdu1[0][3]
    Hcalib = hdu1[0][4]
    # now the filling loop... starting on the 2nd element
    for k in range(0+1,Nlayers,1):   
        Hraw=np.dstack((Hraw, hdu1[k][3] ))
        Hcalib=np.dstack((Hcalib, hdu1[k][4] ))
        Htime.append(hdu1[k][0])
        
    # Hflux and Hraw arrays are filled with the star's measures
    # BUT we must keep in mind that these values are the sum over
    # - 270 frames (of 6.02 sec each) for LONG cadence
    # &
    # - 9 frames fort SHORT cadence
    #
    # To cut above the saturation limit (raw) DN~10,000, we must 
    # normalize. Also NaN or "-1" values mut be deprecated
    # ---> normalized_Hraw = Hraw/N_FRM
    ''' Transform flux arrays from int to float
    '''
    return Hraw.astype(float),Hcalib.astype(float),np.array(Htime),N_FRM



def NearNeighbors_Pix(arr,xx,yy,n):
    ''' Given a 2D-array, returns an nxn (typically n=3) array whose 
    "center" element is arr[xx,yy]
    - np.roll(arr, shift, axis) :roll array elements along a given 
    axis. Shift=number of places by which elements are shifted
    coord_final+coord_inic=1, because the element of interest becames 
    the element (1,1) of the rolled array

    Taken from: http://stackoverflow.com/questions/4148292/how-do-i-select-a-window-from-a-numpy-array-with-periodic-boundary-conditions
    '''
    # Return a nxn stamp
    arr=np.roll( np.roll(arr,shift=1-xx,axis=0), shift=1-yy,axis=1)
    return arr[:n,:n]


    
def Gauss_Profile(im_star,sigma_psf):
    from scipy import ndimage
    ''' Set by hand and carefully the sigma (width) of the gaussian 
    profile to avoid over/under samplig
    '''
    img_gauss = ndimage.filters.gaussian_filter(im_star,sigma=sigma_psf,order=0,mode='nearest')
    return img_gauss



def star_CM(img_soft):
    ''' Called with the gaussian profile
    Mask NaN values
    '''
    img_ma = np.ma.array(img_soft,mask=np.isnan(img_soft))
    # Center of mass
    try:
        from scipy.ndimage.measurements import center_of_mass
        cmass = center_of_mass(img_ma)
        # Maximum of Measurements
        # from scipy.ndimage import measurements
        # max_coord=measurements.maximum_position(img_ma)
        # print cmass,max_coord
        return cmass
    except:
        print '\n\t ### ERROR: center of mass calculation error. \
        No Flux measurement will be performed on this layer'
        return False



def Detect_Outliers(imlayer):
    ''' Estimate Iglewicz and Hoaglin criteria for outlier 
    http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
    Formula:
    Z=0.6745(x_i - median(x)) / MAD
    if abs(z) > 3.5, x_i is a potential outlier

    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect 
    and Handle Outliers", The ASQC Basic References in Quality Control: 
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    '''
    from statsmodels.robust import scale # Alternative way to estimate MAD
    from scipy.stats import norm

    # Percent point function (inverse of cdf -- percentiles) of a normal
    # continous random variable
    cte = norm.ppf(0.75) # aprox 0.6745
    '''
    Flatten the image, to estimate median
    Exclude '-1' values
    '''
    flat_im = imlayer.ravel()
    flat_im = flat_im[flat_im!=-1]
    
    # MAD
    MAD = np.median( np.abs(flat_im-np.median(flat_im)) )
    # alternative: scale.mad(flat_im, c=1, axis=0, center=np.median)
    
    # the  Iglewicz and Hoaglin scorer. If absolute value is
    # higher than 3.5, it is a possible outlier
    Zscore = cte*(flat_im-np.median(flat_im))/MAD
    Zscore = np.abs(Zscore)
    
    ''' 
    search for outliers and if present, replace by the median of neighbors
    '''
    counter_rem = 0
    if len(Zscore[Zscore>3.5])>0:
        for k in range(0,imlayer.shape[0]):
            for m in range(0,imlayer.shape[1]):
                ''' Try 1: If pixel has -1 value, then do not change its value
                '''
                if (np.abs( cte*(imlayer[k,m]-np.median(flat_im))/MAD ) >3.5) and (imlayer[k,m]!=-1):
                    # print '___values of possible outliers (noisy star profile, coords [{1},{2}]): {0:.4f} '.format\
                    #     (imlayer[k,m],k,m)
                    ''' 
                    if outlier coordinate is on the edges, then matrix is replicated
                    function NearNeighbors_Pix() creates an 3x3 array on which the 
                    element located at (1,1) is our pixel of interest
                    '''
                    neighbor = NearNeighbors_Pix(imlayer,k,m,n=3)
                    # set the outlier as NaN, located at [1,1] of the sub-matrix
                    neighbor[1,1] = np.nan
                    '''
                    then, estimate median, EXCLUDING the NaN value
                    (because numpy mess up with NaNs) and set this value
                    to the outlier
                    '''
                    flat_neighbor = neighbor.ravel(); flat_neighbor = flat_neighbor[~np.isnan(flat_neighbor)]
                    # replace the outlier 
                    imlayer[k,m] = np.median(flat_neighbor)
                    
                    # print '___replaced value: {0:.4f} '.format(imlayer[k,m])
                    counter_rem+=1
        print '\nDetect outliers remotion. \n\tPossible={0} , \n\tReplaced={1} ,\n\tTotal (includ. -1)={2}'.format(len(flat_im),counter_rem,len(imlayer.ravel()))
        return imlayer
    else:
        return imlayer



def Calc_Phot(image,xc,yc,radius):
    ''' aperture or PSF measurement of flux
    
    for aperture photometry integrate inside radius
    http://stackoverflow.com/questions/8647024/how-to-apply-a-disc-shaped-mask-to-a-numpy-array
     
    define two 1D arrays based in centroid of distribution (xc,yc) where 
    each one have values: -xc...(Nx-xc)-1  and  -yc...(Ny-yc)-1 
    Nx,Ny being number of elements on each axis
    
    np.ogrid returns meshgrid arrays of only 1 dimension
    '''
    cols,rows = np.ogrid[-xc:image.shape[0]-xc, -yc:image.shape[1]-yc]
    mask_r = rows**2 + cols**2 <= radius**2.   
    ''' 
    then apply this mask composed of True/False to the 2d image and get a submatrix. 
    Then add over 2 axis
    '''
    return np.sum(image[mask_r])

AQUI QUEDE!

def FT_Gaussian(im,sigma_gauss):
    # Fourier Transform of a Gaussian
    from scipy.ndimage.fourier import fourier_gaussian
    return fourier_gaussian(im, sigma=sigma_gauss)



##
## Following function main() not used
##
def main():
    print '\n \t starting, {0}\n\n'.format(time.ctime())
    #
    # ===================
    # S I M U L A T E 
    # S T A R
    # ===================
    #
    ''' call function to create fake star with and w/o noise
    '''
    fk_star,ns_star,x,y=create_Star()
    '''
    IMPORTANT: through the below code, ns_star will be affected and will change
    its profile. 
    '''
    ''' test plot for:
    ideal stellar profile: contour map + 3D profile
    noisy stellar profile: contour map + 3D profile
    '''
    if True:
        plt.close('all')
        fig=plt.figure(figsize=(10,5))
        ax1=fig.add_subplot(221)
        ax2=fig.add_subplot(222,projection='3d')
        ax3=fig.add_subplot(223)
        ax4=fig.add_subplot(224,projection='3d')
        #
        ax1.contourf(x,y,fk_star,cmap=cm.coolwarm)
        ax3.contourf(x,y,ns_star,cmap=cm.coolwarm)
        #
        surf_fk=ax2.plot_surface(x,y,fk_star, rstride=1, cstride=1, cmap=cm.coolwarm,
                                 linewidth=0)
        surf_ns=ax4.plot_surface(x,y,ns_star, rstride=1, cstride=1, cmap=cm.coolwarm,
                                 linewidth=0)
        ''' colorbar
        '''
        #fig.colorbar(surf, shrink=0.5, aspect=18)
        #fig.colorbar(surf_ns, shrink=0.5, aspect=18)
        #
        plt.show()

    #
    # ==================
    # VARIOUS 
    # CALCULATIONS
    # ==================
    #
    ''' See
    http://stackoverflow.com/questions/14765891/image-smoothing-in-python
    lets set as image the values for the noisy star, correct outliers, and then apply a gaussian filter 
    '''
    img_corr=Detect_Outliers(ns_star)
    
    print '\n--->>> PROBLEM: when correcting for outliers, the initial noisy profile is changed, being the same as the corrected one, as can be seen from its maximum value before/after ratio={0}'.format(ns_star.max()/img_corr.max())

    '''gaussian convolution
    there is a parameter to control outliers influence??
    '''
    guess_sigma=0.75
    img_gauss=Gauss_Profile(im_star=img_corr,sigma_psf=guess_sigma)
    print '\n======\nImportant issue: sigma of the PSF will highly influence \nthe goodness of fit between PSF and stellar profile. \nOver/under-sample will depend of it. So use the plots \nof difference between profiles to make your choice.\n\tCurrent sigma = {0}  \n======'.format(guess_sigma)

    # Center of mass for noise-free profile and gaussian profile
    CM,coord_peak=star_CM(img_gauss)
    CM_nsfree,coord_peak_nsfree=star_CM(fk_star)
    print '\n\tCM of gaussian profile {0},\n\tand the maximum of distribution {1}'.format(CM,coord_peak)
    print '\tCM of noise-free stellar profile {0},\n\tand the maximum of distribution {1}\n'.format(CM_nsfree,coord_peak_nsfree)

    # Multiplication of the noisy profile with the Fourier Transform of a Gaussian
    FT_star=FT_Gaussian(img_corr, sigma_gauss=1.)

    ''' Calculate aperture and PSF flux
    '''
    # center coordinates for aperture and PSF
    x_ap,y_ap=int(round(coord_peak[0])),int(round(coord_peak[1]))
    x_psf,y_psf=CM[0],CM[1]
    #
    ap_cal=Calc_Phot(img_corr,x_ap,y_ap,radius=4.)
    psf_cal=Calc_Phot(img_gauss,x_psf,y_psf,radius=4.)
    print '\n\tPhotometry for PSF:{0} and for Aperture:{1}\n'.format(psf_cal, ap_cal)


    #
    # ===============
    # P L O T S
    # ===============
    #
    plotFig=False

    ''' test plot for:
    noisy star contour plot
    gaussian filter to noisy star contour plot + center of mass of gaussian profile AND original star
    original noise-free star contour plot
    noisy star multiplied by Fourier transform of a gaussian
    '''
    if plotFig:
        plt.close('all')
        fig=plt.figure()
        ax1=fig.add_subplot(221)
        ax2=fig.add_subplot(222)
        ax3=fig.add_subplot(223)
        ax4=fig.add_subplot(224)
        #
        ax1.imshow(fk_star); ax1.set_title('original star',fontsize=15)
        ax1.plot(CM_nsfree[1],CM_nsfree[0],'wx',ms=50,alpha=1,lw=5)
        ax2.imshow(img_corr); ax2.set_title('noisy star',fontsize=15)
        ax3.imshow(img_gauss); ax3.set_title('gaussian filter noisy star',fontsize=15)
        ax3.plot(CM[1],CM[0],'w+',ms=50,alpha=1,lw=5)
        ax3.plot(CM_nsfree[1],CM_nsfree[0],'wx',ms=50,alpha=1,lw=5)
        ax4.imshow(FT_star); ax4.set_title('multiplied by FT of gaussian',fontsize=15)
        #
        ax1.set_xlim([0,fk_star.shape[1]-1]); ax1.set_ylim([fk_star.shape[0]-1,0])
        ax3.set_xlim([0,img_gauss.shape[1]-1]); ax3.set_ylim([img_gauss.shape[0]-1,0])
        #
        plt.show()
    
    ''' Checking plots: 3D profiles
    These plots allow us to see if profile fits enough well to data, without oversampling.
    Plots:
    3D noise-free simulated star
    3D noisy simulated star (no outliers are present)
    3D PSF & (PSF-noise free star)
    3D PSF & (PSF-noisy star)
    '''
    if plotFig:
        fig=plt.figure(figsize=(12,5))
        ax1=fig.add_subplot(141,projection='3d')
        ax2=fig.add_subplot(142,projection='3d')
        ax3=fig.add_subplot(143,projection='3d')
        ax4=fig.add_subplot(144,projection='3d')
        # noise free star
        ax1.plot_surface(x,y,fk_star, rstride=1, cstride=1, cmap=cm.coolwarm,
                         linewidth=0); ax1.set_title('Noise free star',fontsize=15)
        # noisy star
        ax2.plot_surface(x,y,img_corr, rstride=1, cstride=1, cmap=cm.coolwarm,
                         linewidth=0); ax2.set_title('Noisy star',fontsize=15)
        # gaussian PSF + gaussian PSF minus fake (noise free) star
        ax3.plot_surface(x,y,img_gauss-fk_star, rstride=1, cstride=1, cmap=cm.gray,
                         linewidth=0)
        ax3.plot_surface(x,y,img_gauss, rstride=1, cstride=1, cmap=cm.coolwarm,
                         linewidth=0,alpha=0.8); ax3.set_title('Gaussian PSF +\nPSF$-$noise-free star',fontsize=15)
        # gaussian PSF + gaussian PSF minus noisy star
        ax4.plot_surface(x,y,img_gauss-img_corr, rstride=1, cstride=1, cmap=cm.gray,
                         linewidth=0)
        ax4.plot_surface(x,y,img_gauss, rstride=1, cstride=1, cmap=cm.coolwarm,
                         linewidth=0,alpha=0.8); ax4.set_title('Gaussian PSF +\nPSF$-$noisy star',fontsize=15)
        ax1.view_init(elev=10., azim=-65.)
        ax2.view_init(elev=10., azim=-65.)
        ax3.view_init(elev=10., azim=-65.)
        ax4.view_init(elev=10., azim=-65.)
        #fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()
    


if __name__=='__main__':
    print '\n\tStarting... {0}'.format(time.ctime())
    t_init=time.clock()
    path_files='targ/'

    #_______________________________________________________
    #_______________________________________________________
    #_______________________________________________________
    #_______________________________________________________
    ''' FITS to numpy format
    '''
    if False:
        for filename in glob.glob(os.path.join(path_files, '*.fits')):
            aux_fn=filename[filename.rfind('/')+1:]
            print '\n\t___file: {0} [t_proc={1}]'.format(aux_fn,time.clock()-t_init)
            #
            # Read and save data
            dc_raw,dc_cal,dc_time,Nframes=Read_DataCube(filename)
            print '\t___[t_proc={0}]'.format(time.clock()-t_init)
            print '\tNumber of frames per bin={0}'.format(Nframes)
            aux_name=aux_fn[:aux_fn.find('-')+7].replace('-','_')
            #
            # Remove outliers by median
            probar=ProgressBar()
            dc_rawAux=Detect_Outliers(dc_raw[:,:,0])
            for dd in probar( range(0,dc_raw.shape[2]) ):
                dc_rawAux=np.dstack(( dc_rawAux,Detect_Outliers(dc_raw[:,:,dd]) ))
            #
            # Write out
            np.save('r_'+aux_name+'.npy',dc_raw)
            np.save('ro_'+aux_name+'.npy',dc_rawAux)
            np.save('c_'+aux_name+'.npy',dc_cal)
            np.save('t_'+aux_name+'.npy',dc_time)

    #_______________________________________________________
    #_______________________________________________________
    #_______________________________________________________
    #_______________________________________________________
    ''' Photometry calculation
    REMEMBER: 2D array posses NaN values which mess with calculations
    '''
    if False:
        ''' Lets use calibrated datacube
        '''
        path_npy='cygAB_npy/'
        for filename in glob.glob(os.path.join(path_npy, 'c_*')):
            dcube=np.load(filename)
        
            t_init, clock_init=time.time(), time.clock()
            # Load array of time
            aux_name=path_npy+'t_'+filename[filename.find('c_')+2:]
            jd_arr=np.load(aux_name)

            # check if there are saturated pixels, assuming all 
            # datacubes being long_cadence 
            if np.any(dcube[~np.isnan(dcube)]/270. >10000.): 
                num_satur= dcube[dcube[~np.isnan(dcube)]/270. >10000.].shape[0]
                print '\t {0} has saturated pixels, N={1}'.format(filename,num_satur )

            # Correct for outliers
            # dcube_2 will be the array corrected for outliers
            print '\n\t___Outliers detection+treatment.'
            pbar=ProgressBar()
            dcube_2=Detect_Outliers(dcube[:,:,0])
            ''' Note we start from the second index, '1'
            '''
            for dd in pbar( range(0+1,dcube.shape[2]) ):
                dcube_2=np.dstack(( dcube_2,Detect_Outliers(dcube[:,:,dd]) )) 
            
            # Set variables for photometry with different aperture radius
            print '\n\t___Calculating photometry values on a sub-image.'
            v_ap=[[],[],[],[]]
            v_psf=[[],[],[],[]]
            pbar2=ProgressBar()
            # time_cal will harbor indices of flux measurements
            time_cal=[]
            print dcube_2.shape,'\n\n'
            for dd in pbar2(xrange(0,dcube_2.shape[2])):
                # Crop image, substracting 4 pix from each border
                dcube_3=dcube_2[4:dcube_2.shape[0]-5,4:dcube_2.shape[1]-4,dd]
                # Gaussian profile
                guess_sigma=1.#0.75
                img_gauss=Gauss_Profile(im_star=dcube_3,sigma_psf=guess_sigma)
                #
                # Center of mass for noise-free profile and gaussian profile
                CM=star_CM(img_gauss)
                if CM !=-1:
                    time_cal.append(dd)
                    x_cm,y_cm=int(round(CM[0])),int(round(CM[1]))
                    '''==========Test Plot. Uncomment if needed============='''
                    # # test plot
                    # x, y=np.mgrid[0:dcube_3.shape[0]:1, 0:dcube_3.shape[1]:1]
                    # centers=np.dstack((x, y))
                    
                    # plt.close('all')
                    # fig=plt.figure(figsize=(10,5))
                    # ax1=fig.add_subplot(221)
                    # ax2=fig.add_subplot(222,projection='3d')
                    # ax3=fig.add_subplot(223)
                    # ax4=fig.add_subplot(224,projection='3d')
                    # #
                    # ax1.contourf(x,y,dcube_3,cmap=cm.coolwarm)
                    # ax3.contourf(x,y,img_gauss,cmap=cm.coolwarm)
                    # #
                    # surf_fk=ax2.plot_surface(x,y,dcube_3, rstride=1, cstride=1, cmap=cm.coolwarm,
                    #                          linewidth=0)
                    # surf_ns=ax4.plot_surface(x,y,img_gauss, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0)

                    # plt.show()
                    '''======================================================'''

                    #
                    # Calculate aperture and PSF flux
                    # Different radius
                    r1,r2,r3,r4=1.25,2.5,3.75,5.0
                    #
                    v_ap[0].append( Calc_Phot(dcube_3,x_cm,y_cm,radius=r1) )
                    v_ap[1].append( Calc_Phot(dcube_3,x_cm,y_cm,radius=r2) )
                    v_ap[2].append( Calc_Phot(dcube_3,x_cm,y_cm,radius=r3) )
                    v_ap[3].append( Calc_Phot(dcube_3,x_cm,y_cm,radius=r4) )
                    #
                    v_psf[0].append( Calc_Phot(img_gauss,x_cm,y_cm,radius=r1) )
                    v_psf[1].append( Calc_Phot(img_gauss,x_cm,y_cm,radius=r2) )
                    v_psf[2].append( Calc_Phot(img_gauss,x_cm,y_cm,radius=r3) )
                    v_psf[3].append( Calc_Phot(img_gauss,x_cm,y_cm,radius=r4) )
                    #
                else:
                    print '\n\tNo flux measurement in layer {0}'.format(dd)
                    continue
            # Save array
            print '\n\t___Saving array of photometry.'
            print '\tElapsed time: {0}min [t_proc={1}]'.format((time.time()-t_init)/60.,
                                                               time.clock()-clock_init)
            #
            out_nm=path_npy+'lc_'+filename[filename.find('c_')+2:]
            # each array will be a column
            out_arr=np.column_stack(( jd_arr[time_cal],np.array(v_ap[0]),np.array(v_ap[1]),np.array(v_ap[2]),np.array(v_ap[3]),np.array(v_psf[0]),np.array(v_psf[1]),np.array(v_psf[2]),np.array(v_psf[3]) ))
            print '\t___Saving LC: {0}'.format(out_nm)
            #
            np.save(out_nm,out_arr)
    #_______________________________________________________
    #_______________________________________________________
    #_______________________________________________________
    #_______________________________________________________

    print '\n\n_______________________\nEnd. {0}'.format(time.ctime())

    # For interpolation
    # http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.ndimage.interpolation.map_coordinates.html


    '''
    ADDITIONAL:
    Use gausKer_v01 to set a gaussian kernel of the profile, row by row and
    then match column by column
    Compare with kernel3D_v01
    Gauss_Profile(im_star,sigma_psf): sigma_psf is the critical value, because
    it undersample/oversample the gaussian profile model.
    
    center of mass is calculated through ndimage routine. If NaN, False is returned
    '''
