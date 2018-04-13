''' 
Create a noisy simulated star, calculate its PSF and centroid. Then perform
aperture and PSF photometry.

The MIT License (MIT)

Copyright (c) 2014-2015 Francisco Paz-Chinchon, UFRN, Brazil, francisco_at_dfte.ufrn.br

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

from __future__ import division
from scipy import stats #multivariate_normal
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import time



def NearNeighbors_Pix(arr,xx,yy,n=3):
    ''' Given a 2D-array, returns an nxn array whose "center" element is arr[x,y]
    Taken from http://stackoverflow.com/questions/4148292/how-do-i-select-a-window-from-a-numpy-array-with-periodic-boundary-conditions

    np.roll(arr, shift, axis) :roll array elements along a given axis. Shift=number of places by which elements are shifted
    coord_final+coord_inic=1, because the element of interest becames the element (1,1) of the rolled array
    '''
    arr=np.roll( np.roll(arr,shift=1-xx,axis=0), shift=1-yy,axis=1)
    return arr[:n,:n]



def create_Star():
    ''' define pixel centroids
    I'll use np.mgrid instead of np.meshgrid, because is faster to initialize. There are 2 ways: 
    [a:b:Nj]   N:length
    [a:b:D]    D:step size
    '''
    x, y=np.mgrid[0:13:1, 0:19:1]
    centers=np.dstack((x, y))
    # dstack() : deep stacking

    ''' create our fake star with noise. Our light profile is tilt using the anti-diagonal of covariance
    '''
    # define multivariate adjust
    mean_pdf=[6, 19/2]
    cov_pdf=[[10, 1.99], [1.99, 30]]
    gProf=stats.multivariate_normal(mean=mean_pdf, cov=cov_pdf) 
    # gProf.pdf([xi,yi]) calls the values of pdf in that point
    noisyProf = gProf.pdf(centers) + 0.5 * gProf.pdf(centers).std() * np.random.random(gProf.pdf(centers).shape)
    # add 2 outliers
    noisyProf[12,3]=10.*gProf.pdf(centers).std()
    noisyProf[5,8]=25.*gProf.pdf(centers).std()
    #
    return  gProf.pdf(centers),noisyProf,x,y


    
def Gauss_Profile(im_star,sigma_psf):
    from scipy import ndimage
    ''' Set by hand and carefully the sigma (width) of the gaussian 
    profile to avoid over/under samplig
    '''
    img_gauss = ndimage.filters.gaussian_filter(im_star,sigma=sigma_psf,order=0,mode='nearest')
    return img_gauss



def star_CM(img_soft):
    ''' Called with the gaussian profile
    '''
    # Center of mass
    from scipy.ndimage.measurements import center_of_mass
    cmass=center_of_mass(img_soft)
    #
    # Maximum of Measurements
    from scipy.ndimage import measurements
    max_coord=measurements.maximum_position(img_soft)
    #
    return cmass,max_coord



def Detect_Outliers(imlayer):
    ''' Estimate Iglewicz and Hoaglin criteria for outlier 
    http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
    Formula:
    Z=0.6745(x_i - median(x)) / MAD
    if abs(z) > 3.5, x_i is a potential outlier

    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and Handle Outliers", The ASQC Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    '''
    from statsmodels.robust import scale #mad
    from scipy.stats import norm

    # Percent point function (inverse of cdf -- percentiles) of a normal continous random variable
    cte=norm.ppf(0.75) #aprox 0.6745
    # Flatten the image, to estimate median
    flat_im=imlayer.ravel()
    
    # MAD
    MAD=np.median( np.abs( flat_im-np.median(flat_im) ) )
    # alternative: scale.mad(flat_im, c=1, axis=0, center=np.median)
    
    # the  Iglewicz and Hoaglin scorer. If absolute value is higher than 3.5 it is a possible outlier
    Zscore=cte*(flat_im-np.median(flat_im))/MAD
    Zscore=np.abs(Zscore)
    
    ''' search for outliers and if present replace by the median of neighbors
    '''
    if len(Zscore[Zscore>3.5])>0:
        for k in range(0,imlayer.shape[0]):
            for m in range(0,imlayer.shape[1]):
                if np.abs( cte*(imlayer[k,m]-np.median(flat_im))/MAD ) >3.5:
                    print '___values of possible outliers (noisy star profile, coords [{1},{2}]): {0:.4f} '.format\
                        (imlayer[k,m],k,m)
                    ''' if coordinate is on the edges then matrix is replicated
                    function NearNeighbors_Pix() creates an 3x3 array on which element (1,1) is our pixel of interest
                    '''
                    neighbor=NearNeighbors_Pix(imlayer,k,m)
                    # set NaN the outlier, located at [1,1] of sub matrix
                    neighbor[1,1]=np.nan
                    # then, estimate median, EXCLUDING the NaN value (because numpy mess up with NaNs)
                    # and set this value to the outlier
                    flat_neighbor=neighbor.ravel(); flat_neighbor=flat_neighbor[~np.isnan(flat_neighbor)]
                    imlayer[k,m]=np.median(flat_neighbor)
                    print '___replaced value: {0:.4f} '.format(imlayer[k,m])
        return imlayer
    else:
        return imlayer



def Calc_Phot(image,xc,yc,radius):
    ''' aperture or PSF measurement of flux!
    '''
    #integrate inside radius
    #http://stackoverflow.com/questions/8647024/how-to-apply-a-disc-shaped-mask-to-a-numpy-array
    #
    ''' define 2 one-dimensional arrays based in centroid of distribution (xc,yc)
     where each one have values:
     -xc...(Nx-xc)-1  and -yc...(Ny-yc)-1 Nx,Ny being number of elements on each axis
     np.ogrid returns meshgrid arrays of only 1 dimension
    '''
    rows,cols = np.ogrid[-yc:image.shape[0]-yc, -xc:image.shape[1]-xc]
    mask_r = cols**2 + rows**2 <= radius**2
    ''' then apply this mask composed of True/False to the 2d image and get a submatrix. 
    Then add over 2 axis
    '''
    return np.sum(image[mask_r])



def FT_Gaussian(im,sigma_gauss):
    # Fourier Transform of a Gaussian
    from scipy.ndimage.fourier import fourier_gaussian
    return fourier_gaussian(im, sigma=sigma_gauss)



def main():
    #
    # ===================
    # S I M U L A T E 
    # S T A R
    # ===================
    #
    ''' call function to create fake star with and w/o noise
    '''
    fk_star,ns_star,x,y=create_Star()

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
    img_corr=Detect_Outliers(imlayer=ns_star)
        
    '''gaussian convolution
    there is a parameter to control outliers influence??
    '''
    guess_sigma=0.75
    img_gauss=Gauss_Profile(ns_star,sigma_psf=guess_sigma)
    print '\n======\nImportant issue: sigma of the PSF will highly influence \nthe goodness of fit between PSF and stellar profile. \nOver/under-sample will depend of it. So use the plots \nof difference between profiles to make your choice.\n\tCurrent sigma = {0}  \n======'.format(guess_sigma)

    # Center of mass for noise-free profile and gaussian profile
    CM,coord_peak=star_CM(img_gauss)
    CM_nsfree,coord_peak_nsfree=star_CM(fk_star)
    print '\n\tCM of gaussian profile {0}, and the maximum of distribution {1}'.format(CM,coord_peak)
    print '\tCM of noise-free stellar profile {0}, and the maximum of distribution {1}\n'.format(CM_nsfree,coord_peak_nsfree)

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
    print psf_cal, ap_cal


    #
    # ===============
    # P L O T S
    # ===============
    #
    
    ''' test plot for:
    ideal stellar profile: contour map + 3D profile
    noisy stellar profile: contour map + 3D profile
    '''
    if False:
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

    ''' test plot for:
    noisy star contour plot
    gaussian filter to noisy star contour plot + center of mass of gaussian profile AND original star
    original noise-free star contour plot
    noisy star multiplied by Fourier transform of a gaussian
    '''
    if False:
        plt.close('all')
        fig=plt.figure()
        ax1=fig.add_subplot(221)
        ax2=fig.add_subplot(222)
        ax3=fig.add_subplot(223)
        ax4=fig.add_subplot(224)
        #
        ax1.imshow(img_corr); ax1.set_title('noisy star')
        ax2.imshow(img_gauss); ax2.set_title('gaussian filter noisy star')
        ax2.plot(CM[1],CM[0],'w+',ms=20,alpha=1)
        ax2.plot(CM_nsfree[1],CM_nsfree[0],'wx',ms=20,alpha=1)
        ax3.imshow(fk_star); ax3.set_title('original star')
        ax4.imshow(FT_star); ax4.set_title('multiplied by FT of gaussian')
        #
        ax2.set_xlim([0,img_gauss.shape[1]-1]); ax2.set_ylim([img_gauss.shape[0]-1,0])
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
    if False:
        fig=plt.figure(figsize=(12,5))
        ax1=fig.add_subplot(141,projection='3d')
        ax2=fig.add_subplot(142,projection='3d')
        ax3=fig.add_subplot(143,projection='3d')
        ax4=fig.add_subplot(144,projection='3d')
        # noise free star
        ax1.plot_surface(x,y,fk_star, rstride=1, cstride=1, cmap=cm.coolwarm,
                         linewidth=0); ax1.set_title('Noise free star')
        # noisy star
        ax2.plot_surface(x,y,img_corr, rstride=1, cstride=1, cmap=cm.coolwarm,
                         linewidth=0); ax2.set_title('Noisy star')
        # gaussian PSF + gaussian PSF minus fake (noise free) star
        ax3.plot_surface(x,y,img_gauss-fk_star, rstride=1, cstride=1, cmap=cm.gray,
                         linewidth=0)
        ax3.plot_surface(x,y,img_gauss, rstride=1, cstride=1, cmap=cm.coolwarm,
                         linewidth=0,alpha=0.8); ax3.set_title('Gaussian PSF +\nPSF$-$noise-free star')
        # gaussian PSF + gaussian PSF minus noisy star
        ax4.plot_surface(x,y,img_gauss-img_corr, rstride=1, cstride=1, cmap=cm.gray,
                         linewidth=0)
        ax4.plot_surface(x,y,img_gauss, rstride=1, cstride=1, cmap=cm.coolwarm,
                         linewidth=0,alpha=0.8); ax4.set_title('Gaussian PSF +\nPSF$-$noisy star')
        ax1.view_init(elev=10., azim=-65.)
        ax2.view_init(elev=10., azim=-65.)
        ax3.view_init(elev=10., azim=-65.)
        ax4.view_init(elev=10., azim=-65.)
        #fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()
    


if __name__=='__main__':
    main()






'''
# does KDE helps in our proposal???
# estimate PDF of the image values

from statsmodels.nonparametric.kernel_density import KDEMultivariate    
kde = KDEMultivariate(data=im,var_type='ccccccccccccc',bw='normal_reference')
#print kde.bw, kde.pdf(ns_star)

exit(0)
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

# bivariate gaussian of dimensions=number of pixels (mxn) where each coordinate 
of the matrix have a value given by a gaussian profile

# __TRY #1
mean = [0,0]
cov = [[1,0],[0,1]]
gx = np.random.multivariate_normal(mean,cov,len(x1)).T
gy = np.random.multivariate_normal(mean,cov,len(y1)).T

gxy = np.random.multivariate_normal(mean,cov,x2.shape)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')#fig.gca(projection='3d')
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)

mean =[0,0] #np.array([[0],[0]])# [0,0]
cov = [[1,0],[0,10]]
Z =construct_Z(X,Y,mu=mean,cov=cov)

print X.shape, Y.shape,Z.shape
surf = ax.plot_surface(X,Y,Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

#add noise with sigma=1/3 per aech dimension

noise_plus_signal = signal_in + numpy.random.normal(loc=0, scale=sigma/numpy.sqrt(3.0), size=(n,3))

(m,n) dimensions create (m,n,N) where N is the sample size

mean = [0,0]
cov = [[1,0],[0,100]]
x = np.random.multivariate_normal(mean,cov,(m,n)).T
plt.plot(x,y,'x'); plt.axis('equal'); plt.show()

#to setup a grid

X, Y = np.meshgrid(X, Y)

ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
'''
