''' 
Quickly LC view coming from TargetPixel photometry (or raw photometry),
quarter by quarter.
I recommend to use not-detrended LCs and thus apply a 3-degree polynomial
for visualization.
'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from progressbar import ProgressBar
import time

def FitPol(xx,yy,degree):
    '''np.poly1d(z)(N) = eval in N
    '''
    z=np.polyfit(xx,yy,degree)
    return np.poly1d(z)

def DetectOutliers(flux):    
    ''' Detect possible outliers through the Iglewicz & Hoaglin criteria
    and return the mask to be applied.
    Input=1D array
    '''
    from statsmodels.robust import scale #mad
    from scipy.stats import norm
    # Percent point function (inverse of cdf -- percentiles) of a normal continous random variable
    cte=norm.ppf(0.75) #aprox 0.6745
    
    MAD=np.median( np.abs( flux-np.median(flux) ) )
    # alternative: scale.mad(flat_im, c=1, axis=0, center=np.median)
    
    # the  Iglewicz and Hoaglin scorer. If absolute value is higher than 3.5 it is a possible outlier
    Zscore=cte*(flux-np.median(flux))/MAD
    Zscore=np.abs(Zscore)
    #
    mask=Zscore<3.5
    return mask
    


if __name__=='__main__':
    print '\n\tStarting... {0}'.format(time.ctime())
    path_lc='path_to_LC/'
    
    for filename in glob.glob(os.path.join(path_lc, 'lc_*')):
        '''Load numpy files
        '''
        lc=np.load(filename)
        
        plt.close('all')
    
        # Outliers remotion
        mask=DetectOutliers(lc[:,6])
    
        # for degree 3 polynomial normalization
        ff=lc[mask,8]-lc[mask,6]
        deg3=FitPol(lc[mask,0],ff,3)
        f_corr=ff/deg3(lc[mask,0])

        fig=plt.figure()
        ax1=fig.add_subplot(211)
        ax2=fig.add_subplot(212)

        #ax1.plot(lc[:,0],lc[:,8],'r.')
        ax1.plot(lc[mask,0],ff,'b.')
        ax2.plot(lc[mask,0],f_corr,'g.')

        plt.show()
