# Script to plot LC and a zoom to it
# source code from: http://matplotlib.org/examples/pylab_examples/axes_zoom_effect.html

from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector, BboxConnectorPatch


def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = prop_lines.copy()
        prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    p = BboxConnectorPatch(bbox1, bbox2,
                           #loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           **prop_patches)
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect01(ax1, ax2, xmin, xmax, **kwargs):
    """
    ax1 : the main axes
    ax2 : the zoomed axes
    (xmin,xmax) : the limits of the colored area in both plot axes.

    connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
    be marked.  The keywords parameters will be used ti create
    patches.

    """

    trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
    trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

    bbox = Bbox.from_extents(xmin, 0, xmax, 1)

    mybbox1 = TransformedBbox(bbox, trans1)
    mybbox2 = TransformedBbox(bbox, trans2)

    prop_patches=kwargs.copy()
    prop_patches["ec"]="none"
    prop_patches["alpha"]=0.15

    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                     prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect02(ax1, ax2, **kwargs):
    """
    ax1 : the main axes
    ax2 : the zoomed axes

    Similar to zoom_effect01.  The xmin & xmax will be taken from the
    ax1.viewLim.
    """

    tt = ax1.transScale + (ax1.transLimits + ax2.transAxes)
    trans = blended_transform_factory(ax2.transData, tt)

    mybbox1 = ax1.bbox
    mybbox2 = TransformedBbox(ax1.viewLim, trans)

    prop_patches=kwargs.copy()
    prop_patches["ec"]="none"
    prop_patches["alpha"]=0.2

    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                     prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p
'''
 Up to here the code was taken from webpage at top (and slightly modified)
======================================================================
'''

#
# FROM HERE I SET MY PARAMETERS
#
import numpy as np
import matplotlib.pyplot as plt
import time
import glob
import os

folder='npy_files/'
path_kcp='sel/'
print '\n{0}\n'.format(time.ctime())

if __name__=='__main__':
    # Plot
    # fig=plt.figure()
    # ax1=fig.add_subplot(111)
    # nrows,ncols=2,1
    # fig, axes = plt.subplots(nrows,ncols)
    # for axrow in axes:
    #     axrow[0].plot(x, color='red')
    #     axrow[1].plot(y, color='green')
    ''' axes[LINE][COL].plot()
    '''

    fig=plt.figure(figsize=(16,8))
    ax1=fig.add_subplot(212)
    ax2=fig.add_subplot(222)
    ax3=fig.add_subplot(241)
    ax4=fig.add_subplot(242)
    # ax3=fig.add_subplot(221)
    
    for filename in glob.glob(os.path.join(folder,'*5792202*')):
        t0=time.clock()
        sam=np.load(filename)
        print 'LC KIC: {0}'.format(int(filename[filename.find('kplr')+4:filename.find('_N1000')]))
        # samples:: Ampl, Freq, lnf

        for fn in glob.glob(os.path.join( path_kcp, '*'+filename[filename.find('kplr'):filename.find('_N1000')]+'*' )):
            lc=np.loadtxt(fn)
            #plt.plot(lc[:,0],lc[:,1],'b.'); plt.show()

        if int(filename[filename.find('kplr')+4:filename.find('_N1000')])==5792202: ini,prot=829.502,14.77
        elif int(filename[filename.find('kplr')+4:filename.find('_N1000')])==6949607: ini,prot=598.538,17.97
        else: prot=None

        t1,t2=ini,ini+3.*prot
    
        for dd in xrange(0,sam.shape[0]/25,1):
            if dd%1000==0: print '\t__processor time N {1}: {0}'.format(time.clock()-t0,dd)
            time_f=np.linspace(t1,t2,100)
            flux_f=sam[dd,0]*np.sin(2.*np.pi*sam[dd,1]*time_f+np.pi)+1.
            # lc on upper panel
            ax1.plot(time_f,flux_f,'-',color='0.5',alpha=0.1,lw=1.5)
            ax1.plot(lc[:,0],lc[:,1],'r.')
            ax2.plot(time_f,flux_f,'-',color='0.5',alpha=0.1,lw=1.5)
            ax2.plot(lc[:,0],lc[:,1],'r.')
        
        # Zoom for the lower panel
        # Call ZOOM
        zoom_effect01(ax2, ax1, 852., 860.)
        # zoom_effect02(ax3, ax1)

    print 'end read tabs... {0}'.format(time.ctime())


    ''' Boxplot    
    '''
    ax3.boxplot(sam[:,0],showmeans=False,notch=True,showfliers=False,labels=['MCMC Amplitude'])#,'g.') # amplitude
    ax4.boxplot(1./sam[:,1],showmeans=False,notch=True,showfliers=False,labels=['MCMC P$_{rot}$'])#1,'')#,'g.') # period

    
    print '\nstarting setup {0}'.format(time.ctime())

    ''' quick plot setup
    '''
    plt.setp( ax1.get_xticklabels(), fontsize=14)
    plt.setp( ax1.get_yticklabels(), fontsize=14)
    plt.setp( ax2.get_xticklabels(), fontsize=14)
    plt.setp( ax2.get_yticklabels(), fontsize=14)#visible=False)
    plt.setp( ax3.get_yticklabels(), fontsize=14)
    plt.setp( ax3.get_xticklabels(), fontsize=14)
    plt.setp( ax4.get_yticklabels(), fontsize=14)
    plt.setp( ax4.get_xticklabels(), fontsize=14)

    ax1.set_xlim([ini,3.*prot+ini])
    ax2.set_xlim([852., 860.])
    ax2.set_ylim([0.96,1.03])

    delta_jd=2454833.0
    ax1.set_xlabel('JD-'+str(delta_jd)+' (d)',fontsize=15)
    ax1.set_ylabel('Flux$_{norm}$ (e$^{-}$ cadence$^{-1}$)',fontsize=15)
    ax2.set_xlabel('JD-'+str(delta_jd)+' (d)',fontsize=15)
    ax2.set_ylabel('Flux$_{norm}$ (e$^{-}$ cadence$^{-1}$)',fontsize=15)
    ax3.set_ylabel('Amplitude (e$^{-}$ cadence$^{-1}$)',fontsize=15)
    ax4.set_ylabel('P$_{rot}$ (d)',fontsize=15)
    

    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    #
    # to format : '%d' or '%1.2f' or '%1.1f cm'
    majorLocator_x1   = MultipleLocator(5)
    majorFormatter_x1 = FormatStrFormatter('%d')
    minorLocator_x1   = MultipleLocator(1)
    #
    majorLocator_x2   = MultipleLocator(1)
    majorFormatter_x2 = FormatStrFormatter('%d')
    minorLocator_x2   = MultipleLocator(.5)
    #
    majorLocator_y1   = MultipleLocator(0.01)
    majorFormatter_y1 = FormatStrFormatter('%.2f')
    minorLocator_y1   = MultipleLocator(0.0025)
    #    
    majorLocator_y2   = MultipleLocator(0.01)
    majorFormatter_y2 = FormatStrFormatter('%.2f')
    minorLocator_y2   = MultipleLocator(0.0025)
    #    
    majorLocator_y3   = MultipleLocator(0.002)
    majorFormatter_y3 = FormatStrFormatter('%.3f')
    minorLocator_y3   = MultipleLocator(0.0005)
    #    
    majorLocator_y4   = MultipleLocator(0.0005)
    majorFormatter_y4 = FormatStrFormatter('%.4f')
    minorLocator_y4   = MultipleLocator(0.00024)
    #
    ax1.xaxis.set_major_locator(majorLocator_x1)
    ax1.xaxis.set_major_formatter(majorFormatter_x1)
    ax2.xaxis.set_major_locator(majorLocator_x2)
    ax2.xaxis.set_major_formatter(majorFormatter_x2)
    #
    ax1.yaxis.set_major_locator(majorLocator_y1)
    ax1.yaxis.set_major_formatter(majorFormatter_y1)
    ax2.yaxis.set_major_locator(majorLocator_y2)
    ax2.yaxis.set_major_formatter(majorFormatter_y2)
    ax3.yaxis.set_major_locator(majorLocator_y3)
    ax3.yaxis.set_major_formatter(majorFormatter_y3)
    ax4.yaxis.set_major_locator(majorLocator_y4)
    ax4.yaxis.set_major_formatter(majorFormatter_y4)
    #
    #for the minor ticks, use no labels; default NullFormatter
    ax1.xaxis.set_minor_locator(minorLocator_x1)
    ax2.xaxis.set_minor_locator(minorLocator_x2)
    #
    ax1.yaxis.set_minor_locator(minorLocator_y1)
    ax2.yaxis.set_minor_locator(minorLocator_y2)
    ax3.yaxis.set_minor_locator(minorLocator_y3)
    ax4.yaxis.set_minor_locator(minorLocator_y4)
    ##


    # adjust
    plt.subplots_adjust(left=0.07, bottom=0.08, right=0.99, top=0.98, wspace=0.45, hspace=0.24)


    #plt.show()
    print 'saving... {0}'.format(time.ctime())

    if True:
        fig.savefig('mcmc_LC_v02.jpg', dpi=100, facecolor='w', edgecolor='w', \
                    orientation='portrait', papertype='letter', format='jpg',\
                    transparent=False, bbox_inches=None, pad_inches=0.1, \
                    frameon=None)

