def WVT_2(data2d):
    '''CONTINOUS wavelet transform
    
    Use this transform to compare to discrete one. Perform significance areas
    Doubts: 
    - may I interpolate? It is not needed since I have continuos data
    - the scales to interpolate are in power of 2, eight? But how translate 
    to len()?
    - how avoid border effects? pad with zeroes?

    IMPORTANT
    ---------
    - as I'll perform 1-dimensional wvt, need to consider that cosmic rays 
    propagate on the CCD following the LONGER axis. On a first approach
    lets make the wvt on the axis where if there is a 'bad' column it will
    only affect the wvt containing this column. So wvt along the longer axis,
    just ast the readout.
    - plot the 1D data. Plot range without ourliers (but keep them in data)
    
    FOR THE ENTIRE FOCAL PLANE
    --------------------------
    Must fill the interCCD space with zeroes or interpolation.
    
    Scales must define the refinement scales I will look. The more scales, the 
    slower calculation and better the resolution.
    '''
    ROWS,COLS = data2d.shape #(4096,2048)
    window = np.arange(2,(COLS)+1)
    for small_ax in xrange(0,COLS):
        data1d = data2d[:,small_ax]
        #cwt_arr = scipy.signal.cwt(data1d,scipy.signal.morlet,window)
        plt.plot(data1d,'b.-')
        plt.show()
    scipy.signal.find_peaks_cwt() #to find peaks in the wavelet map
    data1d = data2d
    L = 100
    window = np.arange(4,L)
    scipy.signal.cwt(data1d,sipy.signal.morlet,window)
    return False
