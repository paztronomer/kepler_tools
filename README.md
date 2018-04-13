# Kepler data tools for fast analysis (old codes from my grad period)
## Tools for Kepler raw data reading, reduction, model, and plotting .The repository also includes codes for other Time Series purposes

The scripts are mainly oriented to analysis of photometric variability in light
curves as tracers of stellar rotation, and for support analysis and display

1. **fits2csv_v05.py**: read FITS files and convert it to tables, to merge
    different quarters using the approach of *Banyai et al. 2013*, also
    removing outliers by the use of *Iglewicz & Hoaglin* criteria.

1. **fastLombScargle.py**: code to perform fast Lomb Scargle (*Press & Ribicky
    1989*). This code is intended to be imported into others.

1. **lscargle_kplr.py**: performs Fast and Generalized (**astroML** routine)
    Lomb Scargle on time series. Saves information for first periodogram peak.

1. **addCol.py**: match tables using nearest neighbor.

1. **zoom_C_v01.py**, **mcmc_LC_v01.py**: simple code to plot a LC with a
    zoom-in of a section, on another subplot.

1. **writeWgetFromList_corrected.py**: code to write the bash wget files to be
    used for get Kepler raw data.

1. **morlet_CWT.py**: simple function to perform Continuous Wavelet Transform,
    using Morlet mother wavelet, on 2D images.

1. **photoKplr_targPix_v03.py**: code to perform aperture photometry on
    Target-Pixel_files from Kepler.

1. **quickAnalysisPhot_v01.py**: code for quickly assessment of Target-Pixel
    Kepler files, quarter by quarter.

1. **kplr_SNsel_v09.py**, **FITStoCSV_AndLS.py**: code for Time Series period
    analysis. Performs period fitting and calculation, returning periodogram
    related information, false-positive estimator, and other useful values.
    Works directly with FITS files.

1. **weibull_v01.py**: simple code to plot Weibull distribution

1. **Mgiants_simbadQuery_v03.py**: code to connect to SIMBAD and retrieve info
    for M giants. Modifiable to other spectral types.

1. **UncertSine_mcmc_v01.py**: early version of the MCMC error estimation,
    based on Fourier decomposition of light curves. Uses Monte Carlo Markov
    Chains approach of *Foreman-Mackey et al. 2013* (emcee project on Github)

1. **kernel3D_v01.py**: create a noisy simulated star, calculate its PSF and
    centroid. Then perform aperture and PSF photometry.

1. **create_colormap.py**: code to customize a colormap for matplotlib.

1. **send_email.py**: simple code to send emails from Python. Can be used to
    alert when a code has finished

1. **HexbAndCMD_v01.py**: code to plot hexbinned color-magnitude diagrams and
    get the overdense tracer, for fitting tracks. Uses R for hexbins.


**Notes**
* these codes does not follows PEP8 in some aspects,
* these codes were written without much optimization in mind
* keep in mind were written with no sharing in mind
