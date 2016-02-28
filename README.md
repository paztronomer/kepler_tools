# kepler_tools
Tools for Kepler raw data reading, reduction and model. 
The scripts are mainly oriented to analysis of photometric variability in light curves being a tracers of stellar rotation.

This repo includes python scripts: 
  - fits2csv_v05.py: read FITS files and convert it to tables, to merge different quarters using the approach of Banyai et al. 2013, erasing outliers by Iglewicz & Hoaglin criteria
  - fastLombScargle.py: performs fast Lomb Scargle (Press & Ribicky 89')
  - lscargle_kplr.py: performs fast and Generalised (astroML) Loms Scargle on time series. Save info for 1st peak (include calculation time)
  - addCol.py: match tables using nearest neighbor.

I will upload (soon, 'cause the motherboard of my notebook fried days ago): useful scripts, among them one to estimate errors by Monte Carlo Markov Chains approach of Foreman-Mackey et al. 2013** (emcee project on github)

Scripts are simpler and commented, please refer to me if use it of your research work.

**visit his project "emcee", on github
