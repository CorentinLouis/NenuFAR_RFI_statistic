# NenuFAR_RFI_statistic
Routines for analyzing the radio frequency interference (RFI) environement at the NenuFAR radiotelescope using LT02 (Exoplanets) observations


`plot_rfi_stat.py` is the main routine. This script performs calculations and plotting of RFI statistical distributions from NenuFAR LT 02 observations.

It includes functions for calculating (`rfi_stat_distribution_functions.py`) and saving (`hdf_functions.py`) RFI statistical distributions and for plotting (`rfi_stat_plot_distributions_functions.py`) these distributions.

The script takes various arguments for controlling the calculations and plotting of the distributions.

Usage:  
    `python plot_rfistat.py [options]`


Options: 
Option | Description
--- | ---  
`-h, --help` | Show this help message and exit.
`--input_sav_file_directory`  | Input sav file directory location to calculate the distributions (default: `"./"`). Useful if calculate_distributions keyword is True
`--input_hdf5_file_directory` | Input hdf5 file directory location to plot the distributions (default: `"./"`). Useful if calculate_distributions keyword is True
`-rfi_level` | Level of RFI mitigation (required).
`--ouput_directory`           | Output directory (default: `"./"`).
`--color_map`                 | Color map for plots (default: `'viridis'`).
`--calculate_distributions`   | Keyword to calculate and save the distributions (default: `False`). If False, and plot_distributions keyword is True, the plotting function will try to load pre-calculated distribution from hdf5 files.
`--calculate_distributions`   | Keyword to plot the distributions (default: `False`).
`--local_time`                | Keyword to plot only on local time (default: `False`).
`--azimuth`                   | Keyword to plot only on azimuth (default: `False`).
`--elevation`                 | Keyword to plot only on elevation (default: `False`).
`--elaz`                      | Keyword to plot only on elevation vs. azimuth plane (default: `False`).

If `--local_time`, `--azimuth`, `--elavation` and `--ezal` `== False`, the plotting functions will plot all these distributions  

Example:   
```    python plot_rfistat.py -rfi_level 0 --calculate_distributions --plot_distributions --local_time --azimuth ```

Requirements:
* python 3.8.10
* pip install -r requirements_packages_exoradio.txt
