# -*- coding: utf-8 -*-
"""
This script performs calculations and plotting of RFI statistical distributions for NenuFAR observations.

It includes functions for calculating and saving RFI statistical distributions and for plotting these distributions.
The script takes various arguments for controlling the calculations and plotting.

Usage:
    python plot_rfistat.py [options]

Options:
    -h, --help                  Show this help message and exit.
    --input_sav_file_directory  Input sav file directory location to calculate the distributions (default: "./"). Useful if calculate_distributions keyword is True
    --input_hdf5_file_directory Input hdf5 file directory location to plot the distributions (default: "./"). Useful if calculate_distributions keyword is True
    -rfi_level                  Level of RFI mitigation (required).
    --ouput_directory           Output directory (default: "./").
    --color_map                 Color map for plots (default: 'viridis').
    --calculate_distributions   Keyword to calculate and save the distributions (default: False). If False, and plot_distributions keyword is True, the plotting function will try to load pre-calculated distribution from hdf5 files.
    --calculate_distributions   Keyword to plot the distributions (default: False).
    --local_time                Keyword to plot only on local time (default: False).
    --azimuth                   Keyword to plot only on azimuth (default: False).
    --elevation                 Keyword to plot only on elevation (default: False).
    --elaz                      Keyword to plot only on elevation vs. azimuth plane (default: False).
If --local_time, --azimuth, --elavation and --ezal == False, the plotting functions will plot all these distributions

Example:
    python plot_rfistat.py -rfi_level 0 --calculate_distributions --plot_distributions --local_time --azimuth
"""


import numpy
from load_functions import *
from rfi_stat_distribution_functions import distribution_flag_ana_beam_all_beam, distribution_flag_nume_beam_all_beam, distribution_flag_elaz_nume_beam_all_beam, distribution_observations_numebeam, distribution_observations_anabeam
from rfi_stat_plot_distributions_functions import plot_distribution_freq, plot_distribution_elaz, plot_1D_time_distribution_observation, plot_1D_distribution_observation, plot_distribution_with_histograms, plot_distribution_with_mean_flag_histo, plot_distribution_mollweide, insert_gaps, plot_flag_freq_versus_time_distribution_with_histograms
from make_gif import make_gif
from doy_to_ymd import doy_specific_year_to_yyyyddd, doy_to_datetime, from_utc_to_local_time, datetime_to_float_since1970
from hdf_functions import save_to_hdf, read_hdf5_file
import argparse

def rebin_distribution_in_frequency(freq: numpy.array,
                                    distribution: numpy.array,
                                    dfreq: numpy.array):
    # Your original frequency array
    #original_freq = numpy.linspace(7.73773193359375, 88.93890380859375, 1664)

    # Create a new frequency array
    new_ffreq = numpy.arange(numpy.min(freq), numpy.max(freq)+dfreq, dfreq)

    # Initialize the rebinned distribution array
    rebinned_distribution = numpy.zeros((distribution.shape[0], distribution.shape[1], len(new_ffreq)-1)) + numpy.nan
    #flag_ = numpy.nanmean(flag, axis = 0)
    flag_ = distribution[:,:,:]
    flag_[flag_ <0] = numpy.nan

    # Loop over the new frequency bins and rebin the data
    for index, i_freq in enumerate(new_ffreq[0:-1]):
        rebinned_distribution[:, :, index] = numpy.nanmean(flag_[:, :, numpy.where((freq >= i_freq) & (freq <new_ffreq[index+1]))], axis=-1).reshape(4, distribution.shape[1])
    new_ffreq=new_ffreq[0:-1]
    # The rebinned_distribution array now contains the rebinned data.

    rebinned_distribution[rebinned_distribution <0] = numpy.nan

    return(new_ffreq, rebinned_distribution)

def read_rfistat_file(intput_directory: str = './',
                   filename: str = "postproc_02_step3_rfistat_rfilevel0.sav"):
    # Create a LazyLoader instance with the data file path and the loader function
    print("###### Restoring temporal RFI mask sav file ######")
    loader = LazyLoader(intput_directory+filename, load_rfistat_data)

    # At this point, no data has been loaded

    # When you access the 'data' attribute, it will trigger lazy loading
    data_dict = loader.data
    # The data is now loaded and available in the 'data_dict' variable
    # print(data_dict.keys())

    # Access specific data fields from the loaded data
    ffreq = data_dict['ffreq']
    time = data_dict['ti']
    beams = data_dict['beams']
    nibeg = data_dict['nibeg']
    nimax = data_dict['nimax']
    nti = data_dict['nti']
    flag = numpy.zeros((4, nti, ffreq.size))
    for ibeam in range(4):
        flag[ibeam,:,:] =  data_dict['flag'][ibeam,:,:].T
    elA = data_dict['elA']
    azA = data_dict['azA']
    elB = data_dict['elB']
    azB = data_dict['azB']

    return(ffreq, time, beams, nibeg, nimax, nti, flag, elA, azA, elB, azB)

def calculate_rfi_stat_distribution(rfilevel,
                                    intput_directory,
                                    output_directory,
                                    save_hdf5
                                    ):
    """
    Calculate RFI statistical distributions from RFI temporal sav file
      and save them into a HDF5 file.

    Parameters:
        rfilevel (int): Level of RFI processing, which determines the flag calculation method:
            - 0: Flag percentage of bad pixels calculated during the processing step (step 1) while integrating over 96 pixels.
            - 1: Flag calculated during the post-processing step (step 2) using PATROL.
            - 2: Flag including level 1 and a threshold at 7.5%.
            - 3: Flag including level 1 and a threshold at 6%.
        intput_directory (str): The directory containing the input data files.
        output_directory (str): The directory where the output data will be saved.

    Returns:
        Tuple: A tuple containing the following elements:
            - rfilevel (int): Level of RFI processing.
            - ffreq (numpy.ndarray): Frequency array of the data (nf = 1664).
            - distribution_flag_data_localtime_all_ana_beam (numpy.ndarray): Distribution of RFI flags vs. local time for analogic beams.
            - distribution_flag_data_azimuth_all_ana_beam (numpy.ndarray): Distribution of RFI flags vs. azimuth for analogic beams.
            - distribution_flag_data_elevation_all_ana_beam (numpy.ndarray): Distribution of RFI flags vs. elevation for analogic beams.
            - distribution_flag_data_elevation_all_nume_beam (numpy.ndarray): Distribution of RFI flags vs. elevation for numeric beams.
            - distribution_flag_data_azimuth_all_nume_beam (numpy.ndarray): Distribution of RFI flags vs. azimuth for numeric beams.
            - distribution_flag_data_elaz_all_nume_beam (numpy.ndarray): Distribution of RFI flags vs. elevation vs. azimuth for numeric beams.

    """    
    (ffreq, time, beams, nibeg, nimax, nti, flag, elA, azA, elB, azB) = read_rfistat_file(intput_directory = intput_directory, filename = "postproc_02_step3_rfistat_rfilevel"+str(rfilevel)+".sav")
    
    
    # time is in days since 2019, January 1st
    time_2019 = time
    # switching to day of Year
    time_doy = doy_specific_year_to_yyyyddd(time_2019, 2019)
    # switching to datetime object
    time_datetime = numpy.array(doy_to_datetime(time_doy))
    # 'time' is now in UT
    # Calculating Local Time (+1 hours in Winter, +2 in Summer)
    local_time_datetime = from_utc_to_local_time(time_datetime)
    local_time = datetime_to_float_since1970(local_time_datetime)

    print("###### Calculating azimuth distribution of 6-minute slices observations for analogic beams######")
    (distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_azimuth) = distribution_observations_numebeam(flag,azB, type="azimuth")
    print("###### Calculating elevation distribution of 6-minute slices observations for analogic beams######")
    (distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_elevation) = distribution_observations_numebeam(flag,elB, type="elevation")
    print("###### Calculating frequency distribution of 6-minute slices observations for analogic beams######")
    (distribution_observation_numebeam_frequency, bins_distribution_observation_numebeam_frequency) = distribution_observations_numebeam(flag,ffreq, type="frequency")
    print("###### Calculating universal time distribution of 6-minute slices observations for analogic beams######")
    (distribution_observation_anabeam_UT, bins_distribution_observation_anabeam_UT) = distribution_observations_anabeam(flag,(time[:]-numpy.fix(time[:]))*24, type='universal time')
    print("###### Calculating local time distribution of 6-minute slices observations for analogic beams######")
    (distribution_observation_anabeam_LT, bins_distribution_observation_anabeam_LT) = distribution_observations_anabeam(flag,(local_time[:]-numpy.fix(local_time[:]))*24, type='local time')

    print("###### Calculating universal time vs. frequency distributions for analogic beams ######")
    (distribution_flag_data_universaltime_all_ana_beam, distribution_observations_data_universaltime_all_ana_beam) = distribution_flag_ana_beam_all_beam(flag, ffreq, (time[:]-numpy.fix(time[:]))*24, type='UT')
    print("###### Calculating local time vs. frequency distributions for analogic beams ######")
    (distribution_flag_data_localtime_all_ana_beam, distribution_observations_data_localtime_all_ana_beam) = distribution_flag_ana_beam_all_beam(flag, ffreq, (local_time[:]-numpy.fix(local_time[:]))*24, type='local time')
    print("###### Calculating azimuth vs. frequency distributions for analogic beams ######")
    (distribution_flag_data_azimuth_all_ana_beam, distribution_observations_data_azimuth_all_ana_beam) = distribution_flag_ana_beam_all_beam(flag, ffreq, azA, type='azimuth')
    print("###### Calculating elevation vs. frequency distributions for analogic beams ######")
    (distribution_flag_data_elevation_all_ana_beam, distribution_observations_data_elevation_all_ana_beam) = distribution_flag_ana_beam_all_beam(flag, ffreq, elA, type='elevation')

    print("###### Calculating elevation vs. frequency distributions for numeric beam ######")
    (distribution_flag_data_elevation_all_nume_beam, distribution_observations_data_elevation_all_nume_beam) = distribution_flag_nume_beam_all_beam(flag, ffreq, elB, type='elevation')
    print("###### Calculating azimuth vs. frequency distributions for numeric beam ######")
    (distribution_flag_data_azimuth_all_nume_beam, distribution_observations_data_azimuth_all_nume_beam) = distribution_flag_nume_beam_all_beam(flag, ffreq, azB, type='azimuth')

    print("###### Calculating azimuth vs. elevation distributions for numeric beam ######")
    (distribution_flag_data_elaz_all_nume_beam, distribution_observations_data_elaz_all_nume_beam) = distribution_flag_elaz_nume_beam_all_beam(flag, elB, azB, ffreq.size)

    
    if save_hdf5:
        print("###### Saving distributions into a HDF5 file ######")
        save_to_hdf(output_directory, rfilevel, time_datetime, ffreq, flag,
                    distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_azimuth,
                    distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_elevation,
                    distribution_observation_numebeam_frequency, bins_distribution_observation_numebeam_frequency,
                    distribution_observation_anabeam_UT, bins_distribution_observation_anabeam_UT,
                    distribution_observation_anabeam_LT, bins_distribution_observation_anabeam_LT,
                    distribution_flag_data_universaltime_all_ana_beam, distribution_observations_data_universaltime_all_ana_beam,
                    distribution_flag_data_localtime_all_ana_beam,distribution_observations_data_localtime_all_ana_beam,
                    distribution_flag_data_azimuth_all_ana_beam,distribution_observations_data_azimuth_all_ana_beam,
                    distribution_flag_data_elevation_all_ana_beam,distribution_observations_data_elevation_all_ana_beam,
                    distribution_flag_data_elevation_all_nume_beam,distribution_observations_data_elevation_all_nume_beam,
                    distribution_flag_data_azimuth_all_nume_beam,distribution_observations_data_azimuth_all_nume_beam,
                    distribution_flag_data_elaz_all_nume_beam,distribution_observations_data_elaz_all_nume_beam
                    )

    return(rfilevel,
           time_datetime, ffreq, flag,
           distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_azimuth,
           distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_elevation,
           distribution_observation_numebeam_frequency, bins_distribution_observation_numebeam_frequency,
           distribution_observation_anabeam_UT, bins_distribution_observation_anabeam_UT,
           distribution_observation_anabeam_LT, bins_distribution_observation_anabeam_LT,
           distribution_flag_data_universaltime_all_ana_beam, distribution_observations_data_universaltime_all_ana_beam,
           distribution_flag_data_localtime_all_ana_beam,distribution_observations_data_localtime_all_ana_beam,
           distribution_flag_data_azimuth_all_ana_beam,distribution_observations_data_azimuth_all_ana_beam,
           distribution_flag_data_elevation_all_ana_beam,distribution_observations_data_elevation_all_ana_beam,
           distribution_flag_data_elevation_all_nume_beam,distribution_observations_data_elevation_all_nume_beam,
           distribution_flag_data_azimuth_all_nume_beam,distribution_observations_data_azimuth_all_nume_beam,
           distribution_flag_data_elaz_all_nume_beam,distribution_observations_data_elaz_all_nume_beam
           )

def plot_rfi_stat(rfilevel,
                   output_directory,
                   universal_time,
                   local_time,
                   azimuth,
                   elevation,
                   elaz_square, 
                   elaz_polar_projection,
                   plot_distributions_rfi,
                   plot_distributions_observations,
                   plot_histograms,
                   time_datetime,
                   ffreq,
                   flag,
                   distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_azimuth,
                   distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_elevation,
                   distribution_observation_numebeam_frequency, bins_distribution_observation_numebeam_frequency,
                   distribution_observation_anabeam_UT, bins_distribution_observation_anabeam_UT,
                   distribution_observation_anabeam_LT, bins_distribution_observation_anabeam_LT,
                   distribution_flag_data_universaltime_all_ana_beam, distribution_observations_data_universaltime_all_ana_beam,
                   distribution_flag_data_localtime_all_ana_beam,distribution_observations_data_localtime_all_ana_beam,
                   distribution_flag_data_azimuth_all_ana_beam,distribution_observations_data_azimuth_all_ana_beam,
                   distribution_flag_data_elevation_all_ana_beam,distribution_observations_data_elevation_all_ana_beam,
                   distribution_flag_data_elevation_all_nume_beam,distribution_observations_data_elevation_all_nume_beam,
                   distribution_flag_data_azimuth_all_nume_beam,distribution_observations_data_azimuth_all_nume_beam,
                   distribution_flag_data_elaz_all_nume_beam,distribution_observations_data_elaz_all_nume_beam,
                   cmap):
    """
    Plot RFI statistical distributions.

    Parameters:
        rfilevel (int): Level of RFI processing:
            - 0: Flag percentage of bad pixels calculated during the processing step (step 1) while integrating over 96 pixels.
            - 1: Flag calculated during the post-processing step (step 2) using PATROL.
            - 2: Flag including level 1 and a threshold at 7.5%.
            - 3: Flag including level 1 and a threshold at 6%.
        output_directory (str): The directory where the output plots will be saved.
        local_time (bool): Whether to plot distributions for local time.
        azimuth (bool): Whether to plot distributions for azimuth.
        elevation (bool): Whether to plot distributions for elevation.
        elaz (bool): Whether to plot distributions for elevation vs. azimuth.
            If local_time == azimuth == elavation == ezal == False, the plotting functions will plot all these distributions
        ffreq (numpy.ndarray): Frequency array of the data (nf = 1664).
        distribution_flag_data_localtime_all_ana_beam (numpy.ndarray): Distribution of RFI flags vs. local time for analogic beams.
        distribution_flag_data_azimuth_all_ana_beam (numpy.ndarray): Distribution of RFI flags vs. azimuth for analogic beams.
        distribution_flag_data_elevation_all_ana_beam (numpy.ndarray): Distribution of RFI flags vs. elevation for analogic beams.
        distribution_flag_data_elevation_all_nume_beam (numpy.ndarray): Distribution of RFI flags vs. elevation for numeric beams.
        distribution_flag_data_azimuth_all_nume_beam (numpy.ndarray): Distribution of RFI flags vs. azimuth for numeric beams.
        distribution_flag_data_elaz_all_nume_beam (numpy.ndarray): Distribution of RFI flags vs. elevation vs. azimuth for numeric beams.

    """

    if plot_distributions_rfi or plot_distributions_observations:
        shading_distrib_freq = "auto"
        if universal_time == local_time == azimuth == elevation:
            universal_time = local_time = azimuth = elevation = True

    name_pdf = output_directory+"distribution_flag_data_"


    if plot_distributions_observations:
        cb_title = 'Occurrence\n(# of 6-minute observations)'
        if universal_time:
            print("###### Plotting universal time vs. frequency distributions for analogic beams ######")
            plot_distribution_freq(ffreq, distribution_observations_data_universaltime_all_ana_beam, type='universal time', shading=shading_distrib_freq, title= 'Frequency vs. universal time distribution of the Observations\n analogic beams', save_image=True, filename = name_pdf+"universaltime_observations_all_ana_beam", cb_title = cb_title)
            if plot_histograms:
                plot_distribution_with_histograms(ffreq,distribution_observations_data_universaltime_all_ana_beam, bins_distribution_observation_anabeam_UT, distribution_observation_anabeam_UT, bins_distribution_observation_numebeam_frequency,distribution_observation_numebeam_frequency,shading_distrib_2D= "auto",figsize = (10, 10),save_image= True,type_x = 'UT',type_y = 'frequency', filename = name_pdf+"universaltime_observations_all_ana_beam_with_histograms_observations")
        if local_time:
            print("###### Plotting local time vs. frequency distributions for analogic beams ######")
            plot_distribution_freq(ffreq, distribution_observations_data_localtime_all_ana_beam, type='local time', shading=shading_distrib_freq, title= 'Frequency vs. local time distribution of the RFI\n analogic beams', save_image=True, filename = name_pdf+"localtime_observations_all_ana_beam", cb_title = cb_title)
            if plot_histograms:
                print("###### Plotting local time vs. frequency distributions for analogic beams with histograms ######")
                plot_distribution_with_histograms(ffreq,distribution_observations_data_localtime_all_ana_beam, bins_distribution_observation_anabeam_LT, distribution_observation_anabeam_LT, bins_distribution_observation_numebeam_frequency,distribution_observation_numebeam_frequency,shading_distrib_2D= "auto",figsize = (10, 10),save_image= True,type_x = 'LT',type_y = 'frequency', filename = name_pdf+"localtime_observations_all_ana_beam_with_histograms_observations")
        if azimuth:
            print("###### Plotting azimuth vs. frequency distributions for analogic beams ######")
            plot_distribution_freq(ffreq, distribution_observations_data_azimuth_all_ana_beam, type='azimuth', shading=shading_distrib_freq, cmap = cmap, title= 'Frequency vs. azimuth distribution of the RFI\n analogic beams', save_image=True, filename = name_pdf+"azimuth_observations_all_ana_beam", cb_title = cb_title)
            print("###### Plotting azimuth vs. frequency distributions for numeric beam ######")
            plot_distribution_freq(ffreq, distribution_observations_data_azimuth_all_nume_beam, type='azimuth', shading=shading_distrib_freq, cmap = cmap, title= 'Frequency vs. azimuth distribution of the RFI\nall numeric beams', save_image=True, filename = name_pdf+"azimuth_observations_all_nume_beam", cb_title = cb_title)
            if plot_histograms:
                print("###### Plotting azimuth vs. frequency distributions for numeric beam with histograms ######")
                plot_distribution_with_histograms(ffreq,distribution_observations_data_azimuth_all_nume_beam, bins_distribution_observation_numebeam_azimuth, distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_frequency,distribution_observation_numebeam_frequency,shading_distrib_2D= "auto",figsize = (10, 10),save_image= True,type_x = 'azimuth',type_y = 'frequency', filename = name_pdf+"azimuth_observations_all_nume_beam_with_histograms_observations")
        if elevation:
            print("###### Plotting elevation vs. frequency distributions for analogic beams ######")
            plot_distribution_freq(ffreq, distribution_observations_data_elevation_all_ana_beam, type='elevation', shading=shading_distrib_freq, cmap = cmap, title= 'Frequency vs. elevation distribution of the RFI\n analogic beams', save_image=True, filename = name_pdf+"elevation_observations_all_ana_beam", cb_title = cb_title)
            print("###### Plotting elevation vs. frequency distributions for numeric beam ######")
            plot_distribution_freq(ffreq, distribution_observations_data_elevation_all_nume_beam, type='elevation', shading=shading_distrib_freq, cmap = cmap, title= 'Frequency vs. elevation distribution of the RFI\nall numeric beams', save_image=True, filename = name_pdf+"elevation_observations_all_nume_beam", cb_title = cb_title)
            if plot_histograms:
                print("###### Plotting elevation vs. frequency distributions for numeric beam with histograms ######")
                plot_distribution_with_histograms(ffreq,distribution_observations_data_elevation_all_nume_beam, bins_distribution_observation_numebeam_elevation, distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_frequency,distribution_observation_numebeam_frequency,shading_distrib_2D= "auto",figsize = (10, 10),save_image= True,type_x = 'elevation',type_y = 'frequency', filename = name_pdf+"elevation_observations_all_nume_beam_with_histograms_observations")
                                  
        #if elaz:
        #    shading_distrib_elaz = "None"
        #    print("###### Plotting azimuth vs. elevation distributions for numeric beam ######")
        #    plot_distribution_elaz(distribution_flag_data_elaz_all_nume_beam[:,:,:], ffreq[:], save_image=True, filename = output_directory+"distribution_flag_elaz/distribution_flag_data_elaz_rfilevel"+str(rfilevel)+"_all_nume_beam", rfilevel=rfilevel, cmap = cmap, shading = shading_distrib_elaz)
        #    filename_png_files = "distribution_flag_data_elaz_rfilevel"+str(rfilevel)+"_all_nume_beam"
        #    print("###### Making GIF for elevation vs. frequency distributions ######")
        #    make_gif(directory_path = 'distribution_flag_elaz/', files_name=f'{filename_png_files}*.png', output_path=f'distribution_flag_elaz/{filename_png_files}_all_freq.gif', fps=60, delete_pngs=False)

    if plot_histograms:
        print("###### Plotting Time distribution of 6-minute slices observations ######")
        plot_1D_time_distribution_observation(time_datetime, save_image = True, filename = output_directory+'histogram_observations_versus_time')
        print("###### Plotting Azimuth distribution of 6-minute slices observations ######")
        plot_1D_distribution_observation(distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_azimuth, type='azimuth', save_image = True, filename = output_directory+'histogram_observations_versus_azimuth_numebeam')
        print("###### Plotting Elevation distribution of 6-minute slices observations ######")
        plot_1D_distribution_observation(distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_elevation, type='elevation', save_image = True, filename = output_directory+'histogram_observations_versus_elevation_numebeam')
        print("###### Plotting Frequency distribution of 6-minute slices observations ######")
        plot_1D_distribution_observation(distribution_observation_numebeam_frequency, bins_distribution_observation_numebeam_frequency, type='frequency', save_image = True, filename = output_directory+'histogram_observations_versus_frequency_numebeam')
        print("###### Plotting Local Time distribution of 6-minute slices observations ######")
        plot_1D_distribution_observation(distribution_observation_anabeam_LT, bins_distribution_observation_anabeam_LT, type='local time', save_image = True, filename = output_directory+'histogram_observations_versus_LT_numebeam', yscale = 'linear')
        print("###### Plotting Universal Time distribution of 6-minute slices observations ######")
        plot_1D_distribution_observation(distribution_observation_anabeam_UT, bins_distribution_observation_anabeam_UT, type='universal time', save_image = True, filename = output_directory+'histogram_observations_versus_UT_numebeam', yscale = 'linear')


    if plot_distributions_rfi:

        print("###### Plotting RFI time vs. frequency distributions for beam 0######")
        new_ffreq, rebinned_distribution = rebin_distribution_in_frequency(ffreq, flag, 0.1953125)
        plot_flag_freq_versus_time_distribution_with_histograms(time_datetime,
                                                         new_ffreq,
                                                         rebinned_distribution[0,:,:],
                                                         rfilevel = rfilevel,
                                                         save_image = True,
                                                         filename = name_pdf+f'RFI_freq_versus_time_with_both_histo_rfilevel{rfilevel}',
                                                         title = f'Frequency vs. Time distribution'+'\n'+f'of the RFI - level {rfilevel}',
                                                         cb_title = r'Flagged data (%)',
                                                         cmap = 'Spectral_r',
                                                         shading=shading_distrib_freq)

        if universal_time:
            print("###### Plotting universal time vs. frequency distributions for analogic beams ######")
            plot_distribution_freq(ffreq, distribution_flag_data_universaltime_all_ana_beam, type='universal time', shading=shading_distrib_freq, title= 'Frequency vs. universal time distribution of the RFI\n analogic beams', save_image=True, filename = name_pdf+"universaltime_rfilevel"+str(rfilevel)+"_all_ana_beam")
            if plot_histograms:
                print("###### Plotting universal time vs. frequency distributions for analogic beams with histograms ######")
                plot_distribution_with_histograms(ffreq,distribution_flag_data_universaltime_all_ana_beam, bins_distribution_observation_anabeam_UT, distribution_observation_anabeam_UT, bins_distribution_observation_numebeam_frequency,distribution_observation_numebeam_frequency,shading_distrib_2D= "auto",figsize = (10, 10),save_image= True,type_x = 'UT',type_y = 'frequency', filename = name_pdf+"universaltime_rfilevel"+str(rfilevel)+"_all_ana_beam_with_histograms_observations")
                plot_distribution_with_mean_flag_histo(ffreq, distribution_flag_data_universaltime_all_ana_beam, shading_distrib_2D= "auto",figsize = (10, 10), save_image= True, type_x = 'UT', type_y = 'frequency', filename = name_pdf+"universaltime_rfilevel"+str(rfilevel)+"_all_ana_beam_with_histograms_flag_mean")
                                  
        if local_time:
            print("###### Plotting local time vs. frequency distributions for analogic beams ######")
            plot_distribution_freq(ffreq, distribution_flag_data_localtime_all_ana_beam, type='local time', shading=shading_distrib_freq, title= 'Frequency vs. local time distribution of the RFI\n analogic beams', save_image=True, filename = name_pdf+"localtime_rfilevel"+str(rfilevel)+"_all_ana_beam")
            if plot_histograms:
                print("###### Plotting local time vs. frequency distributions for analogic beams with histograms ######")
                plot_distribution_with_histograms(ffreq,distribution_flag_data_localtime_all_ana_beam, bins_distribution_observation_anabeam_LT, distribution_observation_anabeam_LT, bins_distribution_observation_numebeam_frequency,distribution_observation_numebeam_frequency,shading_distrib_2D= "auto",figsize = (10, 10),save_image= True,type_x = 'LT',type_y = 'frequency', filename = name_pdf+"localtime_rfilevel"+str(rfilevel)+"_all_ana_beam_with_histograms_observations")
                plot_distribution_with_mean_flag_histo(ffreq, distribution_flag_data_localtime_all_ana_beam, shading_distrib_2D= "auto",figsize = (10, 10), save_image= True, type_x = 'LT',type_y = 'frequency', filename = name_pdf+"localtime_rfilevel"+str(rfilevel)+"_all_ana_beam_with_histograms_flag_mean")
        if azimuth:
            print("###### Plotting azimuth vs. frequency distributions for analogic beams ######")
            plot_distribution_freq(ffreq, distribution_flag_data_azimuth_all_ana_beam, type='azimuth', shading=shading_distrib_freq, cmap = cmap, title= 'Frequency vs. azimuth distribution of the RFI\n analogic beams', save_image=True, filename = name_pdf+"azimuth_rfilevel"+str(rfilevel)+"_all_ana_beam")
            print("###### Plotting azimuth vs. frequency distributions for numeric beam ######")
            plot_distribution_freq(ffreq, distribution_flag_data_azimuth_all_nume_beam, type='azimuth', shading=shading_distrib_freq, cmap = cmap, title= 'Frequency vs. azimuth distribution of the RFI\nall numeric beams', save_image=True, filename = name_pdf+"azimuth_rfilevel"+str(rfilevel)+"_all_nume_beam")
            if plot_histograms:
                print("###### Plotting azimuth vs. frequency distributions for numeric beam with histograms ######")
                plot_distribution_with_histograms(ffreq,distribution_flag_data_azimuth_all_nume_beam, bins_distribution_observation_numebeam_azimuth, distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_frequency,distribution_observation_numebeam_frequency,shading_distrib_2D= "auto",figsize = (10, 10),save_image= True,type_x = 'azimuth',type_y = 'frequency', filename = name_pdf+"azimuth_rfilevel"+str(rfilevel)+"_all_nume_beam_with_histograms_observations")
                plot_distribution_with_mean_flag_histo(ffreq, distribution_flag_data_azimuth_all_nume_beam, shading_distrib_2D= "auto",figsize = (10, 10), save_image= True, type_x = 'azimuth',type_y = 'frequency', filename = name_pdf+"azimuth_rfilevel"+str(rfilevel)+"_all_nume_beam_with_histograms_flag_mean")
        if elevation:
            print("###### Plotting elevation vs. frequency distributions for analogic beams ######")
            plot_distribution_freq(ffreq, distribution_flag_data_elevation_all_ana_beam, type='elevation', shading=shading_distrib_freq, cmap = cmap, title= 'Frequency vs. elevation distribution of the RFI\n analogic beams', save_image=True, filename = name_pdf+"elevation_rfilevel"+str(rfilevel)+"_all_ana_beam")
            print("###### Plotting elevation vs. frequency distributions for numeric beam ######")
            plot_distribution_freq(ffreq, distribution_flag_data_elevation_all_nume_beam, type='elevation', shading=shading_distrib_freq, cmap = cmap, title= 'Frequency vs. elevation distribution of the RFI\nall numeric beams', save_image=True, filename = name_pdf+"elevation_rfilevel"+str(rfilevel)+"_all_nume_beam")
            if plot_histograms:
                print("###### Plotting elevation vs. frequency distributions for numeric beam wit histograms ######")
                plot_distribution_with_histograms(ffreq,distribution_flag_data_elevation_all_nume_beam, bins_distribution_observation_numebeam_elevation, distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_frequency,distribution_observation_numebeam_frequency,shading_distrib_2D= "auto",figsize = (10, 10),save_image= True,type_x = 'elevation',type_y = 'frequency', filename = name_pdf+"elevation_rfilevel"+str(rfilevel)+"_all_nume_beam_with_histograms_observations")
                plot_distribution_with_mean_flag_histo(ffreq, distribution_flag_data_elevation_all_nume_beam, shading_distrib_2D= "auto",figsize = (10, 10), save_image= True, type_x = 'elevation',type_y = 'frequency', filename = name_pdf+"elevation_rfilevel"+str(rfilevel)+"_all_nume_beam_with_histograms_flag_mean")
                                  
        if elaz_square:
            shading_distrib_elaz = "None"
            print("###### Plotting azimuth vs. elevation distributions for numeric beam in a square plot######")
            plot_distribution_elaz(distribution_flag_data_elaz_all_nume_beam[:,:,:], ffreq[:], save_image=True, filename = f'{output_directory}distribution_flag_elaz/rfilevel{rfilevel}/distribution_flag_data_elaz_rfilevel{rfilevel}_all_nume_beam', rfilevel=rfilevel, cmap = cmap, shading = shading_distrib_elaz)
            filename_png_files = "distribution_flag_data_elaz_rfilevel"+str(rfilevel)+"_all_nume_beam"
            print("###### Making GIF for elevation vs. frequency distributions ######")
            make_gif(directory_path = 'distribution_flag_elaz/rfilevel'+str(rfilevel)+'/', files_name=f'{filename_png_files}*.png', output_path=f'distribution_flag_elaz/rfilevel{rfilevel}/{filename_png_files}_all_freq.gif', fps=60, delete_pngs=False)
        if elaz_polar_projection:
            print("###### Plotting azimuth vs. elevation distributions for numeric beam as a polar projection ######")
            plot_distribution_mollweide(distribution_flag_data_elaz_all_nume_beam[:,:,:]*100, ffreq[:], log_colorbar = False,rfilevel = rfilevel, half_sky = True, projection_type = 'orthview',filename = f'{output_directory}distribution_flag_elaz_orthview/rfilevel{rfilevel}/distribution_flag_data_elaz_orthview_rfilevel{rfilevel}_all_nume_beam', save_image = True)
            filename_png_files = "distribution_flag_data_elaz_orthview_rfilevel"+str(rfilevel)+"_all_nume_beam"
            print("###### Making GIF for elevation vs. frequency distributions ######")
            make_gif(directory_path = 'distribution_flag_elaz_orthview/rfilevel'+str(rfilevel)+'/', files_name=f'{filename_png_files}*.png', output_path=f'distribution_flag_elaz_orthview/rfilevel{rfilevel}/{filename_png_files}_all_freq.gif', fps=60, delete_pngs=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Plot RFI flag distribution of NenuFAR observations")
    parser.add_argument('--input_sav_file_directory', dest = 'intput_sav_file_directory', required = False, type = str, default="./", help = "Input sav file directory location to calculate the distributions (default: './'). Useful if calculate_distributions keyword is True.")
    parser.add_argument('--input_hdf5_file_directory', dest = 'input_hdf5_file_directory', required = False, type = str, default = './', help = "Input hdf5 file directory location to calculate the plot (default: './'). Useful if plot_distributions keyword is True.")
    parser.add_argument('-rfi_level', dest = 'rfilevel', required = True, type = int, help = "Level of RFI mitigation (required).")
    parser.add_argument('--ouput_directory', dest = 'output_directory', required = False, type = str, default = './', help = "Output directory (default: './').")
    parser.add_argument('--color_map', dest = 'cmap', required = False, type = str, default = 'viridis', help = "Color map for plots (default: 'viridis').")
    parser.add_argument('--calculate_distributions', dest = 'calculate_distributions', default = False, action = 'store_true',
                        help = 'Keyword to calculate and save the distributions (default: False). If False, and plot_distributions keyword is True, the plotting function will try to load pre-calculated distribution from hdf5 files')
    parser.add_argument('--plot_distributions_rfi', dest = 'plot_distributions_rfi', default = False, action = 'store_true',
                        help = "Keyword to plot the RFI distributions (default: False).")
    parser.add_argument('--plot_distributions_observations', dest = 'plot_distributions_observations', default = False, action = 'store_true',
                        help = "Keyword to plot the Observations 2D distributions (default: False).")
    parser.add_argument('--plot_histograms', dest = 'plot_histograms', default = False, action = 'store_true',
                        help = "Keyword to plot histograms (defaul: False).")
    parser.add_argument('--universal_time', dest = 'universal_time', default = False, action = 'store_true',
                        help = 'Keyword to plot only on the universal time (default: False).')
    parser.add_argument('--local_time', dest = 'local_time', default = False, action = 'store_true',
                        help = 'Keyword to plot only on the local time (default: False).')
    parser.add_argument('--azimuth', dest = 'azimuth', default = False, action = 'store_true',
                        help = 'Keyword to plot only on the azimuth (default: False).')
    parser.add_argument('--elevation', dest = 'elevation', default = False, action = 'store_true',
                        help = 'Keyword to plot only on the elevation (default: False).')
    parser.add_argument('--elaz_square', dest = 'elaz_square', default = False, action = 'store_true',
                        help = 'Keyword to plot only on the elevation vs. azimuth plane, as a square plot (default: False).')
    parser.add_argument('--elaz_polar_projection', dest = 'elaz_polar_projection', default = False, action = 'store_true',
                        help = 'Keyword to plot only on the elevation vs. asimuth, as a polar plot (default: False)')
    parser.add_argument('--save_hdf5', dest = 'save_hdf5', default = False, action = 'store_true',
                        help = 'Keyword to save the calcualted distribution into a hdf5 file (default: False).')
    args = parser.parse_args()

    if args.calculate_distributions == args.plot_distributions_rfi == args.plot_distributions_observations == args.plot_histograms == False:
        raise RuntimeError("Failed to provide any task to perform. Exiting. Please use at least '--calculate_distributions' or '--plot_distributions' keywords")
    
    if args.calculate_distributions:
        print("#### Calculating distributions for RFI mitigation level "+str(args.rfilevel)+" ####")
        (rfilevel, time_datetime, ffreq, flag,
         distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_azimuth,
         distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_elevation,
         distribution_observation_numebeam_frequency, bins_distribution_observation_numebeam_frequency,
         distribution_observation_anabeam_UT, bins_distribution_observation_anabeam_UT,
         distribution_observation_anabeam_LT, bins_distribution_observation_anabeam_LT,
         distribution_flag_data_universaltime_all_ana_beam, distribution_observations_data_universaltime_all_ana_beam,
         distribution_flag_data_localtime_all_ana_beam,distribution_observations_data_localtime_all_ana_beam,
         distribution_flag_data_azimuth_all_ana_beam,distribution_observations_data_azimuth_all_ana_beam,
         distribution_flag_data_elevation_all_ana_beam,distribution_observations_data_elevation_all_ana_beam,
         distribution_flag_data_elevation_all_nume_beam,distribution_observations_data_elevation_all_nume_beam,
         distribution_flag_data_azimuth_all_nume_beam,distribution_observations_data_azimuth_all_nume_beam,
         distribution_flag_data_elaz_all_nume_beam,distribution_observations_data_elaz_all_nume_beam
         ) = calculate_rfi_stat_distribution(args.rfilevel,
                                             args.intput_sav_file_directory,
                                             args.output_directory,
                                             args.save_hdf5
                                            )
        
    if args.calculate_distributions == False:
        print("#### Restoring distributions from hdf5 file ####")
        (time_datetime, ffreq, flag,
         distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_azimuth,
         distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_elevation,
         distribution_observation_numebeam_frequency, bins_distribution_observation_numebeam_frequency,
         distribution_observation_anabeam_UT, bins_distribution_observation_anabeam_UT,
         distribution_observation_anabeam_LT, bins_distribution_observation_anabeam_LT,
         distribution_flag_data_universaltime_all_ana_beam, distribution_observations_data_universaltime_all_ana_beam,
         distribution_flag_data_localtime_all_ana_beam,distribution_observations_data_localtime_all_ana_beam,
         distribution_flag_data_azimuth_all_ana_beam,distribution_observations_data_azimuth_all_ana_beam,
         distribution_flag_data_elevation_all_ana_beam,distribution_observations_data_elevation_all_ana_beam,
         distribution_flag_data_elevation_all_nume_beam,distribution_observations_data_elevation_all_nume_beam,
         distribution_flag_data_azimuth_all_nume_beam,distribution_observations_data_azimuth_all_nume_beam,
         distribution_flag_data_elaz_all_nume_beam,distribution_observations_data_elaz_all_nume_beam
         ) = read_hdf5_file(args.input_hdf5_file_directory, args.rfilevel)

    if args.plot_distributions_rfi or args.plot_histograms:
        print("#### Plotting distributions ####")
        plot_rfi_stat(args.rfilevel,
                      args.output_directory,
                      args.universal_time,
                      args.local_time,
                      args.azimuth,
                      args.elevation,
                      args.elaz_square,
                      args.elaz_polar_projection,
                      args.plot_distributions_rfi,
                      args.plot_distributions_observations,
                      args.plot_histograms,
                      time_datetime,
                      ffreq,
                      flag,
                      distribution_observation_numebeam_azimuth, bins_distribution_observation_numebeam_azimuth,
                      distribution_observation_numebeam_elevation, bins_distribution_observation_numebeam_elevation,
                      distribution_observation_numebeam_frequency, bins_distribution_observation_numebeam_frequency,
                      distribution_observation_anabeam_UT, bins_distribution_observation_anabeam_UT,
                      distribution_observation_anabeam_LT, bins_distribution_observation_anabeam_LT,
                      distribution_flag_data_universaltime_all_ana_beam, distribution_observations_data_universaltime_all_ana_beam,
                      distribution_flag_data_localtime_all_ana_beam,distribution_observations_data_localtime_all_ana_beam,
                      distribution_flag_data_azimuth_all_ana_beam,distribution_observations_data_azimuth_all_ana_beam,
                      distribution_flag_data_elevation_all_ana_beam,distribution_observations_data_elevation_all_ana_beam,
                      distribution_flag_data_elevation_all_nume_beam,distribution_observations_data_elevation_all_nume_beam,
                      distribution_flag_data_azimuth_all_nume_beam,distribution_observations_data_azimuth_all_nume_beam,
                      distribution_flag_data_elaz_all_nume_beam,distribution_observations_data_elaz_all_nume_beam,
                      args.cmap
                     )
        
    print("#### Done ####")
        