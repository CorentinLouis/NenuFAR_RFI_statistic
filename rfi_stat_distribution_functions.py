# -*- coding: utf-8 -*-

"""
RFI_STAT_DISTRIBUTIONS_FUNCTIONS

Author: Corentin K. Louis
Date: 2023-09-21
Description: This module contains useful functions for calculating RFI flag distribution from NenuFAR multi-beam observations, such as:
                - distribution_flag_ana_beam_all_beam determines, for the analogic beam, the:
                    - frequency vs. local time
                    - frequency vs. azimuth
                    - frequency vs. elevation
                distributions
                - distribution_flag_nume_beam_all_beam determines, for the numeric beams, the:
                    - frequency vs. azimuth
                    - frequency vs. elevation
                distributions
                - distribution_flag_elaz_nume_beam_all_beam determines, for the numeric beams, the:
                    - elevation vs. azimuth
                distribution for all frequencies (one plot per frequency channel)
"""


import numpy

def distribution_flag_ana_beam_all_beam(
    flag: numpy.ndarray,  # [nbeam, nfreq, ntime]
    freqs: numpy.ndarray,  # [nfreq]
    ephemeris: numpy.ndarray,  # [ntime]
    type: str = 'azimuth'  # Choices: 'azimuth', 'elevation', 'LT' or 'local time', 'UT' or 'universal time'
) -> numpy.ndarray:
    """
    Calculate the distribution of flag values based on ephemeris conditions.

    Parameters:
        flag (np.ndarray): Flag mask array per beam.
        freqs (np.ndarray): Frequency values.
        ephemeris (np.ndarray): Ephemeris data (azimuth, elevation, or local time).
        type (str): Type of ephemeris data. Choices: 'azimuth', 'elevation', 'UT' or 'universal time', 'LT', or 'local time'.
    
    Returns:
        np.ndarray: Azimuth, Elevation or Local Time vs. Frequency distribution of flag values using the analogic beam values.
    """
    
    # Calculate the bin edges based on the specified type
    if type == 'elevation':
        bins = numpy.linspace(0, 90, 91)
    elif type == 'azimuth':
        bins = numpy.linspace(0, 359, 360)
    elif type == 'LT' or type == 'local time' or type == 'UT' or type == 'universal time':
        bins = numpy.linspace(0, 23, 24)

    # Calculate the indices where ephemeris conditions are met
    ephemeris_indices = numpy.digitize(ephemeris, bins) - 1  # Subtract 1 to match 0-based indexing

    distribution_flag = numpy.zeros((len(bins), freqs.shape[0])) + numpy.nan
    distribution_observation = numpy.zeros((len(bins), freqs.shape[0]))

    for i_freq in range(len(freqs)):
        for i_index in range(len(bins)):
            ephemeris_cond = ephemeris_indices == i_index
            if ephemeris_cond.any():
                values = flag[:, ephemeris_cond, i_freq]
                ephemeris_mask = values > -1
                if ephemeris_mask.any():
                    distribution_flag[i_index, i_freq] = numpy.mean(values[ephemeris_mask])
                    distribution_observation[i_index, i_freq] = ephemeris_mask.sum() 

    distribution_observation[distribution_observation == 0] = numpy.nan
    return distribution_flag, distribution_observation


def distribution_flag_nume_beam_all_beam(
        flag: numpy.ndarray, # [nbeam, nfreq, ntime]
        freqs: numpy.ndarray, # [nfreq]
        ephemeris: numpy.ndarray, # [nbeam, ntime]
        type: str='azimuth' # Choices: 'azimuth' or 'elevation'
) -> numpy.ndarray:
    """
    Calculate the distribution of flag values based on ephemeris conditions.

    Parameters:
        flag (np.ndarray): Flag mask array per beam.
        freqs (np.ndarray): Frequency values.
        ephemeris (np.ndarray): Ephemeris data (azimuth or elevation per beam).
        type (str): Type of ephemeris data. Choices: 'azimuth' or 'elevation'
    
    Returns:
        np.ndarray: Azimuth or Elevation vs. Frequency distribution of flag values using the numeric beams values.
    """
    # Calculate the bin edges based on the specified type
    if type == 'elevation':
        bins = numpy.linspace(0, 90, 91)
    elif type == 'azimuth':
        bins = numpy.linspace(0, 359, 360)

    # Calculate the indices where ephemeris conditions are met
    ephemeris_indices = numpy.digitize(ephemeris, bins) - 1  # Subtract 1 to match 0-based indexing
    distribution_flag = numpy.zeros((len(bins), freqs.shape[0])) + numpy.nan
    distribution_observation = numpy.zeros((len(bins), freqs.shape[0]))

    for i_freq in range(len(freqs)):
        for i_index in range(len(bins)):
            ephemeris_cond = ephemeris_indices == i_index
            if ephemeris_cond.any():
                values = flag[ephemeris_cond, i_freq]
                ephemeris_mask = values > -1
                if ephemeris_mask.any():
                    distribution_flag[i_index, i_freq] = numpy.mean(values[ephemeris_mask])
                    distribution_observation[i_index, i_freq] = ephemeris_mask.sum() 
    
    distribution_observation[distribution_observation == 0] = numpy.nan
    
    return distribution_flag, distribution_observation


def distribution_flag_elaz_nume_beam_all_beam(
        flag: numpy.ndarray, # [nbeam, nfreq, ntime]
        el: numpy.ndarray, # [nbeam, ntime]
        az: numpy.ndarray, # [nbeam, ntime]
        nfreq: numpy.ndarray # [nfreq]
) -> numpy.ndarray:
    
    """
    Calculate the distribution of flag values based on ephemeris conditions.

    Parameters:
        flag (np.ndarray): Flag mask array per beam.
        el (np.ndarray): Elevation values per numeric beam.
        az (np.ndarray): Azimuth values per numeric beam.
        freqs (np.ndarray): Frequency values.    

    Returns:
        np.ndarray: Elevation vs. Azimuth distribution of flag values per frequencu using the numeric beams values.
    """

    # Calculate the bin edges for elevation and azimuth
    bins_el = numpy.linspace(0, 90, 91)
    bins_az = numpy.linspace(0, 359, 360)
    
    # Calculate the indices where elevation and azimuth conditions are met
    el_indices = numpy.digitize(el, bins_el) - 1  # Subtract 1 to match 0-based indexing
    az_indices = numpy.digitize(az, bins_az) - 1

    distribution_flag_elaz = numpy.zeros((len(bins_el), len(bins_az), nfreq)) + numpy.nan
    distribution_observation_elaz = numpy.zeros((len(bins_el), len(bins_az), nfreq))
    
    
    
    for i_freq in range(nfreq):
        for i_el in range(len(bins_el)):
            for i_az in range(len(bins_az)):
                el_cond = el_indices == i_el
                az_cond = az_indices == i_az
                combined_cond = el_cond & az_cond
                if combined_cond.any():
                    values = flag[combined_cond, i_freq]
                    ephemeris_mask = values > -1
                    if ephemeris_mask.any():
                        distribution_flag_elaz[i_el, i_az, i_freq] = numpy.mean(values[ephemeris_mask])
                        distribution_observation_elaz[i_el, i_az, i_freq] = ephemeris_mask.sum() 

    distribution_observation_elaz[distribution_observation_elaz == 0] = numpy.nan
    
    return distribution_flag_elaz, distribution_observation_elaz


def distribution_observations_anabeam(flag,ephemeris, type="None"):
    if type == 'elevation':
        bins = numpy.linspace(0, 90, 91)
        distribution_obs = numpy.zeros(shape=(len(bins)-1))
    elif type == 'azimuth':
        bins = numpy.linspace(0, 359, 360)
        distribution_obs = numpy.zeros(shape=(len(bins)-1))
    elif type == 'LT' or type == 'local time' or type == 'UT' or type == 'universal time':
        bins = numpy.linspace(0, 23, 24)
        distribution_obs = numpy.zeros(shape=(len(bins)-1))
    elif type == 'frequency':
        freq_min = int(ephemeris[0])
        freq_max = int(ephemeris[-1])
        bins = numpy.linspace(freq_min, freq_max, (freq_max-freq_min)*4+1)
        occ_freq_ = numpy.zeros(shape=(len(ephemeris)))
    else:
        print("Please select type='elevation' or type='azimuth' or type='frequency'")
        print("or type ='UT' or type = 'universal time' or type = 'LT' or type = 'local time'")
        return
    
    if type != 'frequency':
        occurrence_mask = numpy.where(flag[0,:,800] > -1, numpy.ones_like(flag[0,:,800]>-1, dtype=bool), False)
        distribution_obs, bins_histo= numpy.histogram(ephemeris[occurrence_mask], bins=bins)
    else:
        distribution_obs = numpy.zeros(len(bins)-1)
        for i_freq in range(len(ephemeris)):
            occurrence_mask_ = numpy.where(flag[0,:,i_freq] > -1, numpy.ones_like(flag[0,:,i_freq]>-1, dtype=bool), False)
            occ_freq_[i_freq] = len(occurrence_mask_[occurrence_mask_ == True])
        for i_index in range(len(bins)-1):
            distribution_obs[i_index] = numpy.mean(occ_freq_[(ephemeris>=bins[i_index]) & (ephemeris <bins[i_index+1])])
        bins_histo = bins
    return distribution_obs, bins_histo


def distribution_observations_numebeam(flag,ephemeris, type="None"):
    if type == 'elevation':
        bins = numpy.linspace(0, 90, 91)
        distribution_obs = numpy.zeros(shape=(4, len(bins)-1))
    elif type == 'azimuth':
        bins = numpy.linspace(0, 359, 360)
        distribution_obs = numpy.zeros(shape=(4, len(bins)-1))
    elif type == 'frequency':
        freq_min = int(ephemeris[0])
        freq_max = int(ephemeris[-1])
        bins = numpy.linspace(freq_min, freq_max, (freq_max-freq_min)*4+1)
        occ_freq_ = numpy.zeros(shape=(4, len(ephemeris)))
    else:
        print("Please select type='elevation' or type='azimuth' or type='frequency'")
        return
    
    if type != 'frequency':
        for ibeam in range(4):
            occurrence_mask = numpy.where(flag[ibeam,:,800] >= 0, numpy.ones_like(flag[ibeam,:,800]>=0, dtype=bool), False)
            distribution_obs[ibeam,:], bins_histo= numpy.histogram(ephemeris[ibeam, occurrence_mask], bins=bins)
        distribution_obs=numpy.sum(distribution_obs,axis=0)
    else:
        distribution_obs = numpy.zeros(shape=(4, len(bins)-1))
        for ibeam in range(4):
            for i_freq in range(len(ephemeris)):
                occurrence_mask_ = numpy.where(flag[ibeam,:,i_freq] >=0, numpy.ones_like(flag[ibeam,:,i_freq]>=0, dtype=bool), False)
                occ_freq_[ibeam, i_freq] = len(occurrence_mask_[occurrence_mask_ == True])
            for i_index in range(len(bins)-1):
                distribution_obs[ibeam, i_index] = numpy.mean(occ_freq_[ibeam, (ephemeris>=bins[i_index]) & (ephemeris <bins[i_index+1])])
        distribution_obs=numpy.sum(distribution_obs,axis=0)
        bins_histo = bins
    return distribution_obs, bins_histo

######################################

def distribution_flag_ana_beam_per_beam_old(flag, freqs, ephemeris, type='azimuth'):
        # flag: array [nfreq, ntime] containing the flag mask
        # freq: array [nfreq] containing frequency values
        # ephemeris: array [ntime] containing either (depending on type) azimuth, elevation or time informations
        # type= 'azimuth' or 'elevation' or 'LT'
        
        if type == 'elevation':
            bins = numpy.linspace(0,90,91)
        if type == 'azimuth':
            bins = numpy.linspace(0,359,360)
        if type == 'LT' or type == 'local time':
            bins = numpy.linspace(0,23,24)

        distribution_flag = numpy.zeros(shape=(len(bins), freqs.shape[0]))

        
        for i_freq in range(len(freqs)):
            for i_index in range(len(bins)):
                ephemeris_mask = numpy.where(flag[(ephemeris>=bins[i_index]) & (ephemeris < bins[i_index]+1), i_freq] > 0, numpy.ones_like(flag[(ephemeris>=bins[i_index]) & (ephemeris < bins[i_index]+1),i_freq], dtype=bool), False)
                distribution_flag[i_index, i_freq] =  numpy.mean(flag[(ephemeris>=bins[i_index]) & (ephemeris < bins[i_index]+1), i_freq][ephemeris_mask])
        
        return distribution_flag



def distribution_flag_ana_beam_all_beam_old(flag, freqs, ephemeris, type='azimuth'):
        # flag: array [nbeam, nfreq, ntime] containing the flag mask
        # freq: array [nfreq] containing frequency values
        # ephemeris: array [ntime] containing either (depending on type) azimuth, elevation or time informations
        # type= 'azimuth' or 'elevation' or 'LT'
        
        if type == 'elevation':
            bins = numpy.linspace(0,90,91)
        if type == 'azimuth':
            bins = numpy.linspace(0,359,360)
        if type == 'LT' or type == 'local time':
            bins = numpy.linspace(0,23,24)

        distribution_flag = numpy.zeros(shape=(len(bins), freqs.shape[0]))

        
        for i_freq in range(len(freqs)):
            for i_index in range(len(bins)):
                ephemeris_mask = numpy.where(flag[:,(ephemeris>=bins[i_index]) & (ephemeris < bins[i_index]+1), i_freq] > -1, numpy.ones_like(flag[:,(ephemeris>=bins[i_index]) & (ephemeris < bins[i_index]+1), i_freq], dtype=bool), False)
                distribution_flag[i_index, i_freq] =  numpy.mean(flag[:, (ephemeris>=bins[i_index]) & (ephemeris < bins[i_index]+1), i_freq][ephemeris_mask])
        
        return distribution_flag



def distribution_flag_nume_beam_all_beam_old(flag, freqs, ephemeris, type='azimuth'):
        # flag: array [nbeam, nfreq, ntime] containing the flag mask
        # freq: array [nfreq] containing frequency values
        # ephemeris: array [nbeam, ntime] containing either (depending on type) azimuth, elevation or time informations
        # type= 'azimuth' or 'elevation'
        
        if type == 'elevation':
            bins = numpy.linspace(0,90,91)
        if type == 'azimuth':
            bins = numpy.linspace(0,359,360)
        #if type == 'LT' or type == 'local time':
        #    bins = numpy.linspace(0,23,24)

        distribution_flag = numpy.zeros(shape=(len(bins), freqs.shape[0]))

        
        for i_freq in range(len(freqs)):
            for i_index in range(len(bins)):
                ephemeris_mask = numpy.where(flag[(ephemeris>=bins[i_index]) & (ephemeris < bins[i_index]+1), i_freq] > -1, numpy.ones_like(flag[(ephemeris>=bins[i_index]) & (ephemeris < bins[i_index]+1), i_freq], dtype=bool), False)
                distribution_flag[i_index, i_freq] =  numpy.mean(flag[(ephemeris>=bins[i_index]) & (ephemeris < bins[i_index]+1), i_freq][ephemeris_mask])
        
        return distribution_flag



def distribution_flag_elaz_nume_beam_all_beam_old(flag, el, az, nfreq):
    # flag: array [nbeam, nfreq, ntime] containing the flag mask
    # freq: array [nfreq] containing frequency values
    # ephemeris: array [nbeam, ntime] containing either (depending on type) azimuth, elevation or time informations
    # type= 'azimuth' or 'elevation'
        
        
    bins_el = numpy.linspace(0,90,91)
    bins_az = numpy.linspace(0,359,360)
        
    distribution_flag_elaz = numpy.zeros(shape=(len(bins_el), len(bins_az),nfreq))

    for i_freq in range(nfreq):
        print("#### Freq "+str(i_freq)+"/"+str(nfreq-1)+" ####")
        for i_el in range(bins_el.size):
            if i_el % 10 ==0:
                print("#### Elevation "+str(bins_el[i_el])+"("+str(i_el)+"/"+str(bins_el.size)+") ####")
            for i_az in range(bins_az.size):
                #if i_az % 100 ==0:
                    #print("#### Azimuth "+str(bins_az[i_az])+"("+str(i_az)+"/"+str(bins_az.size)+") ####")
                ephemeris_mask = numpy.where(flag[(el>=bins_el[i_el]) & (el < bins_el[i_el]+1) & (az>=bins_az[i_az]) & (az < bins_az[i_az]+1), i_freq] > -1, numpy.ones_like(flag[(el>=bins_el[i_el]) & (el < bins_el[i_el]+1) & (az>=bins_az[i_az]) & (az < bins_az[i_az]+1), i_freq], dtype=bool), False)
                distribution_flag_elaz[i_el, i_az, i_freq] =  numpy.mean(flag[(el>=bins_el[i_el]) & (el < bins_el[i_el]+1) & (az>=bins_az[i_az]) & (az < bins_az[i_az]+1), i_freq][ephemeris_mask])
        
    return distribution_flag_elaz

