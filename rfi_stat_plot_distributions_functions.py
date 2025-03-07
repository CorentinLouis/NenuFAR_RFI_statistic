# -*- coding: utf-8 -*-

"""
RFI_STAT_PLOT_DISTRIBUTIONS_FUNCTIONS

Author: Corentin K. Louis
Date: 2023-09-21
Description: This module contains useful functions for plotting RFI flag of number of observations distribution from NenuFAR multi-beam observations, such as:
                - plot_distribution_freq plots, for the analogic beam, the:
                    - frequency vs. local time
                    - frequency vs. azimuth
                    - frequency vs. elevation
                distributions
                - plot_distribution_elaz plots, for the numeric beams, the:
                    - elevation vs. azimuth
                distribution
"""

from tkinter import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import get_cmap
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter,
                               AutoMinorLocator)
import healpy as hp
import numpy
import datetime
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from mpl_toolkits import axes_grid1
#from scipy.stats import binned_statistic_2d


def plot_distribution_freq(
    freqs: numpy.ndarray, 
    D: numpy.ndarray, 
    type: str = 'azimuth', 
    save_image: bool = False, 
    filename: str = "image", 
    title: str = None, 
    cb_title: str = 'Flagged data (%)',
    fig = None,
    ax = None,
    **kwargs
    ) -> tuple:
    """
    Create a pcolormesh plot to visualize the distribution of values.

    Parameters:
        freqs (numpy.ndarray): Array of frequency values.
        D (numpy.ndarray): Array of data values, shape [n_ephem, n_freq].
        type (str): Type of data ('azimuth', 'elevation', 'UT' or 'universal time', or 'LT' or 'local time').
        save_image (bool): Whether to save the plot as an image file (default is False).
        filename (str): Filename for the saved image (default is "image").
        title (str): Title for the plot (default is an empty string).
        **kwargs: Additional keyword arguments to pass to pcolormesh.

    Returns:
        tuple: A tuple containing:
            - fig (matplotlib.figure.Figure): The created figure.
            - ax (matplotlib.axes.Axes): The plot's axes.
            - cb (matplotlib.colorbar.Colorbar): The colorbar associated with the plot.
    """
    fs_labels = 15
    fs_ticks = 12
    
    if fig == None and ax == None:
        fig_, ax_ = plt.subplots(figsize=(5, 5))
    else:
        fig_ = fig
        ax_ = ax
    
    if type == 'elevation':
        bins = numpy.zeros((91, freqs.size))
        ffreqs = numpy.zeros((91, freqs.size))
        for ifreq in range(freqs.size):
            bins[:,ifreq] = numpy.linspace(0, 90, 91)
            ffreqs[:,ifreq] = freqs[ifreq]
        xtitle = r'Elevation (°)'
    elif type == 'azimuth':
        bins = numpy.zeros((360, freqs.size))
        ffreqs = numpy.zeros((360, freqs.size))
        for ifreq in range(freqs.size):
            bins[:, ifreq] = numpy.linspace(0, 359, 360)
            ffreqs[:, ifreq] = freqs[ifreq]
        xtitle = r'Azimuth (°)'
    elif type == 'LT' or type == 'local time' or type == 'UT' or type == 'universal time':
        bins = numpy.zeros((24, freqs.size))
        ffreqs = numpy.zeros((24, freqs.size))
        for ifreq in range(freqs.size):
            bins[:, ifreq] = numpy.linspace(0, 23, 24)
            ffreqs[:, ifreq] = freqs[ifreq]
        if type == 'LT' or type == 'local time':
            xtitle = r'Local Time (hours)'
        elif type == 'UT' or type == 'universal time':
            xtitle = r'UT (hours)'

    if numpy.max(D[~numpy.isnan(D)]) <= 100:
        vmin = 1
        vmax = 100
        scaleZ = colors.Normalize(vmin=vmin, vmax=vmax)
        im = ax_.pcolormesh(bins, ffreqs, D*100, norm=scaleZ, **kwargs)
        if title == None:
            main_title = f'Frequency vs. {type}'+'\n'+'distribution of the RFI'
            ax_.set_title(main_title, fontsize=fs_labels+2)
    else:
        vmin = 0.01
        vmax = numpy.max(D[~numpy.isnan(D)])
        scaleZ = colors.LogNorm(vmin=vmin, vmax=vmax)
        im = ax_.pcolormesh(bins, ffreqs, D, norm=scaleZ, **kwargs)
        if title == None:
            main_title = f'Frequency vs. {type}'+'\n'+'distribution of the Observations'
            ax_.set_title(main_title, fontsize=fs_labels+2)

    # AXIS
    #ax.xaxis.label.set_size(fs_labels)
    #ax.yaxis.label.set_size(fs_labels)  
    if title != None:
        ax_.set_title(title, fontsize=fs_labels+2)

    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.)
    ax2pos = ax_.get_position()
    cax = fig_.add_axes([ax2pos.x1+0.02, ax2pos.y0, 0.02, ax2pos.height])

    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    ax_.xaxis.set_minor_locator(MultipleLocator(1))
    ax_.yaxis.set_minor_locator(MultipleLocator(1))
    ax_.tick_params(labelsize=fs_ticks)

    # CBAR
    cb = fig_.colorbar(im, extend='both', shrink=0.9, cax=cax, ax=ax_)
    cb.ax.tick_params(labelsize=fs_ticks)
    cb.set_label(cb_title, fontsize=fs_labels)
    
    if fig == None and ax == None:
        ax_.set_xlabel(xtitle, fontsize=fs_labels)
        ax_.set_ylabel(r'Frequency (MHz)', fontsize=fs_labels)
        #fig_.tight_layout()
    
    if save_image:
        plt.savefig(filename + '.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    return fig_, ax_, cb




def plot_distribution_elaz(
        D: numpy.ndarray,
        freqs: numpy.ndarray,
        save_image: bool = False,
        log_colorbar: bool = True,
        filename: str = "image",
        rfilevel: int=0,
        fig = None,
        ax = None,
        **kwargs
) -> tuple:
    """
    Create a pcolormesh plot to visualize the distribution of values.

    Parameters:
        D (numpy.ndarray): Array of data values, shape [n_az, n_el].
        freqs (numpy.ndarray): Array of frequency values [n_freq].
        save_image (bool): Whether to save the plot as an image file (default is False).
        filename (str): Filename for the saved image (default is "image").
        **kwargs: Additional keyword arguments to pass to pcolormesh.

    Returns:
        tuple: A tuple containing:
            - fig (matplotlib.figure.Figure): The created figure.
            - ax (matplotlib.axes.Axes): The plot's axes.
            - cb (matplotlib.colorbar.Colorbar): The colorbar associated with the plot.
    """
  
    fs_labels = 15
    fs_ticks = 12


    bins_el_ = numpy.linspace(0, 90, 91)
    bins_az_ = numpy.linspace(0, 359, 360)
    bins_az, bins_el = numpy.meshgrid(bins_az_,bins_el_)

    xtitle = r'Azimuth (°)'
    ytitle = r'Elevation (°)'
    
    if log_colorbar == False:
        vmin = 1
        vmax = 100
        scaleZ = colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        vmin = 10
        vmax = 100
        scaleZ = colors.LogNorm(vmin=vmin, vmax=vmax)

    if fig == None and ax == None:
        fig_, ax_ = plt.subplots(figsize=(5, 5))
    else:
        fig_ = fig
        ax_ = ax

    for i_freq, freq in enumerate(freqs):
        im = ax_.pcolormesh(bins_el, bins_az, D[:, :, i_freq] * 100, norm=scaleZ, **kwargs)

        #To test: plot with hexagonal bins (ax.hexbin())
        # AXIS
        ax_.set_xlabel(xtitle, fontsize=fs_labels)
        ax_.set_ylabel(ytitle, fontsize=fs_labels)

        title_text = f'Azimuth vs. elevation distribution\nof the RFI (level {rfilevel})\nFrequency {freq}'
        ax_.set_title(title_text, fontsize=fs_labels+2)

        divider = make_axes_locatable(ax_)
        cax = divider.append_axes("right", size=0.15, pad=0.2)

        plt.xticks(fontsize=40)
        plt.yticks(fontsize=40)

        # CBAR
        cb = fig_.colorbar(im, extend='both', shrink=0.9, cax=cax, ax=ax_)
        cb.ax.tick_params(labelsize=fs_ticks)
        cb.set_label(r'Flagged data (%)', fontsize=fs_labels)

        #if fig == None and ax == None:
            #fig_.tight_layout()

        if save_image:
            plt.savefig(f"{filename}_{freq}.png", dpi=300, bbox_inches='tight')
            plt.close()


    return fig_, ax_, cb


def plot_distribution_mollweide(
        D: numpy.ndarray,
        freqs: numpy.ndarray,
        vmin: float = 1,
        vmax: float  = 100,
        unit: str = 'Flagged data (%)',
        log_colorbar: str = None, #option is None, 'log', or 'hist'
        save_image: bool = False,
        projection_type: str = 'mollview',
        filename: str = "image",
        cmap: str = 'YlGnBu_r',
        title: str = None,
        rfilevel: int = 0,
        half_sky = True,
        center_projection = (180, 90, 0),
        **kwargs
) -> tuple:
    """
    Create a Mollweide projection plot to visualize the distribution of values.

    Parameters:
        D (numpy.ndarray): Array of data values, shape [n_az, n_el].
        freqs (numpy.ndarray): Array of frequency values [n_freq].
        save_image (bool): Whether to save the plot as an image file (default is False).
        filename (str): Filename for the saved image (default is "image").
        **kwargs: Additional keyword arguments to pass to hp.mollview.

    Returns:
        tuple: A tuple containing:
            - fig (matplotlib.figure.Figure): The created figure.
            - ax (matplotlib.axes.Axes): The plot's axes.
            - cb (matplotlib.colorbar.Colorbar): The colorbar associated with the plot.
    """

    fs_labels = 15
    fs_ticks = 12

    # Calculate the appropriate nside value for HEALPix
    #nside = 2 ** int(np.ceil(np.log2(max(n_az, n_el) / 3)))
    nside = 64
    # Create an empty HEALPix map
    n_pixels = hp.nside2npix(nside)
    hpx_map = numpy.zeros(n_pixels)

    # create as many plot as elements in freqs:
    for ifreq, freq in enumerate(freqs):
        # Flatten the D array for this frequency
        data = D[:, :, ifreq].flatten()  # Assuming the first frequency for demonstration
    
    
        bins_el_ = numpy.linspace(90, 0, 91)
        bins_az_ = numpy.linspace(0, 359, 360)

        azimuth, elevation  = numpy.meshgrid(bins_az_, bins_el_)
        elevation_rad = numpy.radians(elevation.flatten())
        azimuth_rad = numpy.radians(azimuth.flatten())

        # Assign the data to the HEALPix map
        n_pixels = hp.nside2npix(nside)
        hpx_map = numpy.zeros(n_pixels) + numpy.nan
        healpix_indices = hp.ang2pix(nside, elevation_rad, azimuth_rad, nest = False)

        hpx_map[healpix_indices] = data

        #hpx_map = (hp.get_interp_val(hpx_map, elevation_rad, azimuth_rad, nest = False)

        # Create the Mollweide projection
        fig = plt.figure(figsize=(15, 15))
        ax = fig.add_subplot(111, projection='mollweide')

        # Plot the data on the Mollweide projection
        # Set title
        if title == None:
            title_ = f'Azimuth vs. elevation distribution of the RFI (level {rfilevel})'+'\n'+f'Frequency {freq:.4f} MHz'
        if projection_type == 'mollview':
            hp.mollview(hpx_map, cmap=cmap, hold=True, nest=False, min = vmin, max = vmax, cbar=False, flip = 'astro', rot=center_projection, title = title_)
        if projection_type == 'orthview':
            hp.orthview(hpx_map, cmap=cmap, hold=True, nest=False, 
                        min = vmin, max = vmax,
                        cbar=True, unit = unit, norm = log_colorbar, flip = 'astro', rot=center_projection, half_sky=half_sky, title = title_)

        # Create a custom vertical colorbar
        #cax = inset_axes(
        #    ax,
        #    width='3%',
        #    height='100%',
        #    loc='lower left',
        #    bbox_to_anchor=(1.05, 0., 1, 1),
        #    bbox_transform=ax.transAxes,
        #    borderpad=0,
        #)
        #cb = ColorbarBase(
        #    cax,
        #    cmap=get_cmap(name=cmap),  # You can customize the colormap here
        #    orientation='vertical',
        #    norm=Normalize(vmin=vmin, vmax=vmax),
        #    ticks=numpy.linspace(vmin, vmax, 5),  # Adjust the number of ticks as needed
        #    format = "%d"
        #)
        #cb.solids.set_edgecolor("face")
        #colorbar_label = r'Unflagged data (%)'
        #cb.set_label(colorbar_label, fontsize=fs_labels)
        ## Increase the tick label size
        #cb.ax.tick_params(labelsize=fs_ticks)
        ##cb.formatter.set_powerlimits((0, 0))
        
        #create grid
        hp.graticule(dmer = 30, dpar = 10)

        # Add elevation and azimuth annotations
        for i_az, i_va, i_ha in zip([60, 120, 240, 300], ['bottom', 'top', 'top', 'bottom'], ['left','left', 'right', 'right']):
            el_rad = numpy.radians(0)
            az_rad = numpy.radians(i_az)
            # Convert to HEALPix coordinates
            theta, phi = numpy.pi / 2 - el_rad, az_rad
            hp.projtext(theta, phi, f'{i_az}°', lonlat=False, fontsize=fs_labels, ha=i_ha, va = i_va)
        for i_nesw, i_az, i_va, i_ha in zip(numpy.array(["North", "East", "South", "West"]), numpy.array([0, 90, 180, 270]), ['top', 'center', 'top', 'center'], ['center','left', 'center', 'right']):
            el_rad = numpy.radians(0)
            az_rad = numpy.radians(i_az)
            # Convert to HEALPix coordinates
            theta, phi = numpy.pi / 2 - el_rad, az_rad
            hp.projtext(theta, phi, f'{i_nesw}', lonlat=False, fontsize=fs_labels, ha=i_ha, va=i_va)
        az_label = [180]
        #else:
        #    az_label = [0,359] 
        if half_sky:
            el_plot = [20, 40, 60, 80]
        else:
            el_plot = [-80, -60, -40, -20, 20, 40, 60, 80]
        for az_rad in numpy.radians(az_label):
            for el in el_plot:
                el_rad = numpy.radians(el)
                # Convert to HEALPix coordinates
                theta, phi = numpy.pi / 2 - el_rad, az_rad
                hp.projtext(theta, phi, f' {el}°', lonlat=False, fontsize=fs_labels, ha='center', va = 'center')

        #plt.annotate(title, (0.5, 0.1), xycoords='figure fraction', fontsize=fs_labels)
        #fig.suptitle(title, fontsize = fs_labels+2)
    
        if save_image:
            plt.savefig(f'{filename}_{freq:07.4f}.png', dpi=300, bbox_inches='tight')  
            plt.close()

    return



def plot_1D_time_distribution_observation(
        time_datetime: numpy.array,
        time_min: datetime = datetime.datetime(2019,1,1),
        time_max: datetime = datetime.datetime(2023,12,31),
        save_image: bool = False,
        filename: str = "image",
        fig = None,
        ax = None,
        **kwargs
) -> tuple:
    fs_labels = 15
    fs_ticks = 12
    if fig == None and ax == None:
        fig_, ax_ = plt.subplots(figsize = (10,5))
    else:
        fig_ = fig
        ax_ = ax

    ax_.hist(time_datetime,**kwargs)
    main_title = 'Time distribution of 6-minute slices observations'
    ax_.set_title(main_title, fontsize = fs_labels+2) 
    ax_.set_xlabel(r'Time', fontsize=fs_labels)
    ylabel = r'Occurrence'+'\n'+r'(# of 6-minute observations)'
    ax_.set_ylabel(ylabel, fontsize=fs_labels)
    major_locator=mdates.MonthLocator(interval=6)
    minor_locator=mdates.MonthLocator(interval=1)
    ax_.xaxis.set_major_locator(major_locator)
    ax_.xaxis.set_minor_locator(minor_locator)
    dateFmt = mdates.DateFormatter('%Y\n%B')
    ax_.xaxis.set_major_formatter(dateFmt)
    ax_.tick_params(labelsize=fs_ticks)
    ax_.set_xlim(time_min, time_max)



    #ax.set_yscale('log')
    ax_.yaxis.set_major_formatter(FormatStrFormatter("%d"))
    
    #if fig == None and ax == None:
    #    fig_.tight_layout()

    if save_image:
            plt.savefig(f"{filename}.png", dpi=300, bbox_inches='tight')
            plt.close()

    return fig_, ax_
    

def plot_1D_distribution_observation(
        D: numpy.ndarray,
        bins: numpy.ndarray,
        type: str = None, #options: azimuth, elevation, frequency
        yscale: str = 'log',
        save_image: bool = False,
        filename: str = "image",
        title: str = None,
        fig = None,
        ax = None,
        inverse_axis: bool = False,
        **kwargs
        ) -> tuple:

    if fig == None and ax == None:
        fig_, ax_ = plt.subplots(figsize=(10,5))
    else:
        fig_ = fig
        ax_ = ax
    
    fs_labels = 15
    fs_ticks = 12

    if type == "azimuth":
        main_title = "Azimuth"
        units = r'($\degree$)'
    elif type == 'elevation':
        main_title = "Elevation"
        units = r'($\degree$)'
    elif type == 'local time' or type == 'LT':
        main_title = "Local Time"
        units = r'(hours)'
    elif type == 'universal time' or type == 'UT':
        main_title = "Universal Time"
        units = r'(hours)'
    elif type == 'frequency':
        main_title = 'Frequency'
        units = r'(MHz)'
    elif type =='flag':
        main_title="Flag"
        units = r'(%)'
    else:
        print("Please select type='elevation' or type='azimuth' or type='frequency'")
        print("or type ='UT' or type = 'universal time' or type = 'LT' or type = 'local time'")
        return
    
    if inverse_axis == False:
        ax_.bar(bins[:-1], D, width=numpy.diff(bins), align="edge", **kwargs)
        ax_.set_xlim(numpy.min(bins), numpy.max(bins))
        ax_.set_yscale(yscale)
    else:
        ax_.barh(bins[:-1], D, height=numpy.diff(bins), align="edge", **kwargs)
        #ax.set_ylim(numpy.min(bins), numpy.max(bins))
        ax_.set_xscale(yscale)
    if title == None:
        xtitle = main_title + " " + units 

    ax_.tick_params(labelsize=fs_ticks)
    

    if inverse_axis == False:
        ax_.set_xlabel(xtitle, fontsize=fs_labels)
        ylabel = r'Occurrence'+'\n'+r'(# of 6-minute observations)'
        ax_.set_ylabel(ylabel, fontsize=fs_labels)
        #ax.yaxis.set_major_formatter(FormatStrFormatter("%d"))
        
    else:
        ax_.set_ylabel(xtitle, fontsize=fs_labels)
        ylabel = r'Occurrence'+'\n'+r'(# of 6-minute observations)'
        ax_.set_xlabel(ylabel, fontsize=fs_labels)
        #ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
        
    if fig == None and ax == None:
        main_title = main_title + ' distribution of 6-minute slices observations'
        ax_.set_title(main_title, fontsize=fs_labels+2)
        #fig_.tight_layout()
        

    if save_image:
            plt.savefig(f"{filename}.png", dpi=300, bbox_inches='tight')
            plt.close()
    
    return fig_, ax_


def plot_distribution_with_histograms(ffreq: numpy.array,
                                      distribution_2D: numpy.ndarray,
                                      bins_distribution_1D_histogram_x: numpy.array,
                                      distribution_1D_histogram_x: numpy.ndarray,
                                      bins_distribution_1D_histogram_y: numpy.array,
                                      distribution_1D_histogram_y: numpy.ndarray,
                                      shading_distrib_2D: str = "auto",
                                      figsize = (10, 10),
                                      save_image: bool = False,
                                      filename: str = "image",
                                      cmap: str = 'viridis',
                                      type_1D_histograms: str = r'Occurrence'+'\n'+r'(# of 6-minute observations)',
                                      cb_title = "Flagged data (%)",
                                      type_x: str = None,
                                      type_y: str = None
                                      ) -> tuple:
    
    if type_x == None or type_y == None:
        print("######## *type_x* and *type_y* needed to be defined. Need to correspond to the type of distribution_2D. ########")
        print("######## e.g., type_x = 'elevation' and type_y = 'frequency' ########")
        return
    
    fig_main ,ax_main = plt.subplots(2,2,figsize=figsize)

    distribution_2D = 1-distribution_2D
    
    (fig_main, ax_main[0,1], cb) = plot_distribution_freq(ffreq, distribution_2D, type=type_x, shading=shading_distrib_2D, cmap = cmap, fig = fig_main,ax = ax_main[0,1], cb_title = cb_title)
    (fig_main, ax_main[1,1]) = plot_1D_distribution_observation(distribution_1D_histogram_x, bins_distribution_1D_histogram_x, type=type_x, fig = fig_main,ax = ax_main[1,1])
    (fig_main, ax_main[0,0]) = plot_1D_distribution_observation(distribution_1D_histogram_y, bins_distribution_1D_histogram_y, type=type_y, fig = fig_main,ax = ax_main[0,0], inverse_axis = True)


    ax_main[1,0].axis('off')
    
    ax_main[0][1].sharey(ax_main[0][0])

    ax_main[0,1].get_shared_x_axes().join(ax_main[0,1], ax_main[1,1])
    ax_main[0,1].sharex(ax_main[1,1])

    ax200 = ax_main[0,0].twiny()
    ax200.set_xlim(ax_main[0,0].get_xlim())
    ax200.set_xscale('log')
    ax200.set_xticklabels([])
    ax200.yaxis.set_minor_locator(MultipleLocator(1))
    
    ax211 = ax_main[1,1].twinx()
    ax211.set_ylim(ax_main[1,1].get_ylim())
    ax211.set_yscale('log')
    ax211.set_yticklabels([])
    ax211.xaxis.set_minor_locator(MultipleLocator(1))

    # to avoid having to many minor ticks for Azimuth (range [0,360]):
    if type_x.lower() == "azimuth":
        ax_main[0,1].xaxis.set_minor_locator(MultipleLocator(10))
        ax200.yaxis.set_minor_locator(MultipleLocator(10))
        ax211.xaxis.set_minor_locator(MultipleLocator(10))

    ax_main[0,0].set_xlabel('')
    ax_main[1,1].set_ylabel('')

    fs_labels = 15
    extralabel = type_1D_histograms
    
    ax_main[0,0].annotate(extralabel, (0.1, 0.25), xycoords='figure fraction', fontsize=fs_labels)

    plt.subplots_adjust(wspace=0.1,hspace=0.1)
    
    #fig_main.tight_layout()

    if save_image:
            plt.savefig(f"{filename}.png", dpi=300, bbox_inches='tight')
            plt.close()

    return fig_main, ax_main


def plot_distribution_with_mean_flag_histo(ffreq: numpy.array,
                                           distribution_2D: numpy.ndarray,
                                           shading_distrib_2D: str = "auto",
                                           figsize = (10, 10),
                                           save_image: bool = False,
                                           filename: str = "image",
                                           cmap: str = 'viridis',
                                           type_x: str = None,
                                           type_y: str = None
                                      ) -> tuple:
    
    if type_x == None or type_y == None:
        print("######## *type_x* and *type_y* needed to be defined. Need to correspond to the type of distribution_2D. ########")
        print("######## e.g., type_x = 'elevation' and type_y = 'frequency' ########")
        return
    
    distribution_2D = 1-distribution_2D

    if type_x == 'elevation':
        bins_x = numpy.linspace(0, 90, 91)
        xtitle = r'Elevation (°)'
    elif type_x == 'azimuth':
        bins_x = numpy.linspace(0, 359, 360)
        xtitle = r'Azimuth (°)'
    elif type_x == 'LT' or type_x == 'local time' or type_x == 'UT' or type_x == 'universal time':
        bins_x = numpy.linspace(0, 23, 24)
        if type_x == 'LT' or type_x == 'local time':
            xtitle = r'Local Time (hours)'
        elif type_x == 'UT' or type_x == 'universal time':
            xtitle = r'UT (hours)'
    if type_y == 'frequency':
        bins_y = ffreq
        ytitle = r'Frequency (MHz)'
    elif type_y == 'elevation':
        bins_y = numpy.linspace(0, 90, 91)
        ytitle = r'Elevation (°)'
    
    fs_labels = 15
    fs_ticks = 12

    
    fig_main ,ax_main = plt.subplots(2,2,figsize=figsize)
    (fig_main, ax_main[0,1], cb) = plot_distribution_freq(ffreq, distribution_2D, type=type_x, shading=shading_distrib_2D, cmap = cmap, fig = fig_main,ax = ax_main[0,1])
    ax_main[0,0].barh(bins_y, numpy.nanmean(distribution_2D, axis=0)*100)
    ax_main[0,0].set_ylabel(ytitle, fontsize=fs_labels)
    ax_main[1,1].bar(bins_x, numpy.nanmean(distribution_2D, axis=1)*100)
    ax_main[1,1].set_xlabel(xtitle, fontsize=fs_labels)

    if type_x.lower() == "azimuth":
        ax_main[0,1].xaxis.set_minor_locator(MultipleLocator(10))
        ax_main[1,1].xaxis.set_minor_locator(MultipleLocator(10))
    else:
        ax_main[0,1].xaxis.set_minor_locator(MultipleLocator(1))
        ax_main[1,1].xaxis.set_minor_locator(MultipleLocator(1))

    
    ax_main[0,0].xaxis.set_minor_locator(MultipleLocator(1))
    ax_main[0,0].yaxis.set_minor_locator(MultipleLocator(1))
    ax_main[0,1].yaxis.set_minor_locator(MultipleLocator(1))
    ax_main[1,1].yaxis.set_minor_locator(MultipleLocator(1))        

    ax_main[1,0].axis('off')
    
    ax_main[0][1].sharey(ax_main[0][0])

    ax_main[0,1].get_shared_x_axes().join(ax_main[0,1], ax_main[1,1])
    ax_main[0,1].sharex(ax_main[1,1])

    ax200 = ax_main[0,0].twiny()
    ax200.set_xlim(ax_main[0,0].get_xlim())
    #ax200.set_xscale('log')
    ax200.set_xticklabels([])
    

    
    ax211 = ax_main[1,1].twinx()
    ax211.set_ylim(ax_main[1,1].get_ylim())
    #ax211.set_yscale('log')
    ax211.set_yticklabels([])

    ax_main[0,0].set_xlabel('')
    ax_main[1,1].set_ylabel('')

    ax_main[0,0].tick_params(labelsize=fs_ticks)
    ax_main[1,1].tick_params(labelsize=fs_ticks)

    fs_labels = 15
    extralabel = "Mean RFI flag (%)"
    
    ax_main[0,0].annotate(extralabel, (0.1, 0.25), xycoords='figure fraction', fontsize=fs_labels)

    plt.subplots_adjust(wspace=0.1,hspace=0.1)
   # fig_main.tight_layout()

    if save_image:
            plt.savefig(f"{filename}.png", dpi=300, bbox_inches='tight')
            plt.close()

    return fig_main, ax_main
    

    # Function to identify gaps larger than a threshold and insert datetime elements
def insert_gaps(time_array, data_array, gap_threshold_minutes=30):
    new_time_array = []
    new_data_array = []

    for itime in range(1, len(time_array)):
        time_diff = (time_array[itime] - time_array[itime - 1]).total_seconds() / 60  # Time difference in minutes
        if time_diff > gap_threshold_minutes:
            # Insert new datetime element and row of numpy.nan in data array 6 minutes after the start of the gap
            new_time_start = time_array[itime - 1] + datetime.timedelta(minutes=6)
            new_data_row_start = numpy.full(data_array.shape[1], numpy.nan)
            
            new_time_array.append(new_time_start)
            new_data_array.append(new_data_row_start)

            # Insert new datetime element and row of numpy.nan in data array 6 minutes before the end of the gap
            new_time_end = time_array[itime] - datetime.timedelta(minutes=6)
            new_data_row_end = numpy.full(data_array.shape[1], numpy.nan)

            new_time_array.append(new_time_end)
            new_data_array.append(new_data_row_end)

    # Combine the new elements with the original time and data arrays
    combined_time_array = numpy.concatenate([time_array, new_time_array])
    combined_data_array = numpy.concatenate([data_array, numpy.array(new_data_array)])

    # Sort the combined arrays based on time
    sorted_indices = numpy.argsort(combined_time_array)
    sorted_time_array = combined_time_array[sorted_indices]
    sorted_data_array = combined_data_array[sorted_indices]

    return sorted_time_array, sorted_data_array

def plot_flag_freq_versus_time_distribution_with_histograms(time: numpy.array, # 1D datetime object
                                                             freq: numpy.array, # 1D array of frequencies
                                                             distribution: numpy.array, # 2D array shape (ntime, nfreq)
                                                             rfilevel: int = 0,
                                                             fig = None,
                                                             ax = None,
                                                             title: str = None,
                                                             save_image: bool = True,
                                                             filename: str = 'image', #f'RFI_freq_versus_time_with_both_histo_rfilevel{rfilevel}',
                                                             cb_title: str = r'Flagged data (%)',
                                                             cmap: str = 'Spectral_r',
                                                             time_min = None, # datetime object
                                                             time_max = None, # datetime object
                                                             freq_min: float = None,
                                                             freq_max: float = None,
                                                             vmin: float = 0,
                                                             vmax: float = 30,
                                                             gap_threshold_minutes: float = 4320,
                                                             **kwargs
                                                            ):
    if time_min == None:
        time_min = time[0] 
    if time_max == None:
        time_max = time[-1]
    if freq_min == None:
        freq_min =  freq[0]
    if freq_max == None:
        freq_max =  freq[-1]

    
    if fig == None and ax == None:
        fig_, ax_ = plt.subplots(2, 2, gridspec_kw={'height_ratios': [2, 1], 'width_ratios': [5,12]}, figsize=(10, 10))
    else:
        fig_ = fig
        ax_ = ax

    fs_labels = 15
    fs_ticks = 12

    scaleZ = colors.Normalize(vmin=vmin, vmax=vmax)

    bins_with_nans_inserted, D_with_nans_inserted = insert_gaps(time, distribution, gap_threshold_minutes=4320)


    # Converting datetime to numerical values
    num_time_intervals_with_nans = mdates.date2num(bins_with_nans_inserted)
    
    # Ploting using pcolormesh with regularly spaced numerical data
    im = ax_[0, 1].pcolormesh(num_time_intervals_with_nans, freq, (1-D_with_nans_inserted.T)*100, norm=scaleZ, cmap=cmap, **kwargs) #, shading='nearest'

    # Calculating histograms of the specific (time,freq) window asked by the user
    D_mean_per_freq = numpy.nanmean(distribution[(time>=time_min) & (time<= time_max),:], axis = 0)
    D_mean_per_time = numpy.nanmean(D_with_nans_inserted[:,(freq>=freq_min) & (freq<=freq_max)], axis = 1)

    
    ax_[1,1].plot(num_time_intervals_with_nans, (1-D_mean_per_time)*100, 'k-')
    ax_[0,0].plot((1-D_mean_per_freq)*100, freq, color = 'black')
    ax_[0,0].grid(True)
    ax_[1,1].grid(True)
    
    
    # Time format for axis
    major_locator=mdates.MonthLocator(interval=6)
    minor_locator=mdates.MonthLocator(interval=1)
    ax_[1,1].xaxis.set_major_locator(major_locator)
    ax_[1,1].xaxis.set_minor_locator(minor_locator)
    dateFmt = mdates.DateFormatter('%Y\n%B')
    ax_[1,1].xaxis.set_major_formatter(dateFmt)
    ax_[1,1].tick_params(labelsize=fs_ticks)

    ax_[1,1].yaxis.set_minor_locator(MultipleLocator(5))
    ax_[0,0].yaxis.set_minor_locator(MultipleLocator(1))
    ax_[0,0].xaxis.set_minor_locator(MultipleLocator(5))
    ax_[0,1].yaxis.set_minor_locator(MultipleLocator(1))

    # Axis boundaries
    ax_[1,1].set_ylim(vmin, vmax)
    ax_[1,1].set_xlim(time_min, time_max)
    ax_[0,0].set_ylim(freq_min, freq_max)
    ax_[0,0].set_xlim(vmin, vmax)

    # Modifying axes so that axes are correclty shared between panels
    ax_[1,0].axis('off')

    ax_[0,1].sharey(ax_[0,0])

    ax_[0,1].get_shared_x_axes().join(ax_[0,1], ax_[1,1])
    ax_[0,1].sharex(ax_[1,1])

    ax200 = ax_[0,0].twiny()
    ax200.set_xlim(ax_[0,0].get_xlim())
    ax200.set_xticklabels([])

    ax211 = ax_[1,1].twinx()
    ax211.set_ylim(ax_[1,1].get_ylim())
    ax211.set_yticklabels([])

    ax_[0,1].set_xlabel('')
    ax_[0,1].set_ylabel('')  

    # configuring ticks size
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    ax_[0,0].tick_params(labelsize=fs_ticks)
    #ax_[0,1].set_xticklabels([])
    #ax_[0,1].set_yticklabels([])
    ax_[0,1].tick_params(labelsize=fs_ticks)
    ax_[1,1].tick_params(labelsize=fs_ticks)
    
    
    # configuring color bar
    ax2pos = ax_[0,1].get_position()
    cax = fig_.add_axes([ax2pos.x1+0.02, ax2pos.y0, 0.02, ax2pos.height])

    cb = fig_.colorbar(im, extend='both', shrink=0.9, cax=cax, ax=ax_[0])
    cb.ax.tick_params(labelsize=fs_ticks)
    cb.set_label(cb_title, fontsize=fs_labels)

  
    
    # Adding title and labels

    if title == None:
        main_title = f'Frequency vs. Time distribution'+'\n'+f'of the RFI - level {rfilevel}'
        ax_[0,1].set_title(main_title, fontsize=fs_labels+2)
    else:
        ax_[0,1].set_title(title, fontsize=fs_labels+2)

    ## shared label for histograms
    #extralabel = cb_title     
    #ax_[0,0].annotate(extralabel, (0.12, 0.15), xycoords='figure fraction', fontsize=fs_labels)


    if fig == None and ax == None:
        ax_[1,1].set_xlabel(r'Time', fontsize=fs_labels)
        ylabel = r'Frequency (MHz)'
        ax_[0,0].set_ylabel(ylabel, fontsize=fs_labels)
        ax_[0,0].set_xlabel(cb_title, fontsize=fs_labels)
        ax_[1,1].set_ylabel(cb_title, fontsize=fs_labels)
        plt.subplots_adjust(wspace=0.2,hspace=0.2)

        
#        plt.tight_layout()

    if save_image:
        plt.savefig(filename + '.png', dpi=300, bbox_inches='tight')
        plt.close()
    


    return(fig_, ax_)
