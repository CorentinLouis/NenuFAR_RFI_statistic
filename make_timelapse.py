# -*- coding: utf-8 -*-

"""
MAKE_GIF

Author: Corentin K. Louis
Date: 2024-11-08
Description: This module contains a routine that allows to create timelapse from a sequence of PNG files in a specified directory
"""

# Example usage
#  input_directory = "/Users/clouis/Documents/General/Exoradio/RFI_STAT/distribution_flag_elaz_orthview_rebin500KHz/rfilevel3/"
# create_timelapse(input_directory, f'{input_directory}distribution_flag_data_elaz_orthview_rfilevel3_all_nume_beam_all_freq.mp4', framerate=12)

import os
import subprocess


def create_timelapse(input_directory, output_video, framerate=5):
    """
    Renames image files in a directory to a sequential format and creates a video from them using FFmpeg.
    Ensures that the image dimensions are even to avoid encoder errors.
    
    Parameters:
    input_directory (str): Path to the directory containing the image files.
    output_video (str): Path to the output video file (e.g., '/path/to/output_video.mp4').
    framerate (int): The frame rate of the output video (default is 5 fps).
    
    Returns:
    None: The function prints a success message upon completion.
    
    Raises:
    FileNotFoundError: If the input directory does not exist or contains no PNG files.
    OSError: If an error occurs during file renaming or FFmpeg execution.
    """
    # Change to the directory with images
    if not os.path.exists(input_directory):
        raise FileNotFoundError(f"Input directory '{input_directory}' does not exist.")
    
    os.chdir(input_directory)
    
    # Get a list of PNG files in the directory
    files = sorted([f for f in os.listdir(input_directory) if f.endswith('.png')])
    
    if not files:
        raise FileNotFoundError("No PNG files found in the input directory.")
    
    # Rename images to a sequential format
    for i, filename in enumerate(files):
        new_name = f"frame_{i:04d}.png"
        os.rename(filename, new_name)
    
    # FFmpeg command to create the video with an even width and height
    ffmpeg_command = [
        'ffmpeg',
        '-framerate', str(framerate),
        '-i', 'frame_%04d.png',
        '-vf', "pad=ceil(iw/2)*2:ceil(ih/2)*2",  # Ensure width and height are even
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',
        output_video
    ]
    
    # Run the FFmpeg command
    result = subprocess.run(ffmpeg_command, capture_output=True, text=True)
    
    # Check if FFmpeg ran successfully
    if result.returncode == 0:
        print("Video created successfully!")
    else:
        print("An error occurred while creating the video:")
        print(result.stderr)

