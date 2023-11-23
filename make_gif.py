# -*- coding: utf-8 -*-

"""
MAKE_GIF

Author: Corentin K. Louis
Date: 2023-09-21
Description: This module contains a routine that allows to create a GIF animation from a sequence of PNG files in a specified directory
"""


import os
import glob
import imageio

def make_gif(
    directory_path: str = './',
    files_name: str = 'filename*.png',
    output_path: str = 'output.gif',
    fps: int = 10,
    delete_pngs: bool = False
) -> None:
    """
    Create a GIF animation from a sequence of PNG files in a specified directory.

    Parameters:
        directory_path (str): The path to the directory containing the PNG files.
        files_name (str): A wildcard pattern to match the PNG files in the directory.
        output_path (str): The output filename and path for the GIF.
        fps (int): Frames per second (speed) of the GIF animation.
        delete_pngs (bool): Whether to delete the PNG files after creating the GIF.

    Returns:
        None: This function does not return a value.

    Example:
        make_gif(
            directory_path='./images/',
            files_name='frame*.png',
            output_path='animation.gif',
            fps=20,
            delete_pngs=True
        )
    """
    # Use glob to get a list of all the PNG files in the directory
    png_files = sorted(glob.glob(os.path.join(directory_path, files_name)))

    # Use imageio to create the GIF
    with imageio.get_writer(output_path, mode='I', fps=fps) as writer:
        for file in png_files:
            with open(file, 'rb') as f:
                image = imageio.imread(f)
                writer.append_data(image)
            if delete_pngs:
                os.remove(file)

