# -*- coding: utf-8 -*-

"""
MAKE_MOVIE

Author: Corentin K. Louis
Date: 2023-10-26
Description: This module contains a routine that allows to create a movie from a sequence of PNG files in a specified directory
"""

import os
import glob
import subprocess


def make_movie(
    directory_path: str = './',
    files_name: str = 'filename*.png',
    output_path: str = 'output.mp4',
    fps: int = 30,
    delete_pngs: bool = False
) -> None:
    png_files = sorted(glob.glob(os.path.join(directory_path, files_name)))

    # Construct the FFmpeg command to create the movie
    command = [
        "ffmpeg",
        "-i", os.path.join(directory_path, files_name),
        "-y",
        "-acodec aac",
        "-ac 2",
        "-ab 160k",
        "-vcodec libx264",
        "-vf 'scale='bitand(oh*dar,65534)':'min(720,ih)'",
        "-f mp4",
        "-framerate", str(fps),
        #"-pattern_type", "glob",
        #"-c:v", "libx264",
        #"-pix_fmt", "yuv420p",
        output_path
    ]

    # Run the FFmpeg command
    subprocess.run(command)

    if delete_pngs:
        for file in png_files:
            os.remove(file)

