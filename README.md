# DNApaintDeNoise

This repository contains example analysis scripts written in MATLAB for a manuscript currently under review.

### 1. System Requirements: 
MATLAB R2023b (64-bit) or newer.

### 2. Installation Guide: 
Download or clone this repository to a single folder on your computer. Download the example data "ExampleData.mat" and put it in the same local directory. 

### 3. Demo:
Open MATLAB and navigate to the folder with the MATLAB scripts and the example data.<br/>
In the Command Window type 'run("ExampleCode.m")'.<br/>
When the script finishes running, a figure will appear with two scatter plots showing the raw and denoised point cloud data.<br/>
Typical run time for a demo on a normal computer is ~3 min.

### 5. Instructions for use:
The script will run on any data with a four-column format as follows:

x-position (pixel), y-position (pixel), time (frame number), localization precision (pixel)
