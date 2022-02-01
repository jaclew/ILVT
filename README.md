# <img src="https://github.com/jaclew/IVLT/blob/main/ILVT.png" alt="Interactive Long-read Visualization Tool" width="500" height="74" align="middle">

## Table of Contents
* [Introduction](https://github.com/jaclew/IVLT/#introduction)
* [Installation](https://github.com/jaclew/IVLT/#installation)
  * [Linux](https://github.com/jaclew/IVLT/#installing-on-linux)
  * [Windows](https://github.com/jaclew/IVLT/#installing-on-windows)
* [Manual](https://github.com/jaclew/IVLT/#manual)
* [Citation](https://github.com/jaclew/IVLT/#citation)
* [Licence](https://github.com/jaclew/IVLT/#licence)


## Introduction
Interactive Long-Read Visualization Tool

## Installation
ILVT runs on Linux natively and can be run under Windows Subsystem for Linux (WSL).

### Installing on Linux
Follow below instructions to install on Linux.

### Installing on Windows
Follow below instructions to install on Windows.
1. Enable Windows Subsystem for Linux (WSL) and install a Linux-distribution (e.g, Ubuntu)
   * Go to Windows Features
     * Turn on “Windows Subsystem for Linux” 
   * Reboot 
   * Go to Windows store
     * Install Ubuntu (version 20.04)
   * Launch WSL from start-menu 
2. Installing and launching X-server to launch GUI-applications in WSL
   * Download and install <a href="https://sourceforge.net/projects/vcxsrv/" target="_blank">VcXsrv</a> (sourceforge).
   * To enable the display in WSL, type the following into the terminal
     * > Type in terminal: “export DISPLAY=localhost:0.0” (or append to file ~/.bashrc and restart WSL) 
   * To launch the X-server, navigate to installation folder (Default path C:/Program Files/VcXsrv)
     * Launch windows terminal (Shift+F10 -> Open Terminal/PowerShell Here)
     * Type “vcxsrv.exe -multiwindow -clipboard -wgl” into windows terminal to start VcXsrv 

## Manual
Documentation and example use cases are available in the <a href="https://github.com/jaclew/IVLT/raw/main/reference_manual.docx" target="_blank">reference manual</a>.

## Citation
If you use ILVT in your research, please site the following publication:
> Auth1, Auth2, Auth3, Auth4.
> *Journal*, **Volume**:<pages>. [doi:<DOI_link>][doi]

## Licence
GNU General Public License, version 3
