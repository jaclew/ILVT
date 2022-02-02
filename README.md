# <img src="https://github.com/jaclew/IVLT/blob/main/ILVT.png" alt="Interactive Long-read Visualization Tool" width="500" height="74" align="middle">

## Table of Contents
* [Introduction](https://github.com/jaclew/IVLT/#introduction)
* [Installation](https://github.com/jaclew/IVLT/#installation)
  * [Linux](https://github.com/jaclew/IVLT/#installing-on-linux)
    * [Python version and dependencies](https://github.com/jaclew/IVLT/#python-version-and-dependencies)
    * [Installing dependencies via pip](https://github.com/jaclew/IVLT/#installing-dependencies-via-pip)
    * [Installing ILVT](https://github.com/jaclew/IVLT/#installing-ILVT)
  * [Windows](https://github.com/jaclew/IVLT/#installing-on-windows)
* [Launching](https://github.com/jaclew/IVLT/#launching)
* [Manual](https://github.com/jaclew/IVLT/#manual)
* [Citation](https://github.com/jaclew/IVLT/#citation)
* [Licence](https://github.com/jaclew/IVLT/#licence)


## Introduction
Interactive Long-Read Visualization Tool

## Installation
ILVT runs on Linux natively and can be run under Windows Subsystem for Linux (WSL).

### Installing on Linux
Follow below instructions to install on Linux.

#### Python version and dependencies
ILVT runs on Python version 3 (installed by default in Ubuntu 20.04) and depends on the following libraries (available on pip, package installer for Python)
* Pysam 
* PyQt5 
* Matplotlib 

#### Installing dependencies via pip
pip can be installed under Ubuntu by executing the following commands:
  > sudo apt update 
  > 
  > sudo apt install python3-pip 

The dependencies can be installed by executing the following commands:
  > pip install pysam pyqt5 matplotlib

#### Installing ILVT
1. Clone ILVT to your target location
     > git clone https://github.com/jaclew/ILVT
2. Launch ILVT using Python 3
   > python <path_to_ILVT>/read_painter.py

### Installing on Windows
Follow below instructions to install on Windows.
1. Enable Windows Subsystem for Linux (WSL) and install a Linux-distribution (e.g, Ubuntu)
   * Go to Windows Features
     * Turn on “Windows Subsystem for Linux” 
   * Reboot 
   * Go to Windows store
     * Install Ubuntu (version 20.04)
   * Launch WSL (Ubuntu) from the start-menu 
2. Installing and launching X-server to launch GUI-applications in WSL
   * Download and install <a href="https://sourceforge.net/projects/vcxsrv/" target="_blank">VcXsrv</a> (sourceforge).
   * To enable the display in WSL, type the following into the WSL terminal (Note: needs to be done on every new start):
       > export DISPLAY=localhost:0.0
     * The above line can be appended to file ~/.bashrc to avoid exporting the display with the above command. Restart the WSL to apply the command.
   * To launch the X-server, navigate to installation folder (Default path C:/Program Files/VcXsrv):
     * Launch a Windows terminal (Shift+F10 -> Open Terminal/PowerShell Here) and type the following to start VcXsrv
       > vcxsrv.exe -multiwindow -clipboard -wgl
3. Now follow the install instructions for Linux above.

## Launching
Launch ILVT with Python version 3 by executing the following command:
> python <path_to_ILVT>/read_painter.py

If you run on Windows and using WSL, make sure you have started the X-server (for instructions above).

## Manual
Documentation and examples (using a test dataset) are available in the reference manual, located in the file list above and at the following links: (<a href="https://github.com/jaclew/IVLT/raw/main/reference_manual.docx" target="_blank">reference manual</a>, <a href="https://github.com/jaclew/IVLT/tree/main/test_data" target="_blank">test dataset</a>).

## Citation
If you use ILVT in your research, please site the following publication:
> Auth1, Auth2, Auth3, Auth4.
> *Journal*, **Volume**:<pages>. [doi:<DOI_link>][doi]

## Licence
GNU General Public License, version 3
