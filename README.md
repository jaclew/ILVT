# <img src="https://github.com/jaclew/ILVT/blob/main/logo_full.png" alt="Interactive Long-read Visualization Tool" width="387" height="195" align="middle">

## Table of Contents
* [Introduction](https://github.com/jaclew/ILVT/#introduction)
  * [Features](https://github.com/jaclew/ILVT/#features)
  * [Use cases](https://github.com/jaclew/ILVT/#use-cases)
* [Quick-start](https://github.com/jaclew/ILVT/#quick-start)
  * [Setup](https://github.com/jaclew/ILVT/#setup)
  * [Usage](https://github.com/jaclew/ILVT/#usage)
* [Installation](https://github.com/jaclew/ILVT/#installation)
  * [Linux](https://github.com/jaclew/ILVT/#installing-on-linux)
    * [Python version and dependencies](https://github.com/jaclew/ILVT/#python-version-and-dependencies)
    * [Installing dependencies via pip](https://github.com/jaclew/ILVT/#installing-dependencies-via-pip)
    * [Installing ILVT](https://github.com/jaclew/ILVT/#installing-ILVT)
  * [Windows](https://github.com/jaclew/ILVT/#installing-on-windows)
* [Manual](https://github.com/jaclew/ILVT/#manual)
* [Citation](https://github.com/jaclew/ILVT/#citation)
* [Licence](https://github.com/jaclew/ILVT/#licence)


## Introduction
Interactive Long-read Visualization Tool (ILVT) is a software with graphical user interface that enables a user (without bioinformatic expertise) to explore genomic rearrangements from long-read alignments. Rearrangements can be discovered by the user through interaction with individual reads (selected based on various criteria; such as alignment length and overlap to repetitive sequence, or in a specific region or gene) or on a genome-wide scale. It takes a couple of minutes to load the alignments of a Nanopore flowcell of cancer data on a modern laptop, and once loaded, the user can interact with reads within seconds.

ILVT minimally requires an alignment file and additional annotation files (gene GTF/GFF- and repetitive RepeatMasker-file) and variant call file (VCF) can be loaded.

### Features
* Interactive read selection and drawing
* Read lookup in specific regions (or genes)
* Identifying rearrangements genome-wide
* All drawings can be exported (e.g., in vector format)
* Convenient export of read names and (reference and query) coordinates

### Use cases
* Exploration of a new dataset to guide subsequent analysis design
* Rapid investigation of specific regions (e.g., gene fusions)
* Visualization of complex rearrangements (e.g., nested rearrangements in a cancer cell line)

## Quick-start
### Setup
ILVT is written in Python version 3 and depends on Python libraries Pysam, PyQt5 and Matplotlib. Download the software files to your location and launch the software with Python:
> git clone https://github.com/jaclew/ILVT
> 
> python <path_to_ILVT>/ilvt.py

If you run on Windows, see instructions below on how to install Windows Subsystem for Linux (WSL) and X-server (to launch graphical user interfaces under WSL).

### Usage
The minimal requirement to use the software is an alignment file. To import the alignment file, press **File** in the menu bar then **Import BAM**. The software will now load alignments of reads. This action may take some time depending on the size of the imported alignment file but enable fast display of reads once loaded. Optional import of annotation (RepeatMasker and GTF/GFF) and variant-call files (VCF) are loaded similarly.

Reads are first discovered (**Find reads** button in the right-side panel) according to some selection criteria (**Read selection**) to populate the Read-selection list where a read is drawn upon a left click. Right clicking on reads will add these to the current plot and this operation is useful to distinguish a putative rearrangement as heterogeneous or an artefact. (Right-clicked reads are only drawn within the boundaries of the current left-clicked reads and is increased or decreased under **View -> Show settings**).

Drawn reads appear in the plotting area. The user can zoom in by plus (+) or negative (-) keyboard buttons. Shift and alt modifiers control individual X and Y axis, respectively. Navigation is done by dragging the cursor inside the plot area. The view is resettable via menu bar **View -> Reset Zoom**.

To control what is drawn into the plotting area, settings can be modified under the **Show** section in the right panel. To update the current plot, click the **Update** button. Genes and variants can optionally be displayed if these files have been loaded into the software.

The Genome overview module is launched via the menu bar, **Genome Overview -> Launch module** and opens a secondary window. This module is useful to visualize all large-scale rearrangements in the genome and find the reads involved. The module works by binning the genome and querying alignments that are adjacent on the read but which overlap two non-adjacent bins.

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
   > python <path_to_ILVT>/ilvt.py

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

## Manual
Documentation and examples (using a test dataset) are available in the reference manual, located in the file list above and at the following links: (<a href="https://github.com/jaclew/ILVT/raw/main/reference_manual.docx" target="_blank">reference manual</a>, <a href="https://github.com/jaclew/ILVT/tree/main/test_data" target="_blank">test dataset</a>).

## Citation
If you use ILVT in your research, please site the following publication:
> Auth1, Auth2, Auth3, Auth4.
> *Journal*, **Volume**:<pages>. [doi:<DOI_link>][doi]

## Licence
GNU General Public License, version 3
