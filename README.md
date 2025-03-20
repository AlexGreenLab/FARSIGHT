# FARSIGHT
FARSIGHTs are aptamer-based sensors that enable rapid detection of nucleic acid point mutations. For more information on these systems, please read the preprint ["Programmable Fluorescent Aptamer-Based RNA Switches for Rapid Identification of Point Mutations" by Z. Yan et al.](https://www.medrxiv.org/content/10.1101/2025.03.07.25323576v1).

# Requirements
Like the [SNIPR design code](https://github.com/Albert09111/SNIPR/tree/master), the FARSIGHT design code uses Matlab and the multi-objective nucleic acid design and energy calculation functions of NUPACK. Users must download the [NUPACK source code](http://www.nupack.org/downloads) and compile it on their computers. Please follow the NUPACK manual to install NUPACK. 

After installation, it is necessary to add NUPACK to the system’s environmental variables to ensure that it can be located by the FARSIGHT code. 

From a terminal window, open or create a .bash_profile file in your home directory with the command:

`open ~/.bash_profile`

This file will be run whenever a bash shell is opened and can be used to set environment variables. In the resulting text editor environment, add the NUPACK directory as the following example.

`export NUPACKHOME=/home/Documents/nupack3.2.2`

`export PATH=$PATH":${NUPACKHOME}/bin"`

[For nano, save your changes by typing ctrl + o and hit return to save. Exit nano by typing ctrl + x.]

After setting the environment variables, Matlab will need to be opened from a terminal window to use NUPACK based functions. A sample command to open Matlab is below and will vary depending on the Matlab install directory:

`/Applications/MATLAB_R2024b.app/bin/matlab`
