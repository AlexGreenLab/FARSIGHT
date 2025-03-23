# FARSIGHT
FARSIGHTs are aptamer-based sensors that enable rapid detection of nucleic acid point mutations. For more information on these systems, please read the preprint ["Programmable Fluorescent Aptamer-Based RNA Switches for Rapid Identification of Point Mutations"](https://www.medrxiv.org/content/10.1101/2025.03.07.25323576v1) by Zhaoqing Yan et al.

# Requirements
Like the [SNIPR design code](https://github.com/Albert09111/SNIPR/tree/master), the FARSIGHT design code uses Matlab and the multi-objective nucleic acid design and energy calculation functions of NUPACK. Users must download the [NUPACK source code](http://www.nupack.org/downloads) and compile it on their computers. Please follow the NUPACK manual to install NUPACK. 

After installation, it is necessary to add NUPACK to the system’s environmental variables to ensure that it can be located by the FARSIGHT code in Matlab. 

From a terminal window, open or create a .bash_profile file in your home directory with the command:

`open ~/.bash_profile`

For macOS, the equivalent file is named .zshrc. This file will be run whenever a bash shell is opened and can be used to set environment variables. In the resulting text editor environment, add the NUPACK directory as in the following example:

`export PATH=$PATH":/usr/local/bin"`

[For nano, save your changes by pressing ctrl + o and hit return to save. Exit nano by pressing ctrl + x.]

After setting the environment variables, **Matlab will need to be opened from a terminal window to use NUPACK-based functions.** A sample command to open Matlab is below and will vary depending on the Matlab install directory:

`/Applications/MATLAB_R2024b.app/bin/matlab`

The FARSIGHT design code also has several Matlab functions that are needed to call NUPACK, read NUPACK results, and manipulate nucleic acid sequences. These functions are included in the lib folder that should be in your FARSIGHT install folder. This lib folder is added to the Matlab path when you run the FARSIGHT design code.

This code has only been tested on an M-series Mac. Installation details may be different on other computers.

# Running FARSIGHT Design Code
## 1. Specify target sequences and reporter aptamers
Use the file mutant_target_input.csv in the install folder to provide information on the sequences that you want to detect. Provide a name for your target sequence and the wild-type (WT) and mutant (SNP) target sequences. The FARSIGHT code will design sensors for detection of the SNP target.

Use the file input_aptamer_set.m to specify the reporter aptamers that you want to use for the FARSIGHT designs. This file provides names, core sequences, and structures of fluorogenic aptamers like Broccoli, Corn, and Mango in standard and sometimes in rotated configurations. Comment out the appropriate lines of code in input_aptamer_set.m to indicate the specific reporter aptamers to use for the FARSIGHT designs.

## 2. Define the FARSIGHT designs to generate
The next step is to define a set of different FARSIGHT designs with varying domain lengths and binding sites along the target transcript. For each target/reporter aptamer combination, the code will define ~72 FARSIGHT design combinations. To do this, run the following in Matlab:

`FARSIGHT_stage1_define`

## 3. Generate the FARSIGHT designs
This stage is the most time-consuming part of the FARSIGHT design process. The code will take in the set of design definitions compiled in stage 1 and use NUPACK multi-objective design to generate RNA sequences that satisfy the specified secondary structure and sequence constraints. It will also calculate the free energies and defect levels of the FARSIGHT sequences generated. The code will generate `num_designs_overall` designs for each FARSIGHT domain length combination and aptamer/target combination. By default, `num_designs_overall = 2` but it can be changed to reduce design time or increase the range of designs. To begin this process, run the following in Matlab:

`FARSIGHT_stage2_generate`

To reduce the time needed for this stage, you can call this script using a command-line terminal and run it using multiple Matlab instantiations. To do this in macOS, open a Terminal window and navigate to the install folder. Then run the command below:

`/Applications/MATLAB_R2024b.app/bin/matlab < FARSIGHT_stage2_generate.m`

Replace the `/Applications/MATLAB_R2024b.app/bin/matlab` with the appropriate command to start Matlab depending on its install location. You can run this command in multiple terminal windows depending on the number of processor cores and the amount of available system memory. It is also possible to use parallel processing in FARSIGHT_stage2_generate.m to speed up the design process by changing the loops on Lines 74 and 190 from `for` to `parfor`. However, this parallelization approach can lead to Matlab stalling and provide slower performance compared to the Terminal approach with multiple Matlab instances. FARSIGHT_stage2_generate is designed so that it can be terminated at any time and resume operation where it left off.

## 4. Score and select the FARSIGHT designs
In this stage, the FARSIGHT designs are compiled, assigned scores, and selected based on those scores. The FARSIGHT designs are assigned a score based on their ∆∆G values and defect levels. The top `num_designs = 8` for each aptamer/target combination are output to FARSIGHT_design_output/FARSIGHT_top_designs.csv. You can change the value of `num_designs` to select more or fewer FARSIGHTs for testing. The full set of FARSIGHTs generated by the code is also compiled based on the aptamer/target combination in files with the naming format FARSIGHT_design_output/*[aptamer name]*_*[target name]*_full_design_info.csv. These files also contain a set of free energy calculations and defect levels computed for each design for straightforward implementation of alternative scoring functions. To complete this stage, run the following in Matlab:

`/Applications/MATLAB_R2024b.app/bin/matlab < FARSIGHT_stage3_select.m`
