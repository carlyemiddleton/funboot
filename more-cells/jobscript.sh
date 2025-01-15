#!/bin/bash
module unload R/4.4.0
module unload r-spicyr/1.16.4-gcc-8.3.0-r-4.3.2

module load R/4.4.0
Rscript setup.R > output_setup.txt
module unload R/4.4.0

module load r-spicyr/1.16.4-gcc-8.3.0-r-4.3.2 
Rscript spicy.R > output_spicy.txt
module unload r-spicyr/1.16.4-gcc-8.3.0-r-4.3.2

