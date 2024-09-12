![Copy-of-MosAIC-3-1200x737](https://github.com/user-attachments/assets/ce5ea9d1-58da-47bb-b9f5-e1bfe06c07fc)

## Overview 
**MosAIC: Mosquito-Associated Isolate Collection**
MosAIC (Mosquito-Associated Isolate Collection), a resource consisting of 392 bacterial isolates from mosquitoes. These isolates come with extensive metadata and high-quality draft genome assemblies, publicly available for use by the scientific community. Please see https://kcoonlab.bact.wisc.edu/mosaic/ for more information. 

## Analysis Overview 
Scripts for each analysis are written in R. Each directory contains necessary files and code to recreate each figure of the manuscript. To repeat the analysis, clone the repository, and the run each script. Do not `cd` into the cloned repository. 

**Example**
Create a new project in RStudio. To run the scripts in 01_GenomeQC to recreate Figure 1B in the manuscript. 
* Within terminal `git clone https://github.com/MosAIC-Collection/MosAIC_V1`
* Open `checkM_analysis.R` in the files panel window.
* `cmd enter` to run each step

## Recreate the Manuscript Figures
Once the repository has been cloned (above), recreate each figure as follows: 

**Fig 1 Origin of bacterial isolates in MosAIC**
* Script: `04_Sankey_Diagram/Metadata_Snakey_Diagram.R`: run code from line `1` to `127`

## Authors 
* Aidan Foo - aidanfoo96@gmail.com
* Laura Brettell - L.E.Brettell1@salford.ac.uk
* Eva Heinz - eva.heinz@strath.ac.uk
* Kerri Coon - kerri.coon@wisc.edu

## Citation 
Foo A, Brettell LE, Nichols HL, 2022 UW-Madison Capstone in Microbiology Students, Medina Muñoz M, Lysne JA, et al. Establishment and comparative genomics of a high-quality collection of mosquito-associated bacterial isolates - MosAIC (Mosquito-Associated Isolate Collection). 2023 Oct. Available from: http://biorxiv.org/lookup/doi/10.1101/2023.10.04.560816

