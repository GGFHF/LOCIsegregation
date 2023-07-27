###############################################
#                                             #
#               0-Dependencies                #
#                                             #
###############################################
#This software has been developed by:

#  Dpto. Sistemas y Recursos Naturales
#  ETSI Montes, Forestal y del Medio Natural
#  Universidad Politecnica de Madrid
# https://github.com/ggfhf/
  
  
#This script sets the working directories and installs all the dependencies required to run the functions in the scripts
#to assess the seggregation distortion in the SNP from ScnI of López de Heredia et al. (2020)
#https://doi.org/10.3389/fpls.2020.564414

#We set several working directories for input and output data (Modify accordingly)
#WD.input.0 contains the scripts of the software package
WD.input.0 <- "/media/unaimantis/datos3/RCHEMA/Validacion221128/WD3/SCRIPTS-CHEMA-221216"
#WD.input.1 contains the allelic and genotypic frequencies input files, the mother genotyes, the LOC-ATG correspondence
#and the data extracted from Armendariz et al. (submitted) regarding RNAseq in Q. suber, Q. ilex and hybrids
WD.input <- "/media/unaimantis/datos3/RCHEMA/Validacion221128/WD3/SCRIPTS-CHEMA-221216/Input-files"
#WD.output.1 contains the output files with the results of the Chi-square test
WD.output.1 <- "/media/unaimantis/datos3/RCHEMA/Validacion221128/WD3/SCRIPTS-CHEMA-221216/Output-chisq-results"
#WD.output.2 contains the output files with the significant SNP and the Venn diagrams
WD.output.2 <- "/media/unaimantis/datos3/RCHEMA/Validacion221128/WD3/SCRIPTS-CHEMA-221216/Output-significant-venn"
#WD.output.3 contains the position of the significant SNP in the Q. robur chromosomes (linkage groups)
WD.output.3 <- "/media/unaimantis/datos3/RCHEMA/Validacion221128/WD3/SCRIPTS-CHEMA-221216/Output-robur-chromosomes"

setwd(WD.input.0)
#Dependencies installation
#Decomment the installation lines the first time using the software (some errors may arise depending on
#rstudio version and )
#install.packages("devtools")
library(devtools)
#devtools::install_github("bio-services/LinkageMapView", build_vignettes=TRUE)
library(LinkageMapView)
#install.packages("dplyr", dependencies = TRUE)
library(dplyr)
#install.packages("tidyverse", dependecies = TRUE)
library(tidyverse)
library(ggplot2) #already included in tidyverse, may be superfluous
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)







