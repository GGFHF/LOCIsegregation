# LOCIsegregation
R scripts to analyze the distortion of the mendelian segregation patterns in Q. ilex x suber hybrid progenies
####################################################################################################
#                                                                                                  #
#   Scripts to assess and quantify seggregation distortion on SNP from sclerophylous oaks hybrids  #
#                                                                                                  #
####################################################################################################


#This software has been developed by:

	#  Dpto. Sistemas y Recursos Naturales
	#  ETSI Montes, Forestal y del Medio Natural
	#  Universidad Politecnica de Madrid
	# https://github.com/ggfhf/


This software was conceived as a series of scripts aiming to identify, quantify and extract information from the SNP obtained in ScnI (full imputation scenario) for Q. ilex, Q. suber and their hybrids (López de Heredia et al. 2020). The software is part of the Bachelor thesis “ANÁLISIS DE LA SEGREGACIÓN ALÉLICA EN PROGENIES DE Quercus ilex L., Q. suber L. Y SUS HÍBRIDOS” authored by José María Conde González and supervised by prof. Álvaro Soto and Unai López de Heredia

Input files (stored in a single directory, WD.input):
	ScnI-FrecAlelicas-Bialelicos.txt
	ScnI-FrecGenotBialelicos-AL.txt
	ScnI-FrecGenotBialelicos-EFS.txt
	ScnI-FrecGenotBialelicos-EN.txt
	ScnI-MotherGenotypes-biallelic.txt
	A05-FrecGenot-biallelic-sorted3.txt
	A07-FrecGenot-biallelic-sorted3.txt
	A09-FrecGenot-biallelic-sorted3.txt
	A10-FrecGenot-biallelic-sorted3.txt
	E28-FrecGenot-biallelic-sorted3.txt
	E31-FrecGenot-biallelic-sorted3.txt
	E41-FrecGenot-biallelic-sorted3.txt
	E96-FrecGenot-biallelic-sorted3.txt
	FS16_S-FrecGenot-biallelic-sorted3.txt
	FS19_S-FrecGenot-biallelic-sorted3.txt
	FS20_S-FrecGenot-biallelic-sorted3.txt
	FS22_S-FrecGenot-biallelic-sorted3.txt
	CorrespondenceLOC-ATG.txt
	RNASEQ-DEGSa.txt
	input_LMAP_viewer4.txt

R scripts (stored in WD.input.0:
	0-Config.R
	1-Chi2-test-adults.R
	2-Chi2-test-suber-ilex-progenies.R
	3-Chi2-test-hybrid-progenies.R
	4-Significant-SNP-Venn-Diagrams.R
	5-Crossing-RNAseq-data.R
	6-Qrobur-linkage-groups-chromosomes.R



DESCRIPTION OF THE SCRIPTS
0-Config.R
This script contains the information to install and load the following R libraries:
dplyr, ggvenn, tidyverse, LinkageMapView
It also sets the working directories. IMPORTANT: Run this script before running the others.

1-Chi2-test-adults.R
This script computes the Chi-square test for all the biallelic imputed and non-imputed SNP found in the adult hybrids (EFS), Q. Ilex (EN) and Q. Suber (AL) according to the expected genotypic frequencies computed from the allelic frequencies and to the observed genotypic frequencies (0/0, 0/1, 1/1, 0/99, 99/99). Results for all SNP are stored in the directory prompted by the WD.output.1 defined by the user.

2-Chi2-test-suber-ilex-progenies.R
This script computes the Chi-square test for all the biallelic imputed and non-imputed SNP found in the pures species progenies (A05, A07, A09, A10, E28, E31, E41, E96) according to the expected genotypic frequencies computed from the allelic frequencies and to the observed genotypic frequencies (0/0, 0/1, 1/1, 0/99, 99/99). Results for all SNP are stored in the directory prompted by the WD.output.1 defined by the user.

3-Chi2-test-hybrid-progenies.R
This script computes the Chi-square test for all the biallelic imputed and non-imputed SNP found in the hybrid progenies (FS16, FS19, FS20, FS22) according to the expected genotypic frequencies computed from the allelic frequencies and to the observed genotypic frequencies (0/0, 0/1, 1/1, 0/99, 99/99). Results for all SNP are stored in the directory prompted by the WD.output.1 defined by the user.

4-Significant-SNP-Venn-Diagrams.R
This script subsets the files produced with scripts 1, 2 and 3 keeping the SNP with a p-value < 0.05 at the Chi-square test, and generates a stats file with the number of significant SNP. It also computes Venn diagrams for SNP and genes of the different families and adult sets. Results for all SNP are stored in the directory prompted by the WD.output.2 defined by the user.

5-Crossing-RNAseq-data.R
This script crosses the significant SNP shared by pure and hybrid progenies and adults with the information produced in the RNAseq study (Armendariz et al. Submitted). It also identifies the orthologs in Arabidopsis to further run the ShinyGO gene enrichment analysis. Results are stored in the directory prompted by the WD.output.2 defined by the user.

6-Qrobur-linkage-groups-chromosomes.R
This scripr produces graphics showing the position in the linkage groups/chromosomes of Q. robur of the significant SNP share by hybrid adults and by hybrid progenies. Results are stored in the directory prompted by the WD.output.3 defined by the user..

HOW TO RUN THE SCRIPTS

Be sure to run first the config script (0-Dependencies.R). After doing this, it is only necessary to run the last script (6-Qrobur-linkage-groups-chromosomes.R) that has the information to call the rest of scripts.
