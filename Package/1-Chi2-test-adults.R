##########################################################
#                                                        #
#        2-Chi2-test for suber-ilex-hybrid adults        #
#                                                        #
##########################################################
#This software has been developed by:

#  Dpto. Sistemas y Recursos Naturales
#  ETSI Montes, Forestal y del Medio Natural
#  Universidad Politecnica de Madrid
# https://github.com/ggfhf/

#We set the working directory to run the scripts. Be sure all the scripts are within this directory
#setwd(WD.input.0)

#We run the config file
#source('0-Config.R')

#We set the working directory. Be sure all the input files are within this directory
setwd(WD.input)

#The required input files here are:
#ScnI-FrecAlelicas-Bialelicos.txt
#ScnI-FrecGenotBialelicos-EFS.txt
#ScnI-FrecGenotBialelicos-EN.txt
#ScnI-FrecGenotBialelicos-AL.txt


#We scan the above files to have data.frames for allele frequencies (common)
# and for observed genotypic frequencies obtained with NGShelper
#(one data.frame for each species, EFS: hybrids, EN: ilex, AL: suber)
#Here we only consider bi-allelic imputed and non-imputed loci

FrecuenciasAlelicas <- read.table("ScnI-FrecAlelicas-Bialelicos.txt", 
                                  header= TRUE, 
                                  sep = "\t", 
                                  dec = ".")
FrecuenciasGenotipicas.EFS <- read.table("ScnI-FrecGenotBialelicos-EFS.txt", 
                                         header = TRUE, 
                                         sep = "\t", 
                                         dec = ".")
FrecuenciasGenotipicas.EN <- read.table("ScnI-FrecGenotBialelicos-EN.txt", 
                                        header = TRUE, 
                                        sep = "\t", 
                                        dec = ".")
FrecuenciasGenotipicas.AL <- read.table("ScnI-FrecGenotBialelicos-AL.txt", 
                                        header = TRUE, 
                                        sep = "\t", 
                                        dec = ".")

#We set the working directory to produce the output.
setwd(WD.output.1)

############################
#                          #
#     Hybrid adults:EFS    #
#                          #  
############################


#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#using function from dplyr

FrecGeno2.EFS <- FrecuenciasGenotipicas.EFS %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively

bialelic.EFS <- data.frame(FrecuenciasAlelicas$variant_id,
                   FrecuenciasAlelicas$seq_id,
                   FrecuenciasAlelicas$position,
                   FrecuenciasAlelicas$imputations,
                   FrecuenciasAlelicas$allele_id,
                   FrecuenciasAlelicas$AL_frequency,
                   FrecuenciasAlelicas$EN_frequency,
                   FrecuenciasAlelicas$HY_frequency,
                   FrecuenciasAlelicas$AL_mothers_frequency,
                   FrecuenciasAlelicas$EN_mothers_frequency,
                   FrecuenciasAlelicas$HY_mothers_frequency,
                   FrecuenciasAlelicas$AL_progenies_frequency,
                   FrecuenciasAlelicas$EN_progenies_frequency,
                   FrecuenciasAlelicas$HY_progenies_frequency,
                   FrecuenciasAlelicas$genomic_zone,
                   FrecuenciasAlelicas$gene.fragment,
                   FrecuenciasAlelicas$description,
                   FrecuenciasAlelicas$chromosome_id,
                   FrecuenciasAlelicas$bases,
                   FrecGeno2.EFS$X0_0,
                   FrecGeno2.EFS$X0_1,
                   FrecGeno2.EFS$X1_1,
                   FrecGeno2.EFS$X0_99,
                   FrecGeno2.EFS$X1_99,
                   FrecGeno2.EFS$X99_99)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.EFS <- subset(bialelic.EFS, bialelic.EFS$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.EFS <- subset(bialelic.EFS, bialelic.EFS$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2

i <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
FrecuenciasNOimputados.EFS <- list()

while (i <= nrow(bialelicosNOimput.EFS)) {
  
  
  esp1 <- bialelicosNOimput.EFS[i,6] * bialelicosNOimput.EFS[i,7]
  esp2 <- bialelicosNOimput.EFS[i,6] * bialelicosNOimput.EFS[i+1,7] + bialelicosNOimput.EFS[i+1,6]*bialelicosNOimput.EFS[i,7]
  esp3 <- bialelicosNOimput.EFS[i+1,6] * bialelicosNOimput.EFS[i+1,7]
  
  obs1 <- as.double(paste(bialelicosNOimput.EFS[i,20]))
  obs2 <- as.double(paste(bialelicosNOimput.EFS[i,21]))
  obs3 <- as.double(paste(bialelicosNOimput.EFS[i,22]))
  
  FrecuenciasNOimputados.EFS[[a]] <- data.frame(x = c(obs1*22,obs2*22,obs3*22), 
                                                y = c(esp1,esp2,esp3))
  
  
  i <- i+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(22*obs1,22*obs2,22*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.EFS <- 1

while (g < nrow(bialelicosNOimput.EFS)) {
  
  m.EFS <- FrecuenciasNOimputados.EFS[[l]]
  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOimputados.EFS[[l]][,2][1] == 0){
    FrecuenciasNOimputados.EFS[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.EFS <- FrecuenciasNOimputados.EFS[[l]]
  }
  if(FrecuenciasNOimputados.EFS[[l]][,2][2] == 0){
    FrecuenciasNOimputados.EFS[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.EFS <- FrecuenciasNOimputados.EFS[[l]]
  }
  if(FrecuenciasNOimputados.EFS[[l]][,2][3] == 0){
    FrecuenciasNOimputados.EFS[[l]][,2][3] <- 0.000000000000000000000000000000000000001
     m.EFS <- FrecuenciasNOimputados.EFS[[l]]
   }
  
 v.EFS <- chisq.test(x = m.EFS[,1],p = m.EFS[,2], rescale.p = TRUE)

 #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 26 and 27  
  bialelicosNOimput.EFS[g,26] <- v.EFS$statistic
  bialelicosNOimput.EFS[g,27] <- v.EFS$p.value
  bialelicosNOimput.EFS[g+1,26] <- v.EFS$statistic
  bialelicosNOimput.EFS[g+1,27] <- v.EFS$p.value
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.EFS , 
            "BialelicosNOImputados-EFS.txt", 
            sep="\t",
            dec=",")

###################
#                 #
#     Imputed     #
#                 #   
###################

#We repeat the calculations for imputed SNP

#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/99 y 99/99 in column 1 
#and with the expected frequencies for 0/0, 0/99 y 99/99 in column 2
FrecuenciasImputados.EFS <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0

while (h <= nrow(bialelicosimput.EFS)) {
  
  #The uploaded allele frequency file presents consecutive lines for each SNP
  #in an unordered way: sometimes 0 allele comes first, and some others allele 99
  #We consider this issue to perform the calculations

  if(bialelicosimput.EFS[h,5] == 0 && bialelicosimput.EFS[h+1,5] == 99 ){ 
    
    esp1 <- bialelicosimput.EFS[h,6] * bialelicosimput.EFS[h,7]
    esp2 <- bialelicosimput.EFS[h,6] * bialelicosimput.EFS[h+1,7] + bialelicosimput.EFS[h+1,6]*bialelicosimput.EFS[h,7]
    esp3 <- bialelicosimput.EFS[h+1,6] * bialelicosimput.EFS[h+1,7]
    
    obs1 <- as.double(paste(bialelicosimput.EFS[h,20]))
    obs2 <- as.double(paste(bialelicosimput.EFS[h,23]))
    obs3 <- as.double(paste(bialelicosimput.EFS[h,25]))
  }
  if(bialelicosimput.EFS[h,5] == 99 && bialelicosimput.EFS[h+1,5] == 0 ){ 
    
    esp1 <- bialelicosimput.EFS[h+1,6] * bialelicosimput.EFS[h+1,7]
    esp2 <- bialelicosimput.EFS[h+1,6] * bialelicosimput.EFS[h,7] + bialelicosimput.EFS[h,6]*bialelicosimput.EFS[h+1,7]
    esp3 <- bialelicosimput.EFS[h,6] * bialelicosimput.EFS[h,7]
    
    obs1 <- as.double(paste(bialelicosimput.EFS[h+1,20]))
    obs2 <- as.double(paste(bialelicosimput.EFS[h+1,23]))
    obs3 <- as.double(paste(bialelicosimput.EFS[h+1,25]))
  } 
  
  FrecuenciasImputados.EFS[[a]] <- data.frame(x = c(obs1*22,obs2*22,obs3*22), 
                                              y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(22*obs1,22*obs2,22*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.EFS <- 1

while (g < nrow(bialelicosimput.EFS)) {
  
  m.EFS <- FrecuenciasImputados.EFS[[l]]
  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.EFS[[l]][,2][1] == 0){
    FrecuenciasImputados.EFS[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.EFS <- FrecuenciasImputados.EFS[[l]]
  }
  if(FrecuenciasImputados.EFS[[l]][,2][2] == 0){
    FrecuenciasImputados.EFS[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.EFS <- FrecuenciasImputados.EFS[[l]]
  }
  if(FrecuenciasImputados.EFS[[l]][,2][3] == 0){
    FrecuenciasImputados.EFS[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.EFS <- FrecuenciasImputados.EFS[[l]]
  }
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 26 and 27  
  
    v.EFS <- chisq.test(x = m.EFS[,1],p = m.EFS[,2], rescale.p = TRUE)
  
  bialelicosimput.EFS[g,26] <- v.EFS$statistic
  bialelicosimput.EFS[g,27] <- v.EFS$p.value
  bialelicosimput.EFS[g+1,26] <- v.EFS$statistic
  bialelicosimput.EFS[g+1,27] <- v.EFS$p.value
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.EFS , 
            "BialelicosImputados-EFS.txt", 
            sep="\t",
            dec=",")


############################
#                          #
#      Adults Ilex: EN     #
#                          #  
############################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#using function from dplyr

FrecGeno2.EN <- FrecuenciasGenotipicas.EN %>% slice(rep(1:n(),each =2))
dim(FrecGeno2.EN)

#We produce a new data.frame including columns from the allele frequencies file and
#from the observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively

bialelic.EN <- data.frame(FrecuenciasAlelicas$variant_id,
                          FrecuenciasAlelicas$seq_id,
                          FrecuenciasAlelicas$position,
                          FrecuenciasAlelicas$imputations,
                          FrecuenciasAlelicas$allele_id,
                          FrecuenciasAlelicas$AL_frequency,
                          FrecuenciasAlelicas$EN_frequency,
                          FrecuenciasAlelicas$HY_frequency,
                          FrecuenciasAlelicas$AL_mothers_frequency,
                          FrecuenciasAlelicas$EN_mothers_frequency,
                          FrecuenciasAlelicas$HY_mothers_frequency,
                          FrecuenciasAlelicas$AL_progenies_frequency,
                          FrecuenciasAlelicas$EN_progenies_frequency,
                          FrecuenciasAlelicas$HY_progenies_frequency,
                          FrecuenciasAlelicas$genomic_zone,
                          FrecuenciasAlelicas$gene.fragment,
                          FrecuenciasAlelicas$description,
                          FrecuenciasAlelicas$chromosome_id,
                          FrecuenciasAlelicas$bases,
                          FrecGeno2.EN$X0_0,
                          FrecGeno2.EN$X0_1,
                          FrecGeno2.EN$X1_1,
                          FrecGeno2.EN$X0_99,
                          FrecGeno2.EN$X1_99,
                          FrecGeno2.EN$X99_99)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.EN <- subset(bialelic.EN, bialelic.EN$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.EN <- subset(bialelic.EN, bialelic.EN$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2

i <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
FrecuenciasNOimputados.EN <- list()

while (i <= nrow(bialelicosNOimput.EN)) {
  
  
  esp1 <- bialelicosNOimput.EN[i,7] * bialelicosNOimput.EN[i,7]
  esp2 <- bialelicosNOimput.EN[i,7] * bialelicosNOimput.EN[i+1,7] + bialelicosNOimput.EN[i+1,7]*bialelicosNOimput.EN[i,7]
  esp3 <- bialelicosNOimput.EN[i+1,7] * bialelicosNOimput.EN[i+1,7]
  
  obs1 <- as.double(paste(bialelicosNOimput.EN[i,20]))
  obs2 <- as.double(paste(bialelicosNOimput.EN[i,21]))
  obs3 <- as.double(paste(bialelicosNOimput.EN[i,22]))
  
  FrecuenciasNOimputados.EN[[a]] <- data.frame(x = c(obs1*99,obs2*99,obs3*99), 
                                               y = c(esp1,esp2,esp3))
  
  
  i <- i+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(99*obs1,99*obs2,99*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.EN <- 1

while (g < nrow(bialelicosNOimput.EN)) {
  
  m.EN <- FrecuenciasNOimputados.EN[[l]]
  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOimputados.EN[[l]][,2][1] == 0){
    FrecuenciasNOimputados.EN[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.EN <- FrecuenciasNOimputados.EN[[l]]
  }
  if(FrecuenciasNOimputados.EN[[l]][,2][2] == 0){
    FrecuenciasNOimputados.EN[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.EN <- FrecuenciasNOimputados.EN[[l]]
  }
  if(FrecuenciasNOimputados.EN[[l]][,2][3] == 0){
    FrecuenciasNOimputados.EN[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.EN <- FrecuenciasNOimputados.EN[[l]]
  }
  
  v.EN <- chisq.test(x = m.EN[,1],p = m.EN[,2], rescale.p = TRUE)

  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 26 and 27   
  bialelicosNOimput.EN[g,26] <- v.EN$statistic
  bialelicosNOimput.EN[g,27] <- v.EN$p.value
  bialelicosNOimput.EN[g+1,26] <- v.EN$statistic
  bialelicosNOimput.EN[g+1,27] <- v.EN$p.value
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.EN , 
            "BialelicosNOImputados-EN.txt", 
            sep="\t",
            dec=",")

###################
#                 #
#     Imputed     #
#                 #   
###################

#We repeat the calculations for imputed SNP

#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/99 y 99/99 in column 1 
#and with the expected frequencies for 0/0, 0/99 y 99/99 in column 2

FrecuenciasImputados.EN <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0

while (h <= nrow(bialelicosimput.EN)) {
  
  #The uploaded allele frequency file presents consecutive lines for each SNP
  #in an unordered way: sometimes 0 allele comes first, and some others allele 99
  #We consider this issue to perform the calculations
  
  if(bialelicosimput.EN[h,5] == 0 && bialelicosimput.EN[h+1,5] == 99 ){   
  
  esp1 <- bialelicosimput.EN[h,7] * bialelicosimput.EN[h,7]
  esp2 <- bialelicosimput.EN[h,7] * bialelicosimput.EN[h+1,7] + bialelicosimput.EN[h+1,7]*bialelicosimput.EN[h,7]
  esp3 <- bialelicosimput.EN[h+1,7] * bialelicosimput.EN[h+1,7]
  
  obs1 <- as.double(paste(bialelicosimput.EN[h,20]))
  obs2 <- as.double(paste(bialelicosimput.EN[h,23]))
  obs3 <- as.double(paste(bialelicosimput.EN[h,25]))

  }  
  
  if(bialelicosimput.EN[h,5] == 99 && bialelicosimput.EN[h+1,5] == 00 ){   
    
    esp1 <- bialelicosimput.EN[h,7+1] * bialelicosimput.EN[h,7+1]
    esp2 <- bialelicosimput.EN[h,7+1] * bialelicosimput.EN[h,7] + bialelicosimput.EN[h,7]*bialelicosimput.EN[h+1,7]
    esp3 <- bialelicosimput.EN[h,7] * bialelicosimput.EN[h,7]
    
    obs1 <- as.double(paste(bialelicosimput.EN[h,20]))
    obs2 <- as.double(paste(bialelicosimput.EN[h,23]))
    obs3 <- as.double(paste(bialelicosimput.EN[h,25]))
    
  }     
  FrecuenciasImputados.EN[[a]] <- data.frame(x = c(obs1*99,obs2*99,obs3*99), 
                                             y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(22*obs1,22*obs2,22*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.EN <- 1

while (g < nrow(bialelicosimput.EN)) {
  
  m.EN <- FrecuenciasImputados.EN[[l]]
  
  #The uploaded allele frequency file presents consecutive lines for each SNP
  #in an unordered way: sometimes 0 allele comes first, and some others allele 99
  #We consider this issue to perform the calculations
  
  if(FrecuenciasImputados.EN[[l]][,2][1] == 0){
    FrecuenciasImputados.EN[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.EN <- FrecuenciasImputados.EN[[l]]
  }
  if(FrecuenciasImputados.EN[[l]][,2][2] == 0){
    FrecuenciasImputados.EN[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.EN <- FrecuenciasImputados.EN[[l]]
  }
  if(FrecuenciasImputados.EN[[l]][,2][3] == 0){
    FrecuenciasImputados.EN[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.EN <- FrecuenciasImputados.EN[[l]]
  }
  
  v.EN <- chisq.test(x = m.EN[,1],p = m.EN[,2], rescale.p = TRUE)

  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 26 and 27     
  bialelicosimput.EN[g,26] <- v.EN$statistic
  bialelicosimput.EN[g,27] <- v.EN$p.value
  bialelicosimput.EN[g+1,26] <- v.EN$statistic
  bialelicosimput.EN[g+1,27] <- v.EN$p.value
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.EN , 
            "BialelicosImputados-EN.txt", 
            sep="\t",
            dec=",")


############################
#                          #
#     Adults Suber: AL     #
#                          #  
############################
#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#using function from dplyr

FrecGeno2.AL <- FrecuenciasGenotipicas.AL %>% slice(rep(1:n(),each =2))
dim(FrecGeno2.AL)

#We produce a new data.frame including columns from the allele frequencies file and
#from the observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively

bialelic.AL <- data.frame(FrecuenciasAlelicas$variant_id,
                          FrecuenciasAlelicas$seq_id,
                          FrecuenciasAlelicas$position,
                          FrecuenciasAlelicas$imputations,
                          FrecuenciasAlelicas$allele_id,
                          FrecuenciasAlelicas$AL_frequency,
                          FrecuenciasAlelicas$EN_frequency,
                          FrecuenciasAlelicas$HY_frequency,
                          FrecuenciasAlelicas$AL_mothers_frequency,
                          FrecuenciasAlelicas$EN_mothers_frequency,
                          FrecuenciasAlelicas$HY_mothers_frequency,
                          FrecuenciasAlelicas$AL_progenies_frequency,
                          FrecuenciasAlelicas$EN_progenies_frequency,
                          FrecuenciasAlelicas$HY_progenies_frequency,
                          FrecuenciasAlelicas$genomic_zone,
                          FrecuenciasAlelicas$gene.fragment,
                          FrecuenciasAlelicas$description,
                          FrecuenciasAlelicas$chromosome_id,
                          FrecuenciasAlelicas$bases,
                          FrecGeno2.AL$X0_0,
                          FrecGeno2.AL$X0_1,
                          FrecGeno2.AL$X1_1,
                          FrecGeno2.AL$X0_99,
                          FrecGeno2.AL$X1_99,
                          FrecGeno2.AL$X99_99)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.AL <- subset(bialelic.AL, bialelic.AL$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.AL <- subset(bialelic.AL, bialelic.AL$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2

i <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
FrecuenciasNOimputados.AL <- list()

while (i <= nrow(bialelicosNOimput.AL)) {
  
  
  esp1 <- bialelicosNOimput.AL[i,6] * bialelicosNOimput.AL[i,6]
  esp2 <- bialelicosNOimput.AL[i,6] * bialelicosNOimput.AL[i+1,6] + bialelicosNOimput.AL[i+1,6]*bialelicosNOimput.AL[i,6]
  esp3 <- bialelicosNOimput.AL[i+1,6] * bialelicosNOimput.AL[i+1,6]
  
  obs1 <- as.double(paste(bialelicosNOimput.AL[i,20]))
  obs2 <- as.double(paste(bialelicosNOimput.AL[i,21]))
  obs3 <- as.double(paste(bialelicosNOimput.AL[i,22]))
  
  FrecuenciasNOimputados.AL[[a]] <- data.frame(x = c(obs1*98,obs2*98,obs3*98), 
                                               y = c(esp1,esp2,esp3))
  
  
  i <- i+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(98*obs1,98*obs2,98*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.AL <- 1

while (g < nrow(bialelicosNOimput.AL)) {
  
  m.AL <- FrecuenciasNOimputados.AL[[l]]
  
  #The uploaded allele frequency file presents consecutive lines for each SNP
  #in an unordered way: sometimes 0 allele comes first, and some others allele 99
  #We consider this issue to perform the calculations
  
  if(FrecuenciasNOimputados.AL[[l]][,2][1] == 0){
    FrecuenciasNOimputados.AL[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.AL <- FrecuenciasNOimputados.AL[[l]]
  }
  if(FrecuenciasNOimputados.AL[[l]][,2][2] == 0){
    FrecuenciasNOimputados.AL[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.AL <- FrecuenciasNOimputados.AL[[l]]
  }
  if(FrecuenciasNOimputados.AL[[l]][,2][3] == 0){
    FrecuenciasNOimputados.AL[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.AL <- FrecuenciasNOimputados.AL[[l]]
  }

    
  v.AL <- chisq.test(x = m.AL[,1],p = m.AL[,2], rescale.p = TRUE)

    #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 26 and 27     
  bialelicosNOimput.AL[g,26] <- v.AL$statistic
  bialelicosNOimput.AL[g,27] <- v.AL$p.value
  bialelicosNOimput.AL[g+1,26] <- v.AL$statistic
  bialelicosNOimput.AL[g+1,27] <- v.AL$p.value
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.AL , 
            "BialelicosNOImputados-AL.txt", 
            sep="\t",
            dec=",")

###################
#                 #
#     Imputed     #
#                 #   
###################

#We repeat the calculations for imputed SNP

#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/99 y 99/99 in column 1 
#and with the expected frequencies for 0/0, 0/99 y 99/99 in column 2

FrecuenciasImputados.AL <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0

while (h <= nrow(bialelicosimput.AL)) {
  
  #The uploaded allele frequency file presents consecutive lines for each SNP
  #in an unordered way: sometimes 0 allele comes first, and some others allele 99
  #We consider this issue to perform the calculations
  
  if(bialelicosimput.AL[h,5] == 0 && bialelicosimput.AL[h+1,5] == 99 ){     
  
  esp1 <- bialelicosimput.AL[h,6] * bialelicosimput.AL[h,6]
  esp2 <- bialelicosimput.AL[h,6] * bialelicosimput.AL[h+1,6] + bialelicosimput.AL[h+1,6]*bialelicosimput.AL[h,6]
  esp3 <- bialelicosimput.AL[h+1,6] * bialelicosimput.AL[h+1,6]
  
  obs1 <- as.double(paste(bialelicosimput.AL[h,20]))
  obs2 <- as.double(paste(bialelicosimput.AL[h,23]))
  obs3 <- as.double(paste(bialelicosimput.AL[h,25]))
  } 
  
  if(bialelicosimput.AL[h,5] == 99 && bialelicosimput.AL[h+1,5] == 0 ){     
    
    esp1 <- bialelicosimput.AL[h+1,6] * bialelicosimput.AL[h+1,6]
    esp2 <- bialelicosimput.AL[h+1,6] * bialelicosimput.AL[h,6] + bialelicosimput.AL[h,6]*bialelicosimput.AL[h+1,6]
    esp3 <- bialelicosimput.AL[h,6] * bialelicosimput.AL[h,6]
    
    obs1 <- as.double(paste(bialelicosimput.AL[h,20]))
    obs2 <- as.double(paste(bialelicosimput.AL[h,23]))
    obs3 <- as.double(paste(bialelicosimput.AL[h,25]))
  }  
  FrecuenciasImputados.AL[[a]] <- data.frame(x = c(obs1*98,obs2*98,obs3*98), 
                                             y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(98*obs1,98*obs2,98*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.AL <- 1

while (g < nrow(bialelicosimput.AL)) {
  
  m.AL <- FrecuenciasImputados.AL[[l]]
  
  #The uploaded allele frequency file presents consecutive lines for each SNP
  #in an unordered way: sometimes 0 allele comes first, and some others allele 99
  #We consider this issue to perform the calculations
  
  if(FrecuenciasImputados.AL[[l]][,2][1] == 0){
    FrecuenciasImputados.AL[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.AL <- FrecuenciasImputados.AL[[l]]
  }
  if(FrecuenciasImputados.AL[[l]][,2][2] == 0){
    FrecuenciasImputados.AL[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.AL <- FrecuenciasImputados.AL[[l]]
  }
  if(FrecuenciasImputados.AL[[l]][,2][3] == 0){
    FrecuenciasImputados.AL[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.AL <- FrecuenciasImputados.AL[[l]]
  }
  

  v.AL <- chisq.test(x = m.AL[,1],p = m.AL[,2], rescale.p = TRUE)

  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 26 and 27    
  bialelicosimput.AL[g,26] <- v.AL$statistic
  bialelicosimput.AL[g,27] <- v.AL$p.value
  bialelicosimput.AL[g+1,26] <- v.AL$statistic
  bialelicosimput.AL[g+1,27] <- v.AL$p.value
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.AL , 
            "BialelicosImputados-AL.txt", 
            sep="\t",
            dec=",")


#############################################################################################################
#END
#############################################################################################################



