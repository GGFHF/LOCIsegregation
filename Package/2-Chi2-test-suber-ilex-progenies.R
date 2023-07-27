##########################################################
#                                                        #
#          2-Chi2-test for suber-ilex progenies          #
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
#A05-FrecGenot-biallelic-sorted3.txt
#A07-FrecGenot-biallelic-sorted3.txt
#A09-FrecGenot-biallelic-sorted3.txt
#A10-FrecGenot-biallelic-sorted3.txt
#E28-FrecGenot-biallelic-sorted3.txt
#E31-FrecGenot-biallelic-sorted3.txt
#E41-FrecGenot-biallelic-sorted3.txt
#E96-FrecGenot-biallelic-sorted3.txt
#ScnI-MotherGenotypes-biallelic.txt


#We scan the above files to have data.frames for allele frequencies (common)
# and for observed genotypic frequencies obtained with NGShelper
#(one data.frame for each progeny, E28,E31,E41,E96: EN, 
#A05,A07,A09,A10: AL). Many SNP present missing data in all samples
#Here we only consider bi-allelic imputed and non-imputed loci
FrecuenciasAlelicas <- read.table("ScnI-FrecAlelicas-Bialelicos.txt", 
                                  header= TRUE, 
                                  sep = "\t", 
                                  dec = ".")

#Half-sib progenies of Q. suber pollinated by Q. suber
A05 <- read.table("A05-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")
A07 <- read.table("A07-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")
A09 <- read.table("A09-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")
A10 <- read.table("A10-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")
#Half-sib progenies of Q. ilex pollinated by Q. ilex
E28 <- read.table("E28-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")
E31 <- read.table("E31-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")
E41 <- read.table("E41-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")
E96 <- read.table("E96-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")

#We upload the table with the genotypes of the mothers

mother.genotypes<-read.table("ScnI-MotherGenotypes-biallelic.txt", header=TRUE,sep="\t")

#We set the working directory to produce the output.
setwd(WD.output.1)


#We generate a data.frame with sample size for each progeny
sample.size <- data.frame(71,94,101,103,82,84,91,93,46,46,30,67)
colnames(sample.size) <- c("A05","A07","A09","A10","E28","E31","E41","E96","FS16","FS19","FS20","FS22")

#########################
#                       # 
#          A05          #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.A05 <- A05 %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.A05 <- data.frame(FrecuenciasAlelicas$variant_id,
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
                             FrecGeno2.A05$X0.0,
                             FrecGeno2.A05$X0.1,
                             FrecGeno2.A05$X1.1,
                             FrecGeno2.A05$X0.99,
                             FrecGeno2.A05$X1.99,
                             FrecGeno2.A05$X99.99,
                             mother.genotypes2$AL.05)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.A05 <- subset(bialelic.A05, bialelic.A05$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.A05 <- subset(bialelic.A05, bialelic.A05$FrecuenciasAlelicas.imputations == "N")


###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2
FrecuenciasNOImputados.A05 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.A05)) {
if ((bialelicosNOimput.A05[h,20]+bialelicosNOimput.A05[h,21]+bialelicosNOimput.A05[h,22])==1){
  
  if(bialelicosNOimput.A05[h,26] == "0/0"){
  
  esp1 <- bialelicosNOimput.A05[h,6] 
  esp2 <- bialelicosNOimput.A05[h+1,6]
  esp3 <- 0
  
  obs1 <- as.double(paste(bialelicosNOimput.A05[h,20]))
  obs2 <- as.double(paste(bialelicosNOimput.A05[h,21]))
  obs3 <- as.double(paste(bialelicosNOimput.A05[h,22]))
  }
  
  if(bialelicosNOimput.A05[h,26]  == "0/1"){
    
    esp1 <- bialelicosNOimput.A05[h,6] * 0.5 
    esp2 <- bialelicosNOimput.A05[h,6] * 0.5 + bialelicosNOimput.A05[h+1,6] * 0.5
    esp3 <- bialelicosNOimput.A05[h+1,6] * 0.5
    
    obs1 <- as.double(paste(bialelicosNOimput.A05[h,20]))
    obs2 <- as.double(paste(bialelicosNOimput.A05[h,21]))
    obs3 <- as.double(paste(bialelicosNOimput.A05[h,22]))
  }
  
  if(bialelicosNOimput.A05[h,26]  == "1/1"){  
    esp1 <- 0
    esp2 <- bialelicosNOimput.A05[h,6] 
    esp3 <- bialelicosNOimput.A05[h+1,6]
    
    obs1 <- as.double(paste(bialelicosNOimput.A05[h,20]))
    obs2 <- as.double(paste(bialelicosNOimput.A05[h,21]))
    obs3 <- as.double(paste(bialelicosNOimput.A05[h,22]))
  }
  
  if(bialelicosNOimput.A05[h,26] == "./."){
    esp1 <- NA
    esp2 <- NA
    esp3 <- NA 
    
    obs1 <- as.double(paste(bialelicosNOimput.A05[h,20]))
    obs2 <- as.double(paste(bialelicosNOimput.A05[h,21]))
    obs3 <- as.double(paste(bialelicosNOimput.A05[h,22]))
  }
}  
else{
  esp1 <- NA
  esp2 <- NA
  esp3 <- NA
  
  obs1 <- NA
  obs2 <- NA
  obs3 <- NA
} 
   
  FrecuenciasNOImputados.A05[[a]] <- data.frame(x = c(obs1*sample.size$A05,obs2*sample.size$A05,obs3*sample.size$A05), 
                                              y = c(esp1,esp2,esp3))

  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$A05*obs1,sample.size$A05*obs2,sample.size$A05*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.A05 <- 1

while (g < nrow(bialelicosNOimput.A05)) {
  if(is.na(FrecuenciasNOImputados.A05[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.A05[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.A05[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.A05[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.A05[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.A05[[l]][,2][3])){
    
    bialelicosNOimput.A05[g,27] <- NA
    bialelicosNOimput.A05[g,28] <- NA
    bialelicosNOimput.A05[g+1,27] <- NA
    bialelicosNOimput.A05[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.A05 <- FrecuenciasNOImputados.A05[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE

  if(FrecuenciasNOImputados.A05[[l]][,1][1] == 0 && FrecuenciasNOImputados.A05[[l]][,2][1] == 0){
    FrecuenciasNOImputados.A05[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.A05 <- FrecuenciasNOImputados.A05[[l]]
  }
  if(FrecuenciasNOImputados.A05[[l]][,1][2] == 0 && FrecuenciasNOImputados.A05[[l]][,2][2] == 0){
    FrecuenciasNOImputados.A05[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.A05 <- FrecuenciasNOImputados.A05[[l]]
  }
  if(FrecuenciasNOImputados.A05[[l]][,1][3] == 0 && FrecuenciasNOImputados.A05[[l]][,2][3] == 0){
    FrecuenciasNOImputados.A05[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.A05 <- FrecuenciasNOImputados.A05[[l]]
  }
 
    v.A05 <- chisq.test(x = m.A05[,1],p = m.A05[,2], rescale.p = TRUE)
    
    #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28     
    bialelicosNOimput.A05[g,27] <- v.A05$statistic
    bialelicosNOimput.A05[g,28] <- v.A05$p.value
    bialelicosNOimput.A05[g+1,27] <- v.A05$statistic
    bialelicosNOimput.A05[g+1,28] <- v.A05$p.value

  g <- g+2
  l <- l+1
}


#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.A05 , 
            "BialelicosNOImputados-A05.txt", 
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

FrecuenciasImputados.A05 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0

#The uploaded allele frequency file presents consecutive lines for each SNP
#in an unordered way: sometimes 0 allele comes first, and some others allele 99
#We consider this issue to perform the calculations
#We include a condition to be sure that genotypic frequencies of each variant sum up 1

while (h <= nrow(bialelicosimput.A05)) {
  if(bialelicosimput.A05[h,5] == 99 && bialelicosimput.A05[h+1,5] == 0 ){  
    if ((bialelicosimput.A05[h,20]+bialelicosimput.A05[h,23]+bialelicosimput.A05[h,25])==1){
  
    if(bialelicosimput.A05[h,26] == "0/0"){
      
      esp1 <- bialelicosimput.A05[h+1,6] 
      esp2 <- bialelicosimput.A05[h,6]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosimput.A05[h,20]))
      obs2 <- as.double(paste(bialelicosimput.A05[h,23]))
      obs3 <- as.double(paste(bialelicosimput.A05[h,25]))
    }
    
    if(bialelicosimput.A05[h,26] == "0/99"){
      
      esp1 <- 0 
      esp2 <- bialelicosimput.A05[h+1,6] + bialelicosimput.A05[h,6] * 0.5 
      esp3 <- bialelicosimput.A05[h,6] * 0.5
      
      obs1 <- as.double(paste(bialelicosimput.A05[h,20]))
      obs2 <- as.double(paste(bialelicosimput.A05[h,23]))
      obs3 <- as.double(paste(bialelicosimput.A05[h,25]))
    }
    
    if(bialelicosimput.A05[h,26] == "99/99"){  
      esp1 <- 0
      esp2 <- bialelicosimput.A05[h+1,6] 
      esp3 <- bialelicosimput.A05[h,6]
      
      obs1 <- as.double(paste(bialelicosimput.A05[h,20]))
      obs2 <- as.double(paste(bialelicosimput.A05[h,23]))
      obs3 <- as.double(paste(bialelicosimput.A05[h,25]))
    }
    
    if(bialelicosimput.A05[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosimput.A05[h,20]))
      obs2 <- as.double(paste(bialelicosimput.A05[h,23]))
      obs3 <- as.double(paste(bialelicosimput.A05[h,25]))
    }
  }
  else{
    esp1 <- NA
    esp2 <- NA
    esp3 <- NA
    
    obs1 <- NA
    obs2 <- NA
    obs3 <- NA
  } 
}


  if(bialelicosimput.A05[h,5] == 0 && bialelicosimput.A05[h+1,5] == 99 ){  
    if ((bialelicosimput.A05[h,20]+bialelicosimput.A05[h,23]+bialelicosimput.A05[h,25])==1){
      
      if(bialelicosimput.A05[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.A05[h,6] 
        esp2 <- bialelicosimput.A05[h+1,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.A05[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A05[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A05[h,25]))
      }
      
      if(bialelicosimput.A05[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.A05[h,6] + bialelicosimput.A05[h+1,6] * 0.5 
        esp3 <- bialelicosimput.A05[h+1,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.A05[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A05[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A05[h,25]))
      }
      
      if(bialelicosimput.A05[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.A05[h,6] 
        esp3 <- bialelicosimput.A05[h+1,6]
        
        obs1 <- as.double(paste(bialelicosimput.A05[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A05[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A05[h,25]))
      }
      
      if(bialelicosimput.A05[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.A05[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A05[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A05[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
     } 
    }

    
  FrecuenciasImputados.A05[[a]] <- data.frame(x = c(obs1*sample.size$A05,obs2*sample.size$A05,obs3*sample.size$A05), 
                                                y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}


#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$A05*obs1,sample.size$A05*obs2,sample.size$A05*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.A05 <- 1

while (g < nrow(bialelicosimput.A05)) {
  if(is.na(FrecuenciasImputados.A05[[l]][,1][1]) || is.na(FrecuenciasImputados.A05[[l]][,1][2]) || is.na(FrecuenciasImputados.A05[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.A05[[l]][,2][1]) || is.na(FrecuenciasImputados.A05[[l]][,2][2]) || is.na(FrecuenciasImputados.A05[[l]][,2][3])){
    
    bialelicosimput.A05[g,27] <- NA
    bialelicosimput.A05[g,28] <- NA
    bialelicosimput.A05[g+1,27] <- NA
    bialelicosimput.A05[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.A05 <- FrecuenciasImputados.A05[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.A05[[l]][,1][1] == 0 && FrecuenciasImputados.A05[[l]][,2][1] == 0){
    FrecuenciasImputados.A05[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.A05 <- FrecuenciasImputados.A05[[l]]
  }
  if(FrecuenciasImputados.A05[[l]][,1][2] == 0 && FrecuenciasImputados.A05[[l]][,2][2] == 0){
    FrecuenciasImputados.A05[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.A05 <- FrecuenciasImputados.A05[[l]]
  }
  if(FrecuenciasImputados.A05[[l]][,1][3] == 0 && FrecuenciasImputados.A05[[l]][,2][3] == 0){
    FrecuenciasImputados.A05[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.A05 <- FrecuenciasImputados.A05[[l]]
  }
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28    
  
  v.A05 <- chisq.test(x = m.A05[,1],p = m.A05[,2], rescale.p = TRUE)
  
  bialelicosimput.A05[g,27] <- v.A05$statistic
  bialelicosimput.A05[g,28] <- v.A05$p.value
  bialelicosimput.A05[g+1,27] <- v.A05$statistic
  bialelicosimput.A05[g+1,28] <- v.A05$p.value
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.A05 , 
            "BialelicosImputados-A05.txt", 
            sep="\t",
            dec=",")


#####################################################################################


#########################
#                       # 
#          A07          #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.A07 <- A07 %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.A07 <- data.frame(FrecuenciasAlelicas$variant_id,
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
                           FrecGeno2.A07$X0.0,
                           FrecGeno2.A07$X0.1,
                           FrecGeno2.A07$X1.1,
                           FrecGeno2.A07$X0.99,
                           FrecGeno2.A07$X1.99,
                           FrecGeno2.A07$X99.99,
                           mother.genotypes2$AL.07)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.A07 <- subset(bialelic.A07, bialelic.A07$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.A07 <- subset(bialelic.A07, bialelic.A07$FrecuenciasAlelicas.imputations == "N")


###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2
FrecuenciasNOImputados.A07 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.A07)) {
  if ((bialelicosNOimput.A07[h,20]+bialelicosNOimput.A07[h,21]+bialelicosNOimput.A07[h,22])==1){
    
    if(bialelicosNOimput.A07[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.A07[h,6] 
      esp2 <- bialelicosNOimput.A07[h+1,6]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.A07[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A07[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A07[h,22]))
    }
    
    if(bialelicosNOimput.A07[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.A07[h,6] * 0.5 
      esp2 <- bialelicosNOimput.A07[h,6] * 0.5 + bialelicosNOimput.A07[h+1,6] * 0.5
      esp3 <- bialelicosNOimput.A07[h+1,6] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.A07[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A07[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A07[h,22]))
    }
    
    if(bialelicosNOimput.A07[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.A07[h,6] 
      esp3 <- bialelicosNOimput.A07[h+1,6]
      
      obs1 <- as.double(paste(bialelicosNOimput.A07[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A07[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A07[h,22]))
    }
    
    if(bialelicosNOimput.A07[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.A07[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A07[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A07[h,22]))
    }
  }  
  else{
    esp1 <- NA
    esp2 <- NA
    esp3 <- NA
    
    obs1 <- NA
    obs2 <- NA
    obs3 <- NA
  } 
  
  FrecuenciasNOImputados.A07[[a]] <- data.frame(x = c(obs1*sample.size$A07,obs2*sample.size$A07,obs3*sample.size$A07), 
                                                y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$A07*obs1,sample.size$A07*obs2,sample.size$A07*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.A07 <- 1

while (g < nrow(bialelicosNOimput.A07)) {
  if(is.na(FrecuenciasNOImputados.A07[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.A07[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.A07[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.A07[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.A07[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.A07[[l]][,2][3])){
    
    bialelicosNOimput.A07[g,27] <- NA
    bialelicosNOimput.A07[g,28] <- NA
    bialelicosNOimput.A07[g+1,27] <- NA
    bialelicosNOimput.A07[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.A07 <- FrecuenciasNOImputados.A07[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.A07[[l]][,1][1] == 0 && FrecuenciasNOImputados.A07[[l]][,2][1] == 0){
    FrecuenciasNOImputados.A07[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.A07 <- FrecuenciasNOImputados.A07[[l]]
  }
  if(FrecuenciasNOImputados.A07[[l]][,1][2] == 0 && FrecuenciasNOImputados.A07[[l]][,2][2] == 0){
    FrecuenciasNOImputados.A07[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.A07 <- FrecuenciasNOImputados.A07[[l]]
  }
  if(FrecuenciasNOImputados.A07[[l]][,1][3] == 0 && FrecuenciasNOImputados.A07[[l]][,2][3] == 0){
    FrecuenciasNOImputados.A07[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.A07 <- FrecuenciasNOImputados.A07[[l]]
  }
  
  v.A07 <- chisq.test(x = m.A07[,1],p = m.A07[,2], rescale.p = TRUE)
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28     
  bialelicosNOimput.A07[g,27] <- v.A07$statistic
  bialelicosNOimput.A07[g,28] <- v.A07$p.value
  bialelicosNOimput.A07[g+1,27] <- v.A07$statistic
  bialelicosNOimput.A07[g+1,28] <- v.A07$p.value
  
  g <- g+2
  l <- l+1
}


#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.A07 , 
            "BialelicosNOImputados-A07.txt", 
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

FrecuenciasImputados.A07 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0

#The uploaded allele frequency file presents consecutive lines for each SNP
#in an unordered way: sometimes 0 allele comes first, and some others allele 99
#We consider this issue to perform the calculations
#We include a condition to be sure that genotypic frequencies of each variant sum up 1

while (h <= nrow(bialelicosimput.A07)) {
  if(bialelicosimput.A07[h,5] == 99 && bialelicosimput.A07[h+1,5] == 0 ){  
    if ((bialelicosimput.A07[h,20]+bialelicosimput.A07[h,23]+bialelicosimput.A07[h,25])==1){
      
      if(bialelicosimput.A07[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.A07[h+1,6] 
        esp2 <- bialelicosimput.A07[h,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.A07[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A07[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A07[h,25]))
      }
      
      if(bialelicosimput.A07[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.A07[h+1,6] + bialelicosimput.A07[h,6] * 0.5 
        esp3 <- bialelicosimput.A07[h,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.A07[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A07[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A07[h,25]))
      }
      
      if(bialelicosimput.A07[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.A07[h+1,6] 
        esp3 <- bialelicosimput.A07[h,6]
        
        obs1 <- as.double(paste(bialelicosimput.A07[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A07[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A07[h,25]))
      }
      
      if(bialelicosimput.A07[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.A07[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A07[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A07[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  if(bialelicosimput.A07[h,5] == 0 && bialelicosimput.A07[h+1,5] == 99 ){  
    if ((bialelicosimput.A07[h,20]+bialelicosimput.A07[h,23]+bialelicosimput.A07[h,25])==1){
      
      if(bialelicosimput.A07[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.A07[h,6] 
        esp2 <- bialelicosimput.A07[h+1,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.A07[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A07[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A07[h,25]))
      }
      
      if(bialelicosimput.A07[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.A07[h,6] + bialelicosimput.A07[h+1,6] * 0.5 
        esp3 <- bialelicosimput.A07[h+1,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.A07[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A07[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A07[h,25]))
      }
      
      if(bialelicosimput.A07[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.A07[h,6] 
        esp3 <- bialelicosimput.A07[h+1,6]
        
        obs1 <- as.double(paste(bialelicosimput.A07[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A07[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A07[h,25]))
      }
      
      if(bialelicosimput.A07[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.A07[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A07[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A07[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  FrecuenciasImputados.A07[[a]] <- data.frame(x = c(obs1*sample.size$A07,obs2*sample.size$A07,obs3*sample.size$A07), 
                                              y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}


#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$A07*obs1,sample.size$A07*obs2,sample.size$A07*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.A07 <- 1

while (g < nrow(bialelicosimput.A07)) {
  if(is.na(FrecuenciasImputados.A07[[l]][,1][1]) || is.na(FrecuenciasImputados.A07[[l]][,1][2]) || is.na(FrecuenciasImputados.A07[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.A07[[l]][,2][1]) || is.na(FrecuenciasImputados.A07[[l]][,2][2]) || is.na(FrecuenciasImputados.A07[[l]][,2][3])){
    
    bialelicosimput.A07[g,27] <- NA
    bialelicosimput.A07[g,28] <- NA
    bialelicosimput.A07[g+1,27] <- NA
    bialelicosimput.A07[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.A07 <- FrecuenciasImputados.A07[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.A07[[l]][,1][1] == 0 && FrecuenciasImputados.A07[[l]][,2][1] == 0){
    FrecuenciasImputados.A07[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.A07 <- FrecuenciasImputados.A07[[l]]
  }
  if(FrecuenciasImputados.A07[[l]][,1][2] == 0 && FrecuenciasImputados.A07[[l]][,2][2] == 0){
    FrecuenciasImputados.A07[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.A07 <- FrecuenciasImputados.A07[[l]]
  }
  if(FrecuenciasImputados.A07[[l]][,1][3] == 0 && FrecuenciasImputados.A07[[l]][,2][3] == 0){
    FrecuenciasImputados.A07[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.A07 <- FrecuenciasImputados.A07[[l]]
  }
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28    
  
  v.A07 <- chisq.test(x = m.A07[,1],p = m.A07[,2], rescale.p = TRUE)
  
  bialelicosimput.A07[g,27] <- v.A07$statistic
  bialelicosimput.A07[g,28] <- v.A07$p.value
  bialelicosimput.A07[g+1,27] <- v.A07$statistic
  bialelicosimput.A07[g+1,28] <- v.A07$p.value
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.A07 , 
            "BialelicosImputados-A07.txt", 
            sep="\t",
            dec=",")


#####################################################################################

#########################
#                       # 
#          A09          #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.A09 <- A09 %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.A09 <- data.frame(FrecuenciasAlelicas$variant_id,
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
                           FrecGeno2.A09$X0.0,
                           FrecGeno2.A09$X0.1,
                           FrecGeno2.A09$X1.1,
                           FrecGeno2.A09$X0.99,
                           FrecGeno2.A09$X1.99,
                           FrecGeno2.A09$X99.99,
                           mother.genotypes2$AL.09)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.A09 <- subset(bialelic.A09, bialelic.A09$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.A09 <- subset(bialelic.A09, bialelic.A09$FrecuenciasAlelicas.imputations == "N")


###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2
FrecuenciasNOImputados.A09 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.A09)) {
  if ((bialelicosNOimput.A09[h,20]+bialelicosNOimput.A09[h,21]+bialelicosNOimput.A09[h,22])==1){
    
    if(bialelicosNOimput.A09[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.A09[h,6] 
      esp2 <- bialelicosNOimput.A09[h+1,6]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.A09[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A09[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A09[h,22]))
    }
    
    if(bialelicosNOimput.A09[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.A09[h,6] * 0.5 
      esp2 <- bialelicosNOimput.A09[h,6] * 0.5 + bialelicosNOimput.A09[h+1,6] * 0.5
      esp3 <- bialelicosNOimput.A09[h+1,6] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.A09[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A09[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A09[h,22]))
    }
    
    if(bialelicosNOimput.A09[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.A09[h,6] 
      esp3 <- bialelicosNOimput.A09[h+1,6]
      
      obs1 <- as.double(paste(bialelicosNOimput.A09[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A09[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A09[h,22]))
    }
    
    if(bialelicosNOimput.A09[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.A09[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A09[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A09[h,22]))
    }
  }  
  else{
    esp1 <- NA
    esp2 <- NA
    esp3 <- NA
    
    obs1 <- NA
    obs2 <- NA
    obs3 <- NA
  } 
  
  FrecuenciasNOImputados.A09[[a]] <- data.frame(x = c(obs1*sample.size$A09,obs2*sample.size$A09,obs3*sample.size$A09), 
                                                y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$A09*obs1,sample.size$A09*obs2,sample.size$A09*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.A09 <- 1

while (g < nrow(bialelicosNOimput.A09)) {
  if(is.na(FrecuenciasNOImputados.A09[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.A09[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.A09[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.A09[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.A09[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.A09[[l]][,2][3])){
    
    bialelicosNOimput.A09[g,27] <- NA
    bialelicosNOimput.A09[g,28] <- NA
    bialelicosNOimput.A09[g+1,27] <- NA
    bialelicosNOimput.A09[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.A09 <- FrecuenciasNOImputados.A09[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.A09[[l]][,1][1] == 0 && FrecuenciasNOImputados.A09[[l]][,2][1] == 0){
    FrecuenciasNOImputados.A09[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.A09 <- FrecuenciasNOImputados.A09[[l]]
  }
  if(FrecuenciasNOImputados.A09[[l]][,1][2] == 0 && FrecuenciasNOImputados.A09[[l]][,2][2] == 0){
    FrecuenciasNOImputados.A09[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.A09 <- FrecuenciasNOImputados.A09[[l]]
  }
  if(FrecuenciasNOImputados.A09[[l]][,1][3] == 0 && FrecuenciasNOImputados.A09[[l]][,2][3] == 0){
    FrecuenciasNOImputados.A09[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.A09 <- FrecuenciasNOImputados.A09[[l]]
  }
  
  v.A09 <- chisq.test(x = m.A09[,1],p = m.A09[,2], rescale.p = TRUE)
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28     
  bialelicosNOimput.A09[g,27] <- v.A09$statistic
  bialelicosNOimput.A09[g,28] <- v.A09$p.value
  bialelicosNOimput.A09[g+1,27] <- v.A09$statistic
  bialelicosNOimput.A09[g+1,28] <- v.A09$p.value
  
  g <- g+2
  l <- l+1
}


#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.A09 , 
            "BialelicosNOImputados-A09.txt", 
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

FrecuenciasImputados.A09 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0

#The uploaded allele frequency file presents consecutive lines for each SNP
#in an unordered way: sometimes 0 allele comes first, and some others allele 99
#We consider this issue to perform the calculations
#We include a condition to be sure that genotypic frequencies of each variant sum up 1

while (h <= nrow(bialelicosimput.A09)) {
  if(bialelicosimput.A09[h,5] == 99 && bialelicosimput.A09[h+1,5] == 0 ){  
    if ((bialelicosimput.A09[h,20]+bialelicosimput.A09[h,23]+bialelicosimput.A09[h,25])==1){
      
      if(bialelicosimput.A09[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.A09[h+1,6] 
        esp2 <- bialelicosimput.A09[h,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.A09[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A09[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A09[h,25]))
      }
      
      if(bialelicosimput.A09[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.A09[h+1,6] + bialelicosimput.A09[h,6] * 0.5 
        esp3 <- bialelicosimput.A09[h,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.A09[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A09[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A09[h,25]))
      }
      
      if(bialelicosimput.A09[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.A09[h+1,6] 
        esp3 <- bialelicosimput.A09[h,6]
        
        obs1 <- as.double(paste(bialelicosimput.A09[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A09[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A09[h,25]))
      }
      
      if(bialelicosimput.A09[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.A09[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A09[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A09[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  if(bialelicosimput.A09[h,5] == 0 && bialelicosimput.A09[h+1,5] == 99 ){  
    if ((bialelicosimput.A09[h,20]+bialelicosimput.A09[h,23]+bialelicosimput.A09[h,25])==1){
      
      if(bialelicosimput.A09[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.A09[h,6] 
        esp2 <- bialelicosimput.A09[h+1,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.A09[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A09[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A09[h,25]))
      }
      
      if(bialelicosimput.A09[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.A09[h,6] + bialelicosimput.A09[h+1,6] * 0.5 
        esp3 <- bialelicosimput.A09[h+1,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.A09[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A09[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A09[h,25]))
      }
      
      if(bialelicosimput.A09[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.A09[h,6] 
        esp3 <- bialelicosimput.A09[h+1,6]
        
        obs1 <- as.double(paste(bialelicosimput.A09[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A09[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A09[h,25]))
      }
      
      if(bialelicosimput.A09[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.A09[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A09[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A09[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  FrecuenciasImputados.A09[[a]] <- data.frame(x = c(obs1*sample.size$A09,obs2*sample.size$A09,obs3*sample.size$A09), 
                                              y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}


#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$A09*obs1,sample.size$A09*obs2,sample.size$A09*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.A09 <- 1

while (g < nrow(bialelicosimput.A09)) {
  if(is.na(FrecuenciasImputados.A09[[l]][,1][1]) || is.na(FrecuenciasImputados.A09[[l]][,1][2]) || is.na(FrecuenciasImputados.A09[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.A09[[l]][,2][1]) || is.na(FrecuenciasImputados.A09[[l]][,2][2]) || is.na(FrecuenciasImputados.A09[[l]][,2][3])){
    
    bialelicosimput.A09[g,27] <- NA
    bialelicosimput.A09[g,28] <- NA
    bialelicosimput.A09[g+1,27] <- NA
    bialelicosimput.A09[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.A09 <- FrecuenciasImputados.A09[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.A09[[l]][,1][1] == 0 && FrecuenciasImputados.A09[[l]][,2][1] == 0){
    FrecuenciasImputados.A09[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.A09 <- FrecuenciasImputados.A09[[l]]
  }
  if(FrecuenciasImputados.A09[[l]][,1][2] == 0 && FrecuenciasImputados.A09[[l]][,2][2] == 0){
    FrecuenciasImputados.A09[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.A09 <- FrecuenciasImputados.A09[[l]]
  }
  if(FrecuenciasImputados.A09[[l]][,1][3] == 0 && FrecuenciasImputados.A09[[l]][,2][3] == 0){
    FrecuenciasImputados.A09[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.A09 <- FrecuenciasImputados.A09[[l]]
  }
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28    
  
  v.A09 <- chisq.test(x = m.A09[,1],p = m.A09[,2], rescale.p = TRUE)
  
  bialelicosimput.A09[g,27] <- v.A09$statistic
  bialelicosimput.A09[g,28] <- v.A09$p.value
  bialelicosimput.A09[g+1,27] <- v.A09$statistic
  bialelicosimput.A09[g+1,28] <- v.A09$p.value
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.A09 , 
            "BialelicosImputados-A09.txt", 
            sep="\t",
            dec=",")


#####################################################################################


#########################
#                       # 
#          A10          #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.A10 <- A10 %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.A10 <- data.frame(FrecuenciasAlelicas$variant_id,
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
                           FrecGeno2.A10$X0.0,
                           FrecGeno2.A10$X0.1,
                           FrecGeno2.A10$X1.1,
                           FrecGeno2.A10$X0.99,
                           FrecGeno2.A10$X1.99,
                           FrecGeno2.A10$X99.99,
                           mother.genotypes2$AL.10)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.A10 <- subset(bialelic.A10, bialelic.A10$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.A10 <- subset(bialelic.A10, bialelic.A10$FrecuenciasAlelicas.imputations == "N")


###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2
FrecuenciasNOImputados.A10 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.A10)) {
  if ((bialelicosNOimput.A10[h,20]+bialelicosNOimput.A10[h,21]+bialelicosNOimput.A10[h,22])==1){
    
    if(bialelicosNOimput.A10[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.A10[h,6] 
      esp2 <- bialelicosNOimput.A10[h+1,6]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.A10[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A10[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A10[h,22]))
    }
    
    if(bialelicosNOimput.A10[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.A10[h,6] * 0.5 
      esp2 <- bialelicosNOimput.A10[h,6] * 0.5 + bialelicosNOimput.A10[h+1,6] * 0.5
      esp3 <- bialelicosNOimput.A10[h+1,6] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.A10[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A10[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A10[h,22]))
    }
    
    if(bialelicosNOimput.A10[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.A10[h,6] 
      esp3 <- bialelicosNOimput.A10[h+1,6]
      
      obs1 <- as.double(paste(bialelicosNOimput.A10[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A10[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A10[h,22]))
    }
    
    if(bialelicosNOimput.A10[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.A10[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.A10[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.A10[h,22]))
    }
  }  
  else{
    esp1 <- NA
    esp2 <- NA
    esp3 <- NA
    
    obs1 <- NA
    obs2 <- NA
    obs3 <- NA
  } 
  
  FrecuenciasNOImputados.A10[[a]] <- data.frame(x = c(obs1*sample.size$A10,obs2*sample.size$A10,obs3*sample.size$A10), 
                                                y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$A10*obs1,sample.size$A10*obs2,sample.size$A10*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.A10 <- 1

while (g < nrow(bialelicosNOimput.A10)) {
  if(is.na(FrecuenciasNOImputados.A10[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.A10[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.A10[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.A10[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.A10[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.A10[[l]][,2][3])){
    
    bialelicosNOimput.A10[g,27] <- NA
    bialelicosNOimput.A10[g,28] <- NA
    bialelicosNOimput.A10[g+1,27] <- NA
    bialelicosNOimput.A10[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.A10 <- FrecuenciasNOImputados.A10[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.A10[[l]][,1][1] == 0 && FrecuenciasNOImputados.A10[[l]][,2][1] == 0){
    FrecuenciasNOImputados.A10[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.A10 <- FrecuenciasNOImputados.A10[[l]]
  }
  if(FrecuenciasNOImputados.A10[[l]][,1][2] == 0 && FrecuenciasNOImputados.A10[[l]][,2][2] == 0){
    FrecuenciasNOImputados.A10[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.A10 <- FrecuenciasNOImputados.A10[[l]]
  }
  if(FrecuenciasNOImputados.A10[[l]][,1][3] == 0 && FrecuenciasNOImputados.A10[[l]][,2][3] == 0){
    FrecuenciasNOImputados.A10[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.A10 <- FrecuenciasNOImputados.A10[[l]]
  }
  
  v.A10 <- chisq.test(x = m.A10[,1],p = m.A10[,2], rescale.p = TRUE)
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28     
  bialelicosNOimput.A10[g,27] <- v.A10$statistic
  bialelicosNOimput.A10[g,28] <- v.A10$p.value
  bialelicosNOimput.A10[g+1,27] <- v.A10$statistic
  bialelicosNOimput.A10[g+1,28] <- v.A10$p.value
  
  g <- g+2
  l <- l+1
}


#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.A10 , 
            "BialelicosNOImputados-A10.txt", 
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

FrecuenciasImputados.A10 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0

#The uploaded allele frequency file presents consecutive lines for each SNP
#in an unordered way: sometimes 0 allele comes first, and some others allele 99
#We consider this issue to perform the calculations
#We include a condition to be sure that genotypic frequencies of each variant sum up 1

while (h <= nrow(bialelicosimput.A10)) {
  if(bialelicosimput.A10[h,5] == 99 && bialelicosimput.A10[h+1,5] == 0 ){  
    if ((bialelicosimput.A10[h,20]+bialelicosimput.A10[h,23]+bialelicosimput.A10[h,25])==1){
      
      if(bialelicosimput.A10[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.A10[h+1,6] 
        esp2 <- bialelicosimput.A10[h,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.A10[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A10[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A10[h,25]))
      }
      
      if(bialelicosimput.A10[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.A10[h+1,6] + bialelicosimput.A10[h,6] * 0.5 
        esp3 <- bialelicosimput.A10[h,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.A10[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A10[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A10[h,25]))
      }
      
      if(bialelicosimput.A10[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.A10[h+1,6] 
        esp3 <- bialelicosimput.A10[h,6]
        
        obs1 <- as.double(paste(bialelicosimput.A10[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A10[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A10[h,25]))
      }
      
      if(bialelicosimput.A10[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.A10[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A10[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A10[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  if(bialelicosimput.A10[h,5] == 0 && bialelicosimput.A10[h+1,5] == 99 ){  
    if ((bialelicosimput.A10[h,20]+bialelicosimput.A10[h,23]+bialelicosimput.A10[h,25])==1){
      
      if(bialelicosimput.A10[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.A10[h,6] 
        esp2 <- bialelicosimput.A10[h+1,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.A10[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A10[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A10[h,25]))
      }
      
      if(bialelicosimput.A10[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.A10[h,6] + bialelicosimput.A10[h+1,6] * 0.5 
        esp3 <- bialelicosimput.A10[h+1,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.A10[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A10[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A10[h,25]))
      }
      
      if(bialelicosimput.A10[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.A10[h,6] 
        esp3 <- bialelicosimput.A10[h+1,6]
        
        obs1 <- as.double(paste(bialelicosimput.A10[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A10[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A10[h,25]))
      }
      
      if(bialelicosimput.A10[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.A10[h,20]))
        obs2 <- as.double(paste(bialelicosimput.A10[h,23]))
        obs3 <- as.double(paste(bialelicosimput.A10[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  FrecuenciasImputados.A10[[a]] <- data.frame(x = c(obs1*sample.size$A10,obs2*sample.size$A10,obs3*sample.size$A10), 
                                              y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}


#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$A10*obs1,sample.size$A10*obs2,sample.size$A10*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.A10 <- 1

while (g < nrow(bialelicosimput.A10)) {
  if(is.na(FrecuenciasImputados.A10[[l]][,1][1]) || is.na(FrecuenciasImputados.A10[[l]][,1][2]) || is.na(FrecuenciasImputados.A10[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.A10[[l]][,2][1]) || is.na(FrecuenciasImputados.A10[[l]][,2][2]) || is.na(FrecuenciasImputados.A10[[l]][,2][3])){
    
    bialelicosimput.A10[g,27] <- NA
    bialelicosimput.A10[g,28] <- NA
    bialelicosimput.A10[g+1,27] <- NA
    bialelicosimput.A10[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.A10 <- FrecuenciasImputados.A10[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.A10[[l]][,1][1] == 0 && FrecuenciasImputados.A10[[l]][,2][1] == 0){
    FrecuenciasImputados.A10[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.A10 <- FrecuenciasImputados.A10[[l]]
  }
  if(FrecuenciasImputados.A10[[l]][,1][2] == 0 && FrecuenciasImputados.A10[[l]][,2][2] == 0){
    FrecuenciasImputados.A10[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.A10 <- FrecuenciasImputados.A10[[l]]
  }
  if(FrecuenciasImputados.A10[[l]][,1][3] == 0 && FrecuenciasImputados.A10[[l]][,2][3] == 0){
    FrecuenciasImputados.A10[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.A10 <- FrecuenciasImputados.A10[[l]]
  }
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28    
  
  v.A10 <- chisq.test(x = m.A10[,1],p = m.A10[,2], rescale.p = TRUE)
  
  bialelicosimput.A10[g,27] <- v.A10$statistic
  bialelicosimput.A10[g,28] <- v.A10$p.value
  bialelicosimput.A10[g+1,27] <- v.A10$statistic
  bialelicosimput.A10[g+1,28] <- v.A10$p.value
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.A10 , 
            "BialelicosImputados-A10.txt", 
            sep="\t",
            dec=",")


#####################################################################################

#########################
#                       # 
#          E28          #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.E28 <- E28 %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.E28 <- data.frame(FrecuenciasAlelicas$variant_id,
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
                           FrecGeno2.E28$X0.0,
                           FrecGeno2.E28$X0.1,
                           FrecGeno2.E28$X1.1,
                           FrecGeno2.E28$X0.99,
                           FrecGeno2.E28$X1.99,
                           FrecGeno2.E28$X99.99,
                           mother.genotypes2$EN.28)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.E28 <- subset(bialelic.E28, bialelic.E28$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.E28 <- subset(bialelic.E28, bialelic.E28$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2

FrecuenciasNOImputados.E28 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.E28)) {
  if ((bialelicosNOimput.E28[h,20]+bialelicosNOimput.E28[h,21]+bialelicosNOimput.E28[h,22])==1){
    
    if(bialelicosNOimput.E28[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.E28[h,7] 
      esp2 <- bialelicosNOimput.E28[h+1,7]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.E28[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E28[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E28[h,22]))
    }
    
    if(bialelicosNOimput.E28[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.E28[h,7] * 0.5 
      esp2 <- bialelicosNOimput.E28[h,7] * 0.5 + bialelicosNOimput.E28[h+1,7] * 0.5
      esp3 <- bialelicosNOimput.E28[h+1,7] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.E28[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E28[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E28[h,22]))
    }
    
    if(bialelicosNOimput.E28[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.E28[h,7] 
      esp3 <- bialelicosNOimput.E28[h+1,7]
      
      obs1 <- as.double(paste(bialelicosNOimput.E28[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E28[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E28[h,22]))
    }
    
    if(bialelicosNOimput.E28[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.E28[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E28[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E28[h,22]))
    }
  }  
  else{
    esp1 <- NA
    esp2 <- NA
    esp3 <- NA
    
    obs1 <- NA
    obs2 <- NA
    obs3 <- NA
  } 
  
  FrecuenciasNOImputados.E28[[a]] <- data.frame(x = c(obs1*sample.size$E28,obs2*sample.size$E28,obs3*sample.size$E28), 
                                                y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$E28*obs1,sample.size$E28*obs2,sample.size$E28*obs3)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.E28 <- 1

while (g < nrow(bialelicosNOimput.E28)) {

  if(is.na(FrecuenciasNOImputados.E28[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.E28[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.E28[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.E28[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.E28[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.E28[[l]][,2][3])){
    
    bialelicosNOimput.E28[g,27] <- NA
    bialelicosNOimput.E28[g,28] <- NA
    bialelicosNOimput.E28[g+1,27] <- NA
    bialelicosNOimput.E28[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.E28 <- FrecuenciasNOImputados.E28[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.E28[[l]][,1][1] == 0 && FrecuenciasNOImputados.E28[[l]][,2][1] == 0){
    FrecuenciasNOImputados.E28[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.E28 <- FrecuenciasNOImputados.E28[[l]]
  }
  if(FrecuenciasNOImputados.E28[[l]][,1][2] == 0 && FrecuenciasNOImputados.E28[[l]][,2][2] == 0){
    FrecuenciasNOImputados.E28[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.E28 <- FrecuenciasNOImputados.E28[[l]]
  }
  if(FrecuenciasNOImputados.E28[[l]][,1][3] == 0 && FrecuenciasNOImputados.E28[[l]][,2][3] == 0){
    FrecuenciasNOImputados.E28[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.E28 <- FrecuenciasNOImputados.E28[[l]]
  }
  v.E28 <- chisq.test(x = m.E28[,1],p = m.E28[,2], rescale.p = TRUE)

  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28     
  bialelicosNOimput.E28[g,27] <- v.E28$statistic
  bialelicosNOimput.E28[g,28] <- v.E28$p.value
  bialelicosNOimput.E28[g+1,27] <- v.E28$statistic
  bialelicosNOimput.E28[g+1,28] <- v.E28$p.value
  
  
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.E28 , 
            "BialelicosNOImputados-E28.txt", 
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

FrecuenciasImputados.E28 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#The uploaded allele frequency file presents consecutive lines for each SNP
#in an unordered way: sometimes 0 allele comes first, and some others allele 99
#We consider this issue to perform the calculations
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosimput.E28)) {
  if(bialelicosimput.E28[h,5] == 99 && bialelicosimput.E28[h+1,5] == 0 ){  
    if ((bialelicosimput.E28[h,20]+bialelicosimput.E28[h,23]+bialelicosimput.E28[h,25])==1){
      
      if(bialelicosimput.E28[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.E28[h+1,7] 
        esp2 <- bialelicosimput.E28[h,7]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.E28[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E28[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E28[h,25]))
      }
      
      if(bialelicosimput.E28[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.E28[h+1,7] + bialelicosimput.E28[h,7] * 0.5 
        esp3 <- bialelicosimput.E28[h,7] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.E28[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E28[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E28[h,25]))
      }
      
      if(bialelicosimput.E28[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.E28[h+1,7] 
        esp3 <- bialelicosimput.E28[h,7]
        
        obs1 <- as.double(paste(bialelicosimput.E28[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E28[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E28[h,25]))
      }
      
      if(bialelicosimput.E28[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.E28[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E28[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E28[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  if(bialelicosimput.E28[h,5] == 0 && bialelicosimput.E28[h+1,5] == 99 ){  
    if ((bialelicosimput.E28[h,20]+bialelicosimput.E28[h,23]+bialelicosimput.E28[h,25])==1){
      
      if(bialelicosimput.E28[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.E28[h,7] 
        esp2 <- bialelicosimput.E28[h+1,7]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.E28[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E28[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E28[h,25]))
      }
      
      if(bialelicosimput.E28[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.E28[h,7] + bialelicosimput.E28[h+1,7] * 0.5 
        esp3 <- bialelicosimput.E28[h+1,7] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.E28[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E28[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E28[h,25]))
      }
      
      if(bialelicosimput.E28[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.E28[h,7] 
        esp3 <- bialelicosimput.E28[h+1,7]
        
        obs1 <- as.double(paste(bialelicosimput.E28[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E28[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E28[h,25]))
      }
      
      if(bialelicosimput.E28[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.E28[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E28[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E28[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  FrecuenciasImputados.E28[[a]] <- data.frame(x = c(obs1*sample.size$E28,obs2*sample.size$E28,obs3*sample.size$E28), 
                                              y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$E28*obs1,sample.size$E28*obs2,sample.size$E28*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.E28 <- 1

while (g < nrow(bialelicosimput.E28)) {
  
  
  
  
  if(is.na(FrecuenciasImputados.E28[[l]][,1][1]) || is.na(FrecuenciasImputados.E28[[l]][,1][2]) || is.na(FrecuenciasImputados.E28[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.E28[[l]][,2][1]) || is.na(FrecuenciasImputados.E28[[l]][,2][2]) || is.na(FrecuenciasImputados.E28[[l]][,2][3])){
    
    bialelicosimput.E28[g,27] <- NA
    bialelicosimput.E28[g,28] <- NA
    bialelicosimput.E28[g+1,27] <- NA
    bialelicosimput.E28[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.E28 <- FrecuenciasImputados.E28[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.E28[[l]][,1][1] == 0 && FrecuenciasImputados.E28[[l]][,2][1] == 0){
    FrecuenciasImputados.E28[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.E28 <- FrecuenciasImputados.E28[[l]]
  }
  if(FrecuenciasImputados.E28[[l]][,1][2] == 0 && FrecuenciasImputados.E28[[l]][,2][2] == 0){
    FrecuenciasImputados.E28[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.E28 <- FrecuenciasImputados.E28[[l]]
  }
  if(FrecuenciasImputados.E28[[l]][,1][3] == 0 && FrecuenciasImputados.E28[[l]][,2][3] == 0){
    FrecuenciasImputados.E28[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.E28 <- FrecuenciasImputados.E28[[l]]
  }

  v.E28 <- chisq.test(x = m.E28[,1],p = m.E28[,2], rescale.p = TRUE)
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28    
  bialelicosimput.E28[g,27] <- v.E28$statistic
  bialelicosimput.E28[g,28] <- v.E28$p.value
  bialelicosimput.E28[g+1,27] <- v.E28$statistic
  bialelicosimput.E28[g+1,28] <- v.E28$p.value
  
  
  
  g <- g+2
  l <- l+1
}

write.table(bialelicosimput.E28 , "BialelicosImputados-E28.txt")
#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.E28 , 
            "BialelicosImputados-E28..txt", 
            sep="\t",
            dec=",")

#####################################################################################


#########################
#                       # 
#          E31          #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.E31 <- E31 %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.E31 <- data.frame(FrecuenciasAlelicas$variant_id,
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
                           FrecGeno2.E31$X0.0,
                           FrecGeno2.E31$X0.1,
                           FrecGeno2.E31$X1.1,
                           FrecGeno2.E31$X0.99,
                           FrecGeno2.E31$X1.99,
                           FrecGeno2.E31$X99.99,
                           mother.genotypes2$EN.31)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.E31 <- subset(bialelic.E31, bialelic.E31$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.E31 <- subset(bialelic.E31, bialelic.E31$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2

FrecuenciasNOImputados.E31 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.E31)) {
  if ((bialelicosNOimput.E31[h,20]+bialelicosNOimput.E31[h,21]+bialelicosNOimput.E31[h,22])==1){
    
    if(bialelicosNOimput.E31[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.E31[h,7] 
      esp2 <- bialelicosNOimput.E31[h+1,7]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.E31[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E31[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E31[h,22]))
    }
    
    if(bialelicosNOimput.E31[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.E31[h,7] * 0.5 
      esp2 <- bialelicosNOimput.E31[h,7] * 0.5 + bialelicosNOimput.E31[h+1,7] * 0.5
      esp3 <- bialelicosNOimput.E31[h+1,7] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.E31[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E31[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E31[h,22]))
    }
    
    if(bialelicosNOimput.E31[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.E31[h,7] 
      esp3 <- bialelicosNOimput.E31[h+1,7]
      
      obs1 <- as.double(paste(bialelicosNOimput.E31[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E31[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E31[h,22]))
    }
    
    if(bialelicosNOimput.E31[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.E31[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E31[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E31[h,22]))
    }
  }  
  else{
    esp1 <- NA
    esp2 <- NA
    esp3 <- NA
    
    obs1 <- NA
    obs2 <- NA
    obs3 <- NA
  } 
  
  FrecuenciasNOImputados.E31[[a]] <- data.frame(x = c(obs1*sample.size$E31,obs2*sample.size$E31,obs3*sample.size$E31), 
                                                y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$E31*obs1,sample.size$E31*obs2,sample.size$E31*obs3)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.E31 <- 1

while (g < nrow(bialelicosNOimput.E31)) {
  
  if(is.na(FrecuenciasNOImputados.E31[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.E31[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.E31[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.E31[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.E31[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.E31[[l]][,2][3])){
    
    bialelicosNOimput.E31[g,27] <- NA
    bialelicosNOimput.E31[g,28] <- NA
    bialelicosNOimput.E31[g+1,27] <- NA
    bialelicosNOimput.E31[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.E31 <- FrecuenciasNOImputados.E31[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.E31[[l]][,1][1] == 0 && FrecuenciasNOImputados.E31[[l]][,2][1] == 0){
    FrecuenciasNOImputados.E31[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.E31 <- FrecuenciasNOImputados.E31[[l]]
  }
  if(FrecuenciasNOImputados.E31[[l]][,1][2] == 0 && FrecuenciasNOImputados.E31[[l]][,2][2] == 0){
    FrecuenciasNOImputados.E31[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.E31 <- FrecuenciasNOImputados.E31[[l]]
  }
  if(FrecuenciasNOImputados.E31[[l]][,1][3] == 0 && FrecuenciasNOImputados.E31[[l]][,2][3] == 0){
    FrecuenciasNOImputados.E31[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.E31 <- FrecuenciasNOImputados.E31[[l]]
  }
  v.E31 <- chisq.test(x = m.E31[,1],p = m.E31[,2], rescale.p = TRUE)
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28     
  bialelicosNOimput.E31[g,27] <- v.E31$statistic
  bialelicosNOimput.E31[g,28] <- v.E31$p.value
  bialelicosNOimput.E31[g+1,27] <- v.E31$statistic
  bialelicosNOimput.E31[g+1,28] <- v.E31$p.value
  
  
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.E31 , 
            "BialelicosNOImputados-E31.txt", 
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

FrecuenciasImputados.E31 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#The uploaded allele frequency file presents consecutive lines for each SNP
#in an unordered way: sometimes 0 allele comes first, and some others allele 99
#We consider this issue to perform the calculations
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosimput.E31)) {
  if(bialelicosimput.E31[h,5] == 99 && bialelicosimput.E31[h+1,5] == 0 ){  
    if ((bialelicosimput.E31[h,20]+bialelicosimput.E31[h,23]+bialelicosimput.E31[h,25])==1){
      
      if(bialelicosimput.E31[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.E31[h+1,7] 
        esp2 <- bialelicosimput.E31[h,7]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.E31[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E31[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E31[h,25]))
      }
      
      if(bialelicosimput.E31[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.E31[h+1,7] + bialelicosimput.E31[h,7] * 0.5 
        esp3 <- bialelicosimput.E31[h,7] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.E31[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E31[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E31[h,25]))
      }
      
      if(bialelicosimput.E31[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.E31[h+1,7] 
        esp3 <- bialelicosimput.E31[h,7]
        
        obs1 <- as.double(paste(bialelicosimput.E31[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E31[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E31[h,25]))
      }
      
      if(bialelicosimput.E31[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.E31[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E31[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E31[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  if(bialelicosimput.E31[h,5] == 0 && bialelicosimput.E31[h+1,5] == 99 ){  
    if ((bialelicosimput.E31[h,20]+bialelicosimput.E31[h,23]+bialelicosimput.E31[h,25])==1){
      
      if(bialelicosimput.E31[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.E31[h,7] 
        esp2 <- bialelicosimput.E31[h+1,7]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.E31[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E31[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E31[h,25]))
      }
      
      if(bialelicosimput.E31[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.E31[h,7] + bialelicosimput.E31[h+1,7] * 0.5 
        esp3 <- bialelicosimput.E31[h+1,7] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.E31[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E31[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E31[h,25]))
      }
      
      if(bialelicosimput.E31[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.E31[h,7] 
        esp3 <- bialelicosimput.E31[h+1,7]
        
        obs1 <- as.double(paste(bialelicosimput.E31[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E31[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E31[h,25]))
      }
      
      if(bialelicosimput.E31[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.E31[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E31[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E31[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  FrecuenciasImputados.E31[[a]] <- data.frame(x = c(obs1*sample.size$E31,obs2*sample.size$E31,obs3*sample.size$E31), 
                                              y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$E31*obs1,sample.size$E31*obs2,sample.size$E31*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.E31 <- 1

while (g < nrow(bialelicosimput.E31)) {
  
  
  
  
  if(is.na(FrecuenciasImputados.E31[[l]][,1][1]) || is.na(FrecuenciasImputados.E31[[l]][,1][2]) || is.na(FrecuenciasImputados.E31[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.E31[[l]][,2][1]) || is.na(FrecuenciasImputados.E31[[l]][,2][2]) || is.na(FrecuenciasImputados.E31[[l]][,2][3])){
    
    bialelicosimput.E31[g,27] <- NA
    bialelicosimput.E31[g,28] <- NA
    bialelicosimput.E31[g+1,27] <- NA
    bialelicosimput.E31[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.E31 <- FrecuenciasImputados.E31[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.E31[[l]][,1][1] == 0 && FrecuenciasImputados.E31[[l]][,2][1] == 0){
    FrecuenciasImputados.E31[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.E31 <- FrecuenciasImputados.E31[[l]]
  }
  if(FrecuenciasImputados.E31[[l]][,1][2] == 0 && FrecuenciasImputados.E31[[l]][,2][2] == 0){
    FrecuenciasImputados.E31[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.E31 <- FrecuenciasImputados.E31[[l]]
  }
  if(FrecuenciasImputados.E31[[l]][,1][3] == 0 && FrecuenciasImputados.E31[[l]][,2][3] == 0){
    FrecuenciasImputados.E31[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.E31 <- FrecuenciasImputados.E31[[l]]
  }
  
  v.E31 <- chisq.test(x = m.E31[,1],p = m.E31[,2], rescale.p = TRUE)
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28    
  bialelicosimput.E31[g,27] <- v.E31$statistic
  bialelicosimput.E31[g,28] <- v.E31$p.value
  bialelicosimput.E31[g+1,27] <- v.E31$statistic
  bialelicosimput.E31[g+1,28] <- v.E31$p.value
  
  
  
  g <- g+2
  l <- l+1
}

write.table(bialelicosimput.E31 , "BialelicosImputados-E31.txt")
#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.E31 , 
            "BialelicosImputados-E31..txt", 
            sep="\t",
            dec=",")

#####################################################################################

#########################
#                       # 
#          E41          #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.E41 <- E41 %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.E41 <- data.frame(FrecuenciasAlelicas$variant_id,
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
                           FrecGeno2.E41$X0.0,
                           FrecGeno2.E41$X0.1,
                           FrecGeno2.E41$X1.1,
                           FrecGeno2.E41$X0.99,
                           FrecGeno2.E41$X1.99,
                           FrecGeno2.E41$X99.99,
                           mother.genotypes2$EN.41)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.E41 <- subset(bialelic.E41, bialelic.E41$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.E41 <- subset(bialelic.E41, bialelic.E41$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2

FrecuenciasNOImputados.E41 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.E41)) {
  if ((bialelicosNOimput.E41[h,20]+bialelicosNOimput.E41[h,21]+bialelicosNOimput.E41[h,22])==1){
    
    if(bialelicosNOimput.E41[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.E41[h,7] 
      esp2 <- bialelicosNOimput.E41[h+1,7]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.E41[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E41[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E41[h,22]))
    }
    
    if(bialelicosNOimput.E41[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.E41[h,7] * 0.5 
      esp2 <- bialelicosNOimput.E41[h,7] * 0.5 + bialelicosNOimput.E41[h+1,7] * 0.5
      esp3 <- bialelicosNOimput.E41[h+1,7] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.E41[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E41[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E41[h,22]))
    }
    
    if(bialelicosNOimput.E41[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.E41[h,7] 
      esp3 <- bialelicosNOimput.E41[h+1,7]
      
      obs1 <- as.double(paste(bialelicosNOimput.E41[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E41[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E41[h,22]))
    }
    
    if(bialelicosNOimput.E41[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.E41[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E41[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E41[h,22]))
    }
  }  
  else{
    esp1 <- NA
    esp2 <- NA
    esp3 <- NA
    
    obs1 <- NA
    obs2 <- NA
    obs3 <- NA
  } 
  
  FrecuenciasNOImputados.E41[[a]] <- data.frame(x = c(obs1*sample.size$E41,obs2*sample.size$E41,obs3*sample.size$E41), 
                                                y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$E41*obs1,sample.size$E41*obs2,sample.size$E41*obs3)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.E41 <- 1

while (g < nrow(bialelicosNOimput.E41)) {
  
  if(is.na(FrecuenciasNOImputados.E41[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.E41[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.E41[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.E41[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.E41[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.E41[[l]][,2][3])){
    
    bialelicosNOimput.E41[g,27] <- NA
    bialelicosNOimput.E41[g,28] <- NA
    bialelicosNOimput.E41[g+1,27] <- NA
    bialelicosNOimput.E41[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.E41 <- FrecuenciasNOImputados.E41[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.E41[[l]][,1][1] == 0 && FrecuenciasNOImputados.E41[[l]][,2][1] == 0){
    FrecuenciasNOImputados.E41[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.E41 <- FrecuenciasNOImputados.E41[[l]]
  }
  if(FrecuenciasNOImputados.E41[[l]][,1][2] == 0 && FrecuenciasNOImputados.E41[[l]][,2][2] == 0){
    FrecuenciasNOImputados.E41[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.E41 <- FrecuenciasNOImputados.E41[[l]]
  }
  if(FrecuenciasNOImputados.E41[[l]][,1][3] == 0 && FrecuenciasNOImputados.E41[[l]][,2][3] == 0){
    FrecuenciasNOImputados.E41[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.E41 <- FrecuenciasNOImputados.E41[[l]]
  }
  v.E41 <- chisq.test(x = m.E41[,1],p = m.E41[,2], rescale.p = TRUE)
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28     
  bialelicosNOimput.E41[g,27] <- v.E41$statistic
  bialelicosNOimput.E41[g,28] <- v.E41$p.value
  bialelicosNOimput.E41[g+1,27] <- v.E41$statistic
  bialelicosNOimput.E41[g+1,28] <- v.E41$p.value
  
  
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.E41 , 
            "BialelicosNOImputados-E41.txt", 
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

FrecuenciasImputados.E41 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#The uploaded allele frequency file presents consecutive lines for each SNP
#in an unordered way: sometimes 0 allele comes first, and some others allele 99
#We consider this issue to perform the calculations
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosimput.E41)) {
  if(bialelicosimput.E41[h,5] == 99 && bialelicosimput.E41[h+1,5] == 0 ){  
    if ((bialelicosimput.E41[h,20]+bialelicosimput.E41[h,23]+bialelicosimput.E41[h,25])==1){
      
      if(bialelicosimput.E41[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.E41[h+1,7] 
        esp2 <- bialelicosimput.E41[h,7]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.E41[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E41[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E41[h,25]))
      }
      
      if(bialelicosimput.E41[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.E41[h+1,7] + bialelicosimput.E41[h,7] * 0.5 
        esp3 <- bialelicosimput.E41[h,7] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.E41[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E41[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E41[h,25]))
      }
      
      if(bialelicosimput.E41[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.E41[h+1,7] 
        esp3 <- bialelicosimput.E41[h,7]
        
        obs1 <- as.double(paste(bialelicosimput.E41[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E41[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E41[h,25]))
      }
      
      if(bialelicosimput.E41[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.E41[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E41[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E41[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  if(bialelicosimput.E41[h,5] == 0 && bialelicosimput.E41[h+1,5] == 99 ){  
    if ((bialelicosimput.E41[h,20]+bialelicosimput.E41[h,23]+bialelicosimput.E41[h,25])==1){
      
      if(bialelicosimput.E41[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.E41[h,7] 
        esp2 <- bialelicosimput.E41[h+1,7]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.E41[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E41[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E41[h,25]))
      }
      
      if(bialelicosimput.E41[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.E41[h,7] + bialelicosimput.E41[h+1,7] * 0.5 
        esp3 <- bialelicosimput.E41[h+1,7] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.E41[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E41[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E41[h,25]))
      }
      
      if(bialelicosimput.E41[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.E41[h,7] 
        esp3 <- bialelicosimput.E41[h+1,7]
        
        obs1 <- as.double(paste(bialelicosimput.E41[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E41[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E41[h,25]))
      }
      
      if(bialelicosimput.E41[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.E41[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E41[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E41[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  FrecuenciasImputados.E41[[a]] <- data.frame(x = c(obs1*sample.size$E41,obs2*sample.size$E41,obs3*sample.size$E41), 
                                              y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$E41*obs1,sample.size$E41*obs2,sample.size$E41*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.E41 <- 1

while (g < nrow(bialelicosimput.E41)) {
  
  
  
  
  if(is.na(FrecuenciasImputados.E41[[l]][,1][1]) || is.na(FrecuenciasImputados.E41[[l]][,1][2]) || is.na(FrecuenciasImputados.E41[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.E41[[l]][,2][1]) || is.na(FrecuenciasImputados.E41[[l]][,2][2]) || is.na(FrecuenciasImputados.E41[[l]][,2][3])){
    
    bialelicosimput.E41[g,27] <- NA
    bialelicosimput.E41[g,28] <- NA
    bialelicosimput.E41[g+1,27] <- NA
    bialelicosimput.E41[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.E41 <- FrecuenciasImputados.E41[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.E41[[l]][,1][1] == 0 && FrecuenciasImputados.E41[[l]][,2][1] == 0){
    FrecuenciasImputados.E41[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.E41 <- FrecuenciasImputados.E41[[l]]
  }
  if(FrecuenciasImputados.E41[[l]][,1][2] == 0 && FrecuenciasImputados.E41[[l]][,2][2] == 0){
    FrecuenciasImputados.E41[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.E41 <- FrecuenciasImputados.E41[[l]]
  }
  if(FrecuenciasImputados.E41[[l]][,1][3] == 0 && FrecuenciasImputados.E41[[l]][,2][3] == 0){
    FrecuenciasImputados.E41[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.E41 <- FrecuenciasImputados.E41[[l]]
  }
  
  v.E41 <- chisq.test(x = m.E41[,1],p = m.E41[,2], rescale.p = TRUE)
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28    
  bialelicosimput.E41[g,27] <- v.E41$statistic
  bialelicosimput.E41[g,28] <- v.E41$p.value
  bialelicosimput.E41[g+1,27] <- v.E41$statistic
  bialelicosimput.E41[g+1,28] <- v.E41$p.value
  
  
  
  g <- g+2
  l <- l+1
}

write.table(bialelicosimput.E41 , "BialelicosImputados-E41.txt")
#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.E41 , 
            "BialelicosImputados-E41..txt", 
            sep="\t",
            dec=",")

#####################################################################################

#########################
#                       # 
#          E96          #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.E96 <- E96 %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.E96 <- data.frame(FrecuenciasAlelicas$variant_id,
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
                           FrecGeno2.E96$X0.0,
                           FrecGeno2.E96$X0.1,
                           FrecGeno2.E96$X1.1,
                           FrecGeno2.E96$X0.99,
                           FrecGeno2.E96$X1.99,
                           FrecGeno2.E96$X99.99,
                           mother.genotypes2$EN.96)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.E96 <- subset(bialelic.E96, bialelic.E96$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.E96 <- subset(bialelic.E96, bialelic.E96$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################


#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2

FrecuenciasNOImputados.E96 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.E96)) {
  if ((bialelicosNOimput.E96[h,20]+bialelicosNOimput.E96[h,21]+bialelicosNOimput.E96[h,22])==1){
    
    if(bialelicosNOimput.E96[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.E96[h,7] 
      esp2 <- bialelicosNOimput.E96[h+1,7]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.E96[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E96[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E96[h,22]))
    }
    
    if(bialelicosNOimput.E96[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.E96[h,7] * 0.5 
      esp2 <- bialelicosNOimput.E96[h,7] * 0.5 + bialelicosNOimput.E96[h+1,7] * 0.5
      esp3 <- bialelicosNOimput.E96[h+1,7] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.E96[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E96[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E96[h,22]))
    }
    
    if(bialelicosNOimput.E96[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.E96[h,7] 
      esp3 <- bialelicosNOimput.E96[h+1,7]
      
      obs1 <- as.double(paste(bialelicosNOimput.E96[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E96[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E96[h,22]))
    }
    
    if(bialelicosNOimput.E96[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.E96[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.E96[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.E96[h,22]))
    }
  }  
  else{
    esp1 <- NA
    esp2 <- NA
    esp3 <- NA
    
    obs1 <- NA
    obs2 <- NA
    obs3 <- NA
  } 
  
  FrecuenciasNOImputados.E96[[a]] <- data.frame(x = c(obs1*sample.size$E96,obs2*sample.size$E96,obs3*sample.size$E96), 
                                                y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$E96*obs1,sample.size$E96*obs2,sample.size$E96*obs3)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.E96 <- 1

while (g < nrow(bialelicosNOimput.E96)) {
  
  if(is.na(FrecuenciasNOImputados.E96[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.E96[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.E96[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.E96[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.E96[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.E96[[l]][,2][3])){
    
    bialelicosNOimput.E96[g,27] <- NA
    bialelicosNOimput.E96[g,28] <- NA
    bialelicosNOimput.E96[g+1,27] <- NA
    bialelicosNOimput.E96[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.E96 <- FrecuenciasNOImputados.E96[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.E96[[l]][,1][1] == 0 && FrecuenciasNOImputados.E96[[l]][,2][1] == 0){
    FrecuenciasNOImputados.E96[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.E96 <- FrecuenciasNOImputados.E96[[l]]
  }
  if(FrecuenciasNOImputados.E96[[l]][,1][2] == 0 && FrecuenciasNOImputados.E96[[l]][,2][2] == 0){
    FrecuenciasNOImputados.E96[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.E96 <- FrecuenciasNOImputados.E96[[l]]
  }
  if(FrecuenciasNOImputados.E96[[l]][,1][3] == 0 && FrecuenciasNOImputados.E96[[l]][,2][3] == 0){
    FrecuenciasNOImputados.E96[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.E96 <- FrecuenciasNOImputados.E96[[l]]
  }
  v.E96 <- chisq.test(x = m.E96[,1],p = m.E96[,2], rescale.p = TRUE)
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28     
  bialelicosNOimput.E96[g,27] <- v.E96$statistic
  bialelicosNOimput.E96[g,28] <- v.E96$p.value
  bialelicosNOimput.E96[g+1,27] <- v.E96$statistic
  bialelicosNOimput.E96[g+1,28] <- v.E96$p.value
  
  
  
  g <- g+2
  l <- l+1
}

#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.E96 , 
            "BialelicosNOImputados-E96.txt", 
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

FrecuenciasImputados.E96 <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#The uploaded allele frequency file presents consecutive lines for each SNP
#in an unordered way: sometimes 0 allele comes first, and some others allele 99
#We consider this issue to perform the calculations
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosimput.E96)) {
  if(bialelicosimput.E96[h,5] == 99 && bialelicosimput.E96[h+1,5] == 0 ){  
    if ((bialelicosimput.E96[h,20]+bialelicosimput.E96[h,23]+bialelicosimput.E96[h,25])==1){
      
      if(bialelicosimput.E96[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.E96[h+1,7] 
        esp2 <- bialelicosimput.E96[h,7]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.E96[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E96[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E96[h,25]))
      }
      
      if(bialelicosimput.E96[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.E96[h+1,7] + bialelicosimput.E96[h,7] * 0.5 
        esp3 <- bialelicosimput.E96[h,7] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.E96[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E96[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E96[h,25]))
      }
      
      if(bialelicosimput.E96[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.E96[h+1,7] 
        esp3 <- bialelicosimput.E96[h,7]
        
        obs1 <- as.double(paste(bialelicosimput.E96[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E96[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E96[h,25]))
      }
      
      if(bialelicosimput.E96[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.E96[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E96[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E96[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  if(bialelicosimput.E96[h,5] == 0 && bialelicosimput.E96[h+1,5] == 99 ){  
    if ((bialelicosimput.E96[h,20]+bialelicosimput.E96[h,23]+bialelicosimput.E96[h,25])==1){
      
      if(bialelicosimput.E96[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.E96[h,7] 
        esp2 <- bialelicosimput.E96[h+1,7]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.E96[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E96[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E96[h,25]))
      }
      
      if(bialelicosimput.E96[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.E96[h,7] + bialelicosimput.E96[h+1,7] * 0.5 
        esp3 <- bialelicosimput.E96[h+1,7] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.E96[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E96[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E96[h,25]))
      }
      
      if(bialelicosimput.E96[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.E96[h,7] 
        esp3 <- bialelicosimput.E96[h+1,7]
        
        obs1 <- as.double(paste(bialelicosimput.E96[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E96[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E96[h,25]))
      }
      
      if(bialelicosimput.E96[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.E96[h,20]))
        obs2 <- as.double(paste(bialelicosimput.E96[h,23]))
        obs3 <- as.double(paste(bialelicosimput.E96[h,25]))
      }
    }
    else{
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA
      
      obs1 <- NA
      obs2 <- NA
      obs3 <- NA
    } 
  }
  
  
  FrecuenciasImputados.E96[[a]] <- data.frame(x = c(obs1*sample.size$E96,obs2*sample.size$E96,obs3*sample.size$E96), 
                                              y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$E96*obs1,sample.size$E96*obs2,sample.size$E96*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.E96 <- 1

while (g < nrow(bialelicosimput.E96)) {
  
  
  
  
  if(is.na(FrecuenciasImputados.E96[[l]][,1][1]) || is.na(FrecuenciasImputados.E96[[l]][,1][2]) || is.na(FrecuenciasImputados.E96[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.E96[[l]][,2][1]) || is.na(FrecuenciasImputados.E96[[l]][,2][2]) || is.na(FrecuenciasImputados.E96[[l]][,2][3])){
    
    bialelicosimput.E96[g,27] <- NA
    bialelicosimput.E96[g,28] <- NA
    bialelicosimput.E96[g+1,27] <- NA
    bialelicosimput.E96[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.E96 <- FrecuenciasImputados.E96[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.E96[[l]][,1][1] == 0 && FrecuenciasImputados.E96[[l]][,2][1] == 0){
    FrecuenciasImputados.E96[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.E96 <- FrecuenciasImputados.E96[[l]]
  }
  if(FrecuenciasImputados.E96[[l]][,1][2] == 0 && FrecuenciasImputados.E96[[l]][,2][2] == 0){
    FrecuenciasImputados.E96[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.E96 <- FrecuenciasImputados.E96[[l]]
  }
  if(FrecuenciasImputados.E96[[l]][,1][3] == 0 && FrecuenciasImputados.E96[[l]][,2][3] == 0){
    FrecuenciasImputados.E96[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.E96 <- FrecuenciasImputados.E96[[l]]
  }
  
  v.E96 <- chisq.test(x = m.E96[,1],p = m.E96[,2], rescale.p = TRUE)
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28    
  bialelicosimput.E96[g,27] <- v.E96$statistic
  bialelicosimput.E96[g,28] <- v.E96$p.value
  bialelicosimput.E96[g+1,27] <- v.E96$statistic
  bialelicosimput.E96[g+1,28] <- v.E96$p.value
  
  
  
  g <- g+2
  l <- l+1
}

write.table(bialelicosimput.E96 , "BialelicosImputados-E96.txt")
#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.E96 , 
            "BialelicosImputados-E96..txt", 
            sep="\t",
            dec=",")

#####################################################################################
##################
#                #
#       END      #
#                #
##################




