##########################################################
#                                                        #
#          3-Chi2-test for hybrid progenies          #
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
#FS16_S-FrecGenot-biallelic-sorted3.txt
#FS19_S-FrecGenot-biallelic-sorted3.txt
#FS20_S-FrecGenot-biallelic-sorted3.txt
#FS22_S-FrecGenot-biallelic-sorted3.txt
#ScnI-MotherGenotypes-biallelic.txt


#We scan the above files to have data.frames for allele frequencies (common)
# and for observed genotypic frequencies obtained with NGShelper
#(one data.frame for each progeny, FS16,FS19,FS20,FS22: hybrids 
#Many SNP present missing data in all samples
#Here we only consider bi-allelic imputed and non-imputed loci

FrecuenciasAlelicas <- read.table("ScnI-FrecAlelicas-Bialelicos.txt", 
                                  header= TRUE, 
                                  sep = "\t", 
                                  dec = ".")
#Half-sib progenies of hybrids pollinated by Q. suber
FS16.S <- read.table("FS16_S-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")
FS19.S <- read.table("FS19_S-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")
FS20.S <- read.table("FS20_S-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")
FS22.S <- read.table("FS22_S-FrecGenot-biallelic-sorted3.txt", header= TRUE, sep = "\t", dec = ".")

#We upload the table with the genotypes of the mothers
mother.genotypes<-read.table("ScnI-MotherGenotypes-biallelic.txt", header=TRUE,sep="\t")


#We set the working directory to produce the output.
setwd(WD.output.1)

#We generate a data.frame with sample size for each progeny
sample.size <- data.frame(71,94,101,103,82,84,91,93,46,46,30,67)
colnames(sample.size) <- c("A05","A07","A09","A10","E28","E31","E41","E96","FS16","FS19","FS20","FS22")

#########################
#                       # 
#          FS16.S       #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.FS16.S <- FS16.S %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.FS16.S <- data.frame(FrecuenciasAlelicas$variant_id,
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
                              FrecGeno2.FS16.S$X0.0,
                              FrecGeno2.FS16.S$X0.1,
                              FrecGeno2.FS16.S$X1.1,
                              FrecGeno2.FS16.S$X0.99,
                              FrecGeno2.FS16.S$X1.99,
                              FrecGeno2.FS16.S$X99.99,
                              mother.genotypes2$EFS.16)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.FS16.S <- subset(bialelic.FS16.S, bialelic.FS16.S$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.FS16.S <- subset(bialelic.FS16.S, bialelic.FS16.S$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################

#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2
FrecuenciasNOImputados.FS16.S <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.FS16.S)) {
  if ((bialelicosNOimput.FS16.S[h,20]+bialelicosNOimput.FS16.S[h,21]+bialelicosNOimput.FS16.S[h,22])==1){
    
    if(bialelicosNOimput.FS16.S[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.FS16.S[h,6] 
      esp2 <- bialelicosNOimput.FS16.S[h+1,6]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.FS16.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS16.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS16.S[h,22]))
    }
    
    if(bialelicosNOimput.FS16.S[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.FS16.S[h,6] * 0.5 
      esp2 <- bialelicosNOimput.FS16.S[h,6] * 0.5 + bialelicosNOimput.FS16.S[h+1,6] * 0.5
      esp3 <- bialelicosNOimput.FS16.S[h+1,6] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.FS16.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS16.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS16.S[h,22]))
    }
    
    if(bialelicosNOimput.FS16.S[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.FS16.S[h,6] 
      esp3 <- bialelicosNOimput.FS16.S[h+1,6]
      
      obs1 <- as.double(paste(bialelicosNOimput.FS16.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS16.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS16.S[h,22]))
    }
    
    if(bialelicosNOimput.FS16.S[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.FS16.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS16.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS16.S[h,22]))
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
  
  FrecuenciasNOImputados.FS16.S[[a]] <- data.frame(x = c(obs1*sample.size$FS16,obs2*sample.size$FS16,obs3*sample.size$FS16), 
                                                   y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$FS16*obs1,sample.size$FS16*obs2,sample.size$FS16*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.FS16.S <- 1

while (g < nrow(bialelicosNOimput.FS16.S)) {
 
  if(is.na(FrecuenciasNOImputados.FS16.S[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.FS16.S[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.FS16.S[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.FS16.S[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.FS16.S[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.FS16.S[[l]][,2][3])){
    
    bialelicosNOimput.FS16.S[g,27] <- NA
    bialelicosNOimput.FS16.S[g,28] <- NA
    bialelicosNOimput.FS16.S[g+1,27] <- NA
    bialelicosNOimput.FS16.S[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.FS16.S <- FrecuenciasNOImputados.FS16.S[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.FS16.S[[l]][,1][1] == 0 && FrecuenciasNOImputados.FS16.S[[l]][,2][1] == 0){
    FrecuenciasNOImputados.FS16.S[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.FS16.S <- FrecuenciasNOImputados.FS16.S[[l]]
  }
  if(FrecuenciasNOImputados.FS16.S[[l]][,1][2] == 0 && FrecuenciasNOImputados.FS16.S[[l]][,2][2] == 0){
    FrecuenciasNOImputados.FS16.S[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.FS16.S <- FrecuenciasNOImputados.FS16.S[[l]]
  }
  if(FrecuenciasNOImputados.FS16.S[[l]][,1][3] == 0 && FrecuenciasNOImputados.FS16.S[[l]][,2][3] == 0){
    FrecuenciasNOImputados.FS16.S[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.FS16.S <- FrecuenciasNOImputados.FS16.S[[l]]
  }

  v.FS16.S <- chisq.test(x = m.FS16.S[,1],p = m.FS16.S[,2], rescale.p = TRUE)

  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28   
  bialelicosNOimput.FS16.S[g,27] <- v.FS16.S$statistic
  bialelicosNOimput.FS16.S[g,28] <- v.FS16.S$p.value
  bialelicosNOimput.FS16.S[g+1,27] <- v.FS16.S$statistic
  bialelicosNOimput.FS16.S[g+1,28] <- v.FS16.S$p.value

  g <- g+2
  l <- l+1
}


#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.FS16.S, 
            "BialelicosNOImputados-FS16-S.txt", 
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
FrecuenciasImputados.FS16.S <- list()
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

while (h <= nrow(bialelicosimput.FS16.S)) {
  if(bialelicosimput.FS16.S[h,5] == 99 && bialelicosimput.FS16.S[h+1,5] == 0 ){  
    if ((bialelicosimput.FS16.S[h,20]+bialelicosimput.FS16.S[h,23]+bialelicosimput.FS16.S[h,25])==1){
      
      if(bialelicosimput.FS16.S[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.FS16.S[h+1,6] 
        esp2 <- bialelicosimput.FS16.S[h,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.FS16.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS16.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS16.S[h,25]))
      }
      
      if(bialelicosimput.FS16.S[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.FS16.S[h+1,6] + bialelicosimput.FS16.S[h,6] * 0.5 
        esp3 <- bialelicosimput.FS16.S[h,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.FS16.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS16.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS16.S[h,25]))
      }
      
      if(bialelicosimput.FS16.S[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.FS16.S[h+1,6] 
        esp3 <- bialelicosimput.FS16.S[h,6]
        
        obs1 <- as.double(paste(bialelicosimput.FS16.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS16.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS16.S[h,25]))
      }
      
      if(bialelicosimput.FS16.S[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.FS16.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS16.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS16.S[h,25]))
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
  
  
  if(bialelicosimput.FS16.S[h,5] == 0 && bialelicosimput.FS16.S[h+1,5] == 99 ){  
    if ((bialelicosimput.FS16.S[h,20]+bialelicosimput.FS16.S[h,23]+bialelicosimput.FS16.S[h,25])==1){
      
      if(bialelicosimput.FS16.S[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.FS16.S[h,6] 
        esp2 <- bialelicosimput.FS16.S[h+1,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.FS16.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS16.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS16.S[h,25]))
      }
      
      if(bialelicosimput.FS16.S[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.FS16.S[h,6] + bialelicosimput.FS16.S[h+1,6] * 0.5 
        esp3 <- bialelicosimput.FS16.S[h+1,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.FS16.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS16.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS16.S[h,25]))
      }
      
      if(bialelicosimput.FS16.S[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.FS16.S[h,6] 
        esp3 <- bialelicosimput.FS16.S[h+1,6]
        
        obs1 <- as.double(paste(bialelicosimput.FS16.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS16.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS16.S[h,25]))
      }
      
      if(bialelicosimput.FS16.S[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.FS16.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS16.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS16.S[h,25]))
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
  
  
  FrecuenciasImputados.FS16.S[[a]] <- data.frame(x = c(obs1*sample.size$FS16,obs2*sample.size$FS16,obs3*sample.size$FS16), 
                                                 y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$FS16*obs1,sample.size$FS16*obs2,sample.size$FS16*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.FS16.S <- 1

while (g < nrow(bialelicosimput.FS16.S)) {
  
  if(is.na(FrecuenciasImputados.FS16.S[[l]][,1][1]) 
     || is.na(FrecuenciasImputados.FS16.S[[l]][,1][2]) 
     || is.na(FrecuenciasImputados.FS16.S[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.FS16.S[[l]][,2][1]) 
     || is.na(FrecuenciasImputados.FS16.S[[l]][,2][2]) 
     || is.na(FrecuenciasImputados.FS16.S[[l]][,2][3])){
    
    bialelicosimput.FS16.S[g,27] <- NA
    bialelicosimput.FS16.S[g,28] <- NA
    bialelicosimput.FS16.S[g+1,27] <- NA
    bialelicosimput.FS16.S[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.FS16.S <- FrecuenciasImputados.FS16.S[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.FS16.S[[l]][,1][1] == 0 && FrecuenciasImputados.FS16.S[[l]][,2][1] == 0){
    FrecuenciasImputados.FS16.S[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.FS16.S <- FrecuenciasImputados.FS16.S[[l]]
  }
  if(FrecuenciasImputados.FS16.S[[l]][,1][2] == 0 && FrecuenciasImputados.FS16.S[[l]][,2][2] == 0){
    FrecuenciasImputados.FS16.S[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.FS16.S <- FrecuenciasImputados.FS16.S[[l]]
  }
  if(FrecuenciasImputados.FS16.S[[l]][,1][3] == 0 && FrecuenciasImputados.FS16.S[[l]][,2][3] == 0){
    FrecuenciasImputados.FS16.S[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.FS16.S <- FrecuenciasImputados.FS16.S[[l]]
  }
  
  v.FS16.S <- chisq.test(x = m.FS16.S[,1],p = m.FS16.S[,2], rescale.p = TRUE)
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28  
  
  bialelicosimput.FS16.S[g,27] <- v.FS16.S$statistic
  bialelicosimput.FS16.S[g,28] <- v.FS16.S$p.value
  bialelicosimput.FS16.S[g+1,27] <- v.FS16.S$statistic
  bialelicosimput.FS16.S[g+1,28] <- v.FS16.S$p.value
  
  
  
  g <- g+2
  l <- l+1
}

write.table(bialelicosimput.FS16.S , "BialelicosImputados-FS16-S.txt")
#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.FS16.S , 
            "BialelicosImputados-FS16-S.txt", 
            sep="\t",
            dec=",")

#####################################################################################



#########################
#                       # 
#          FS19.S       #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.FS19.S <- FS19.S %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.FS19.S <- data.frame(FrecuenciasAlelicas$variant_id,
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
                              FrecGeno2.FS19.S$X0.0,
                              FrecGeno2.FS19.S$X0.1,
                              FrecGeno2.FS19.S$X1.1,
                              FrecGeno2.FS19.S$X0.99,
                              FrecGeno2.FS19.S$X1.99,
                              FrecGeno2.FS19.S$X99.99,
                              mother.genotypes2$EFS.19)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.FS19.S <- subset(bialelic.FS19.S, bialelic.FS19.S$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.FS19.S <- subset(bialelic.FS19.S, bialelic.FS19.S$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################

#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2
FrecuenciasNOImputados.FS19.S <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.FS19.S)) {
  if ((bialelicosNOimput.FS19.S[h,20]+bialelicosNOimput.FS19.S[h,21]+bialelicosNOimput.FS19.S[h,22])==1){
    
    if(bialelicosNOimput.FS19.S[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.FS19.S[h,6] 
      esp2 <- bialelicosNOimput.FS19.S[h+1,6]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.FS19.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS19.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS19.S[h,22]))
    }
    
    if(bialelicosNOimput.FS19.S[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.FS19.S[h,6] * 0.5 
      esp2 <- bialelicosNOimput.FS19.S[h,6] * 0.5 + bialelicosNOimput.FS19.S[h+1,6] * 0.5
      esp3 <- bialelicosNOimput.FS19.S[h+1,6] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.FS19.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS19.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS19.S[h,22]))
    }
    
    if(bialelicosNOimput.FS19.S[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.FS19.S[h,6] 
      esp3 <- bialelicosNOimput.FS19.S[h+1,6]
      
      obs1 <- as.double(paste(bialelicosNOimput.FS19.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS19.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS19.S[h,22]))
    }
    
    if(bialelicosNOimput.FS19.S[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.FS19.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS19.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS19.S[h,22]))
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
  
  FrecuenciasNOImputados.FS19.S[[a]] <- data.frame(x = c(obs1*sample.size$FS19,obs2*sample.size$FS19,obs3*sample.size$FS19), 
                                                   y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$FS19*obs1,sample.size$FS19*obs2,sample.size$FS19*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.FS19.S <- 1

while (g < nrow(bialelicosNOimput.FS19.S)) {
  
  if(is.na(FrecuenciasNOImputados.FS19.S[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.FS19.S[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.FS19.S[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.FS19.S[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.FS19.S[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.FS19.S[[l]][,2][3])){
    
    bialelicosNOimput.FS19.S[g,27] <- NA
    bialelicosNOimput.FS19.S[g,28] <- NA
    bialelicosNOimput.FS19.S[g+1,27] <- NA
    bialelicosNOimput.FS19.S[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.FS19.S <- FrecuenciasNOImputados.FS19.S[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.FS19.S[[l]][,1][1] == 0 && FrecuenciasNOImputados.FS19.S[[l]][,2][1] == 0){
    FrecuenciasNOImputados.FS19.S[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.FS19.S <- FrecuenciasNOImputados.FS19.S[[l]]
  }
  if(FrecuenciasNOImputados.FS19.S[[l]][,1][2] == 0 && FrecuenciasNOImputados.FS19.S[[l]][,2][2] == 0){
    FrecuenciasNOImputados.FS19.S[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.FS19.S <- FrecuenciasNOImputados.FS19.S[[l]]
  }
  if(FrecuenciasNOImputados.FS19.S[[l]][,1][3] == 0 && FrecuenciasNOImputados.FS19.S[[l]][,2][3] == 0){
    FrecuenciasNOImputados.FS19.S[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.FS19.S <- FrecuenciasNOImputados.FS19.S[[l]]
  }
  
  v.FS19.S <- chisq.test(x = m.FS19.S[,1],p = m.FS19.S[,2], rescale.p = TRUE)
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28   
  bialelicosNOimput.FS19.S[g,27] <- v.FS19.S$statistic
  bialelicosNOimput.FS19.S[g,28] <- v.FS19.S$p.value
  bialelicosNOimput.FS19.S[g+1,27] <- v.FS19.S$statistic
  bialelicosNOimput.FS19.S[g+1,28] <- v.FS19.S$p.value
  
  g <- g+2
  l <- l+1
}


#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.FS19.S, 
            "BialelicosNOImputados-FS19-S.txt", 
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
FrecuenciasImputados.FS19.S <- list()
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

while (h <= nrow(bialelicosimput.FS19.S)) {
  if(bialelicosimput.FS19.S[h,5] == 99 && bialelicosimput.FS19.S[h+1,5] == 0 ){  
    if ((bialelicosimput.FS19.S[h,20]+bialelicosimput.FS19.S[h,23]+bialelicosimput.FS19.S[h,25])==1){
      
      if(bialelicosimput.FS19.S[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.FS19.S[h+1,6] 
        esp2 <- bialelicosimput.FS19.S[h,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.FS19.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS19.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS19.S[h,25]))
      }
      
      if(bialelicosimput.FS19.S[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.FS19.S[h+1,6] + bialelicosimput.FS19.S[h,6] * 0.5 
        esp3 <- bialelicosimput.FS19.S[h,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.FS19.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS19.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS19.S[h,25]))
      }
      
      if(bialelicosimput.FS19.S[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.FS19.S[h+1,6] 
        esp3 <- bialelicosimput.FS19.S[h,6]
        
        obs1 <- as.double(paste(bialelicosimput.FS19.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS19.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS19.S[h,25]))
      }
      
      if(bialelicosimput.FS19.S[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.FS19.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS19.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS19.S[h,25]))
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
  
  
  if(bialelicosimput.FS19.S[h,5] == 0 && bialelicosimput.FS19.S[h+1,5] == 99 ){  
    if ((bialelicosimput.FS19.S[h,20]+bialelicosimput.FS19.S[h,23]+bialelicosimput.FS19.S[h,25])==1){
      
      if(bialelicosimput.FS19.S[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.FS19.S[h,6] 
        esp2 <- bialelicosimput.FS19.S[h+1,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.FS19.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS19.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS19.S[h,25]))
      }
      
      if(bialelicosimput.FS19.S[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.FS19.S[h,6] + bialelicosimput.FS19.S[h+1,6] * 0.5 
        esp3 <- bialelicosimput.FS19.S[h+1,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.FS19.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS19.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS19.S[h,25]))
      }
      
      if(bialelicosimput.FS19.S[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.FS19.S[h,6] 
        esp3 <- bialelicosimput.FS19.S[h+1,6]
        
        obs1 <- as.double(paste(bialelicosimput.FS19.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS19.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS19.S[h,25]))
      }
      
      if(bialelicosimput.FS19.S[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.FS19.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS19.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS19.S[h,25]))
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
  
  
  FrecuenciasImputados.FS19.S[[a]] <- data.frame(x = c(obs1*sample.size$FS19,obs2*sample.size$FS19,obs3*sample.size$FS19), 
                                                 y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$FS19*obs1,sample.size$FS19*obs2,sample.size$FS19*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.FS19.S <- 1

while (g < nrow(bialelicosimput.FS19.S)) {
  
  if(is.na(FrecuenciasImputados.FS19.S[[l]][,1][1]) 
     || is.na(FrecuenciasImputados.FS19.S[[l]][,1][2]) 
     || is.na(FrecuenciasImputados.FS19.S[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.FS19.S[[l]][,2][1]) 
     || is.na(FrecuenciasImputados.FS19.S[[l]][,2][2]) 
     || is.na(FrecuenciasImputados.FS19.S[[l]][,2][3])){
    
    bialelicosimput.FS19.S[g,27] <- NA
    bialelicosimput.FS19.S[g,28] <- NA
    bialelicosimput.FS19.S[g+1,27] <- NA
    bialelicosimput.FS19.S[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.FS19.S <- FrecuenciasImputados.FS19.S[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.FS19.S[[l]][,1][1] == 0 && FrecuenciasImputados.FS19.S[[l]][,2][1] == 0){
    FrecuenciasImputados.FS19.S[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.FS19.S <- FrecuenciasImputados.FS19.S[[l]]
  }
  if(FrecuenciasImputados.FS19.S[[l]][,1][2] == 0 && FrecuenciasImputados.FS19.S[[l]][,2][2] == 0){
    FrecuenciasImputados.FS19.S[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.FS19.S <- FrecuenciasImputados.FS19.S[[l]]
  }
  if(FrecuenciasImputados.FS19.S[[l]][,1][3] == 0 && FrecuenciasImputados.FS19.S[[l]][,2][3] == 0){
    FrecuenciasImputados.FS19.S[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.FS19.S <- FrecuenciasImputados.FS19.S[[l]]
  }
  
  v.FS19.S <- chisq.test(x = m.FS19.S[,1],p = m.FS19.S[,2], rescale.p = TRUE)
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28  
  
  bialelicosimput.FS19.S[g,27] <- v.FS19.S$statistic
  bialelicosimput.FS19.S[g,28] <- v.FS19.S$p.value
  bialelicosimput.FS19.S[g+1,27] <- v.FS19.S$statistic
  bialelicosimput.FS19.S[g+1,28] <- v.FS19.S$p.value
  
  
  
  g <- g+2
  l <- l+1
}

write.table(bialelicosimput.FS19.S , "BialelicosImputados-FS19-S.txt")
#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.FS19.S , 
            "BialelicosImputados-FS19-S.txt", 
            sep="\t",
            dec=",")

#####################################################################################

#########################
#                       # 
#          FS20.S       #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.FS20.S <- FS20.S %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.FS20.S <- data.frame(FrecuenciasAlelicas$variant_id,
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
                              FrecGeno2.FS20.S$X0.0,
                              FrecGeno2.FS20.S$X0.1,
                              FrecGeno2.FS20.S$X1.1,
                              FrecGeno2.FS20.S$X0.99,
                              FrecGeno2.FS20.S$X1.99,
                              FrecGeno2.FS20.S$X99.99,
                              mother.genotypes2$EFS.20)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.FS20.S <- subset(bialelic.FS20.S, bialelic.FS20.S$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.FS20.S <- subset(bialelic.FS20.S, bialelic.FS20.S$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################

#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2
FrecuenciasNOImputados.FS20.S <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.FS20.S)) {
  if ((bialelicosNOimput.FS20.S[h,20]+bialelicosNOimput.FS20.S[h,21]+bialelicosNOimput.FS20.S[h,22])==1){
    
    if(bialelicosNOimput.FS20.S[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.FS20.S[h,6] 
      esp2 <- bialelicosNOimput.FS20.S[h+1,6]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.FS20.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS20.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS20.S[h,22]))
    }
    
    if(bialelicosNOimput.FS20.S[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.FS20.S[h,6] * 0.5 
      esp2 <- bialelicosNOimput.FS20.S[h,6] * 0.5 + bialelicosNOimput.FS20.S[h+1,6] * 0.5
      esp3 <- bialelicosNOimput.FS20.S[h+1,6] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.FS20.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS20.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS20.S[h,22]))
    }
    
    if(bialelicosNOimput.FS20.S[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.FS20.S[h,6] 
      esp3 <- bialelicosNOimput.FS20.S[h+1,6]
      
      obs1 <- as.double(paste(bialelicosNOimput.FS20.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS20.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS20.S[h,22]))
    }
    
    if(bialelicosNOimput.FS20.S[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.FS20.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS20.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS20.S[h,22]))
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
  
  FrecuenciasNOImputados.FS20.S[[a]] <- data.frame(x = c(obs1*sample.size$FS20,obs2*sample.size$FS20,obs3*sample.size$FS20), 
                                                   y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$FS20*obs1,sample.size$FS20*obs2,sample.size$FS20*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.FS20.S <- 1

while (g < nrow(bialelicosNOimput.FS20.S)) {
  
  if(is.na(FrecuenciasNOImputados.FS20.S[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.FS20.S[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.FS20.S[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.FS20.S[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.FS20.S[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.FS20.S[[l]][,2][3])){
    
    bialelicosNOimput.FS20.S[g,27] <- NA
    bialelicosNOimput.FS20.S[g,28] <- NA
    bialelicosNOimput.FS20.S[g+1,27] <- NA
    bialelicosNOimput.FS20.S[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.FS20.S <- FrecuenciasNOImputados.FS20.S[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.FS20.S[[l]][,1][1] == 0 && FrecuenciasNOImputados.FS20.S[[l]][,2][1] == 0){
    FrecuenciasNOImputados.FS20.S[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.FS20.S <- FrecuenciasNOImputados.FS20.S[[l]]
  }
  if(FrecuenciasNOImputados.FS20.S[[l]][,1][2] == 0 && FrecuenciasNOImputados.FS20.S[[l]][,2][2] == 0){
    FrecuenciasNOImputados.FS20.S[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.FS20.S <- FrecuenciasNOImputados.FS20.S[[l]]
  }
  if(FrecuenciasNOImputados.FS20.S[[l]][,1][3] == 0 && FrecuenciasNOImputados.FS20.S[[l]][,2][3] == 0){
    FrecuenciasNOImputados.FS20.S[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.FS20.S <- FrecuenciasNOImputados.FS20.S[[l]]
  }
  
  v.FS20.S <- chisq.test(x = m.FS20.S[,1],p = m.FS20.S[,2], rescale.p = TRUE)
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28   
  bialelicosNOimput.FS20.S[g,27] <- v.FS20.S$statistic
  bialelicosNOimput.FS20.S[g,28] <- v.FS20.S$p.value
  bialelicosNOimput.FS20.S[g+1,27] <- v.FS20.S$statistic
  bialelicosNOimput.FS20.S[g+1,28] <- v.FS20.S$p.value
  
  g <- g+2
  l <- l+1
}


#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.FS20.S, 
            "BialelicosNOImputados-FS20-S.txt", 
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
FrecuenciasImputados.FS20.S <- list()
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

while (h <= nrow(bialelicosimput.FS20.S)) {
  if(bialelicosimput.FS20.S[h,5] == 99 && bialelicosimput.FS20.S[h+1,5] == 0 ){  
    if ((bialelicosimput.FS20.S[h,20]+bialelicosimput.FS20.S[h,23]+bialelicosimput.FS20.S[h,25])==1){
      
      if(bialelicosimput.FS20.S[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.FS20.S[h+1,6] 
        esp2 <- bialelicosimput.FS20.S[h,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.FS20.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS20.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS20.S[h,25]))
      }
      
      if(bialelicosimput.FS20.S[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.FS20.S[h+1,6] + bialelicosimput.FS20.S[h,6] * 0.5 
        esp3 <- bialelicosimput.FS20.S[h,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.FS20.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS20.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS20.S[h,25]))
      }
      
      if(bialelicosimput.FS20.S[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.FS20.S[h+1,6] 
        esp3 <- bialelicosimput.FS20.S[h,6]
        
        obs1 <- as.double(paste(bialelicosimput.FS20.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS20.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS20.S[h,25]))
      }
      
      if(bialelicosimput.FS20.S[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.FS20.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS20.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS20.S[h,25]))
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
  
  
  if(bialelicosimput.FS20.S[h,5] == 0 && bialelicosimput.FS20.S[h+1,5] == 99 ){  
    if ((bialelicosimput.FS20.S[h,20]+bialelicosimput.FS20.S[h,23]+bialelicosimput.FS20.S[h,25])==1){
      
      if(bialelicosimput.FS20.S[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.FS20.S[h,6] 
        esp2 <- bialelicosimput.FS20.S[h+1,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.FS20.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS20.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS20.S[h,25]))
      }
      
      if(bialelicosimput.FS20.S[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.FS20.S[h,6] + bialelicosimput.FS20.S[h+1,6] * 0.5 
        esp3 <- bialelicosimput.FS20.S[h+1,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.FS20.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS20.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS20.S[h,25]))
      }
      
      if(bialelicosimput.FS20.S[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.FS20.S[h,6] 
        esp3 <- bialelicosimput.FS20.S[h+1,6]
        
        obs1 <- as.double(paste(bialelicosimput.FS20.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS20.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS20.S[h,25]))
      }
      
      if(bialelicosimput.FS20.S[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.FS20.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS20.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS20.S[h,25]))
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
  
  
  FrecuenciasImputados.FS20.S[[a]] <- data.frame(x = c(obs1*sample.size$FS20,obs2*sample.size$FS20,obs3*sample.size$FS20), 
                                                 y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$FS20*obs1,sample.size$FS20*obs2,sample.size$FS20*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.FS20.S <- 1

while (g < nrow(bialelicosimput.FS20.S)) {
  
  if(is.na(FrecuenciasImputados.FS20.S[[l]][,1][1]) 
     || is.na(FrecuenciasImputados.FS20.S[[l]][,1][2]) 
     || is.na(FrecuenciasImputados.FS20.S[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.FS20.S[[l]][,2][1]) 
     || is.na(FrecuenciasImputados.FS20.S[[l]][,2][2]) 
     || is.na(FrecuenciasImputados.FS20.S[[l]][,2][3])){
    
    bialelicosimput.FS20.S[g,27] <- NA
    bialelicosimput.FS20.S[g,28] <- NA
    bialelicosimput.FS20.S[g+1,27] <- NA
    bialelicosimput.FS20.S[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.FS20.S <- FrecuenciasImputados.FS20.S[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.FS20.S[[l]][,1][1] == 0 && FrecuenciasImputados.FS20.S[[l]][,2][1] == 0){
    FrecuenciasImputados.FS20.S[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.FS20.S <- FrecuenciasImputados.FS20.S[[l]]
  }
  if(FrecuenciasImputados.FS20.S[[l]][,1][2] == 0 && FrecuenciasImputados.FS20.S[[l]][,2][2] == 0){
    FrecuenciasImputados.FS20.S[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.FS20.S <- FrecuenciasImputados.FS20.S[[l]]
  }
  if(FrecuenciasImputados.FS20.S[[l]][,1][3] == 0 && FrecuenciasImputados.FS20.S[[l]][,2][3] == 0){
    FrecuenciasImputados.FS20.S[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.FS20.S <- FrecuenciasImputados.FS20.S[[l]]
  }
  
  v.FS20.S <- chisq.test(x = m.FS20.S[,1],p = m.FS20.S[,2], rescale.p = TRUE)
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28  
  
  bialelicosimput.FS20.S[g,27] <- v.FS20.S$statistic
  bialelicosimput.FS20.S[g,28] <- v.FS20.S$p.value
  bialelicosimput.FS20.S[g+1,27] <- v.FS20.S$statistic
  bialelicosimput.FS20.S[g+1,28] <- v.FS20.S$p.value
  
  
  
  g <- g+2
  l <- l+1
}

write.table(bialelicosimput.FS20.S , "BialelicosImputados-FS20-S.txt")
#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.FS20.S , 
            "BialelicosImputados-FS20-S.txt", 
            sep="\t",
            dec=",")


#########################
#                       # 
#          FS22.S       #
#                       #
#########################

#We produce a new data.frame with duplicate rows from the genotypic frequencies files
#and the mother genotypes files using a function from dplyr
FrecGeno2.FS22.S <- FS22.S %>% slice(rep(1:n(),each =2))
mother.genotypes2 <- mother.genotypes %>% slice(rep(1:n(),each =2))

#We produce a new data.frame including columns from the allele frequencies file and
#from the mother genotypes file and observed genotypic frequencies for each variant(SNP):
# (0/0),(0/1),(1/1) in columns 20, 21 y 22, respectively. 
# (0/99), (1/99) y (99/99) in columns 23, 24 y 25, respectively
bialelic.FS22.S <- data.frame(FrecuenciasAlelicas$variant_id,
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
                              FrecGeno2.FS22.S$X0.0,
                              FrecGeno2.FS22.S$X0.1,
                              FrecGeno2.FS22.S$X1.1,
                              FrecGeno2.FS22.S$X0.99,
                              FrecGeno2.FS22.S$X1.99,
                              FrecGeno2.FS22.S$X99.99,
                              mother.genotypes2$EFS.22)

##We subset the bialelic data.frame to separate imputed from non-imputed SNP
bialelicosimput.FS22.S <- subset(bialelic.FS22.S, bialelic.FS22.S$FrecuenciasAlelicas.imputations == "Y")
bialelicosNOimput.FS22.S <- subset(bialelic.FS22.S, bialelic.FS22.S$FrecuenciasAlelicas.imputations == "N")

###################
#                 #
#   Non-Imputed   #
#                 #   
###################

#Expected frequency calculation. We build a list iteratively containing a 
#data.frame with the observed genotypic frequencies for 0/0, 0/1 y 1/1 in column 1 
#and with the expected frequencies for 0/0, 0/1 y 1/1 in column 2
FrecuenciasNOImputados.FS22.S <- list()
h <- 1
a <- 1
obs1 <- 0
obs2 <- 0
obs3 <- 0
esp1 <- 0
esp2 <- 0
esp3 <- 0
#We include a condition to be sure that genotypic frequencies of each variant sum up 1
while (h <= nrow(bialelicosNOimput.FS22.S)) {
  if ((bialelicosNOimput.FS22.S[h,20]+bialelicosNOimput.FS22.S[h,21]+bialelicosNOimput.FS22.S[h,22])==1){
    
    if(bialelicosNOimput.FS22.S[h,26] == "0/0"){
      
      esp1 <- bialelicosNOimput.FS22.S[h,6] 
      esp2 <- bialelicosNOimput.FS22.S[h+1,6]
      esp3 <- 0
      
      obs1 <- as.double(paste(bialelicosNOimput.FS22.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS22.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS22.S[h,22]))
    }
    
    if(bialelicosNOimput.FS22.S[h,26]  == "0/1"){
      
      esp1 <- bialelicosNOimput.FS22.S[h,6] * 0.5 
      esp2 <- bialelicosNOimput.FS22.S[h,6] * 0.5 + bialelicosNOimput.FS22.S[h+1,6] * 0.5
      esp3 <- bialelicosNOimput.FS22.S[h+1,6] * 0.5
      
      obs1 <- as.double(paste(bialelicosNOimput.FS22.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS22.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS22.S[h,22]))
    }
    
    if(bialelicosNOimput.FS22.S[h,26]  == "1/1"){  
      esp1 <- 0
      esp2 <- bialelicosNOimput.FS22.S[h,6] 
      esp3 <- bialelicosNOimput.FS22.S[h+1,6]
      
      obs1 <- as.double(paste(bialelicosNOimput.FS22.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS22.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS22.S[h,22]))
    }
    
    if(bialelicosNOimput.FS22.S[h,26] == "./."){
      esp1 <- NA
      esp2 <- NA
      esp3 <- NA 
      
      obs1 <- as.double(paste(bialelicosNOimput.FS22.S[h,20]))
      obs2 <- as.double(paste(bialelicosNOimput.FS22.S[h,21]))
      obs3 <- as.double(paste(bialelicosNOimput.FS22.S[h,22]))
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
  
  FrecuenciasNOImputados.FS22.S[[a]] <- data.frame(x = c(obs1*sample.size$FS22,obs2*sample.size$FS22,obs3*sample.size$FS22), 
                                                   y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$FS22*obs1,sample.size$FS22*obs2,sample.size$FS22*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.FS22.S <- 1

while (g < nrow(bialelicosNOimput.FS22.S)) {
  
  if(is.na(FrecuenciasNOImputados.FS22.S[[l]][,1][1]) 
     || is.na(FrecuenciasNOImputados.FS22.S[[l]][,1][2]) 
     || is.na(FrecuenciasNOImputados.FS22.S[[l]][,1][3]) 
     || is.na(FrecuenciasNOImputados.FS22.S[[l]][,2][1]) 
     || is.na(FrecuenciasNOImputados.FS22.S[[l]][,2][2]) 
     || is.na(FrecuenciasNOImputados.FS22.S[[l]][,2][3])){
    
    bialelicosNOimput.FS22.S[g,27] <- NA
    bialelicosNOimput.FS22.S[g,28] <- NA
    bialelicosNOimput.FS22.S[g+1,27] <- NA
    bialelicosNOimput.FS22.S[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.FS22.S <- FrecuenciasNOImputados.FS22.S[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasNOImputados.FS22.S[[l]][,1][1] == 0 && FrecuenciasNOImputados.FS22.S[[l]][,2][1] == 0){
    FrecuenciasNOImputados.FS22.S[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.FS22.S <- FrecuenciasNOImputados.FS22.S[[l]]
  }
  if(FrecuenciasNOImputados.FS22.S[[l]][,1][2] == 0 && FrecuenciasNOImputados.FS22.S[[l]][,2][2] == 0){
    FrecuenciasNOImputados.FS22.S[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.FS22.S <- FrecuenciasNOImputados.FS22.S[[l]]
  }
  if(FrecuenciasNOImputados.FS22.S[[l]][,1][3] == 0 && FrecuenciasNOImputados.FS22.S[[l]][,2][3] == 0){
    FrecuenciasNOImputados.FS22.S[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.FS22.S <- FrecuenciasNOImputados.FS22.S[[l]]
  }
  
  v.FS22.S <- chisq.test(x = m.FS22.S[,1],p = m.FS22.S[,2], rescale.p = TRUE)
  
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28   
  bialelicosNOimput.FS22.S[g,27] <- v.FS22.S$statistic
  bialelicosNOimput.FS22.S[g,28] <- v.FS22.S$p.value
  bialelicosNOimput.FS22.S[g+1,27] <- v.FS22.S$statistic
  bialelicosNOimput.FS22.S[g+1,28] <- v.FS22.S$p.value
  
  g <- g+2
  l <- l+1
}


#We output the results to the working directory in Tab separated text format
write.table(bialelicosNOimput.FS22.S, 
            "BialelicosNOImputados-FS22-S.txt", 
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
FrecuenciasImputados.FS22.S <- list()
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

while (h <= nrow(bialelicosimput.FS22.S)) {
  if(bialelicosimput.FS22.S[h,5] == 99 && bialelicosimput.FS22.S[h+1,5] == 0 ){  
    if ((bialelicosimput.FS22.S[h,20]+bialelicosimput.FS22.S[h,23]+bialelicosimput.FS22.S[h,25])==1){
      
      if(bialelicosimput.FS22.S[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.FS22.S[h+1,6] 
        esp2 <- bialelicosimput.FS22.S[h,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.FS22.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS22.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS22.S[h,25]))
      }
      
      if(bialelicosimput.FS22.S[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.FS22.S[h+1,6] + bialelicosimput.FS22.S[h,6] * 0.5 
        esp3 <- bialelicosimput.FS22.S[h,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.FS22.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS22.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS22.S[h,25]))
      }
      
      if(bialelicosimput.FS22.S[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.FS22.S[h+1,6] 
        esp3 <- bialelicosimput.FS22.S[h,6]
        
        obs1 <- as.double(paste(bialelicosimput.FS22.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS22.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS22.S[h,25]))
      }
      
      if(bialelicosimput.FS22.S[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.FS22.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS22.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS22.S[h,25]))
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
  
  
  if(bialelicosimput.FS22.S[h,5] == 0 && bialelicosimput.FS22.S[h+1,5] == 99 ){  
    if ((bialelicosimput.FS22.S[h,20]+bialelicosimput.FS22.S[h,23]+bialelicosimput.FS22.S[h,25])==1){
      
      if(bialelicosimput.FS22.S[h,26] == "0/0"){
        
        esp1 <- bialelicosimput.FS22.S[h,6] 
        esp2 <- bialelicosimput.FS22.S[h+1,6]
        esp3 <- 0
        
        obs1 <- as.double(paste(bialelicosimput.FS22.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS22.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS22.S[h,25]))
      }
      
      if(bialelicosimput.FS22.S[h,26] == "0/99"){
        
        esp1 <- 0 
        esp2 <- bialelicosimput.FS22.S[h,6] + bialelicosimput.FS22.S[h+1,6] * 0.5 
        esp3 <- bialelicosimput.FS22.S[h+1,6] * 0.5
        
        obs1 <- as.double(paste(bialelicosimput.FS22.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS22.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS22.S[h,25]))
      }
      
      if(bialelicosimput.FS22.S[h,26] == "99/99"){  
        esp1 <- 0
        esp2 <- bialelicosimput.FS22.S[h,6] 
        esp3 <- bialelicosimput.FS22.S[h+1,6]
        
        obs1 <- as.double(paste(bialelicosimput.FS22.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS22.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS22.S[h,25]))
      }
      
      if(bialelicosimput.FS22.S[h,26] == "./."){
        esp1 <- NA
        esp2 <- NA
        esp3 <- NA 
        
        obs1 <- as.double(paste(bialelicosimput.FS22.S[h,20]))
        obs2 <- as.double(paste(bialelicosimput.FS22.S[h,23]))
        obs3 <- as.double(paste(bialelicosimput.FS22.S[h,25]))
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
  
  
  FrecuenciasImputados.FS22.S[[a]] <- data.frame(x = c(obs1*sample.size$FS22,obs2*sample.size$FS22,obs3*sample.size$FS22), 
                                                 y = c(esp1,esp2,esp3))
  
  
  h <- h+2
  a <- a+1
  
}

#Iterative calculation of the Chi2 statistics and the p-value to detect
#SNP deviating from the expectations
#The observed genotypic frequencies are multiplied by the number of individuals
#because they must be in absolute terms x=c(sample.size$FS22*obs1,sample.size$FS22*obs2,sample.size$FS22*obs3,)
#The expected frequencies are taken as an expected probability p=c(esp1,esp2,esp3)
#With a short number of observations, Warnings will arise (neglectable)

g <- 1
l <- 1
m.FS22.S <- 1

while (g < nrow(bialelicosimput.FS22.S)) {
  
  if(is.na(FrecuenciasImputados.FS22.S[[l]][,1][1]) 
     || is.na(FrecuenciasImputados.FS22.S[[l]][,1][2]) 
     || is.na(FrecuenciasImputados.FS22.S[[l]][,1][3]) 
     || is.na(FrecuenciasImputados.FS22.S[[l]][,2][1]) 
     || is.na(FrecuenciasImputados.FS22.S[[l]][,2][2]) 
     || is.na(FrecuenciasImputados.FS22.S[[l]][,2][3])){
    
    bialelicosimput.FS22.S[g,27] <- NA
    bialelicosimput.FS22.S[g,28] <- NA
    bialelicosimput.FS22.S[g+1,27] <- NA
    bialelicosimput.FS22.S[g+1,28] <- NA
    
    g <- g+2
    l <- l+1
    
    next()
    
  }    
  
  m.FS22.S <- FrecuenciasImputados.FS22.S[[l]]  
  #We add a small number 10⁻39 to the elements of the vector of expected frequencies that are 0
  # to prevent the malfunction of the chisq.test() function
  #This step is optional when using rescale.p = TRUE
  
  if(FrecuenciasImputados.FS22.S[[l]][,1][1] == 0 && FrecuenciasImputados.FS22.S[[l]][,2][1] == 0){
    FrecuenciasImputados.FS22.S[[l]][,2][1] <- 0.000000000000000000000000000000000000001
    m.FS22.S <- FrecuenciasImputados.FS22.S[[l]]
  }
  if(FrecuenciasImputados.FS22.S[[l]][,1][2] == 0 && FrecuenciasImputados.FS22.S[[l]][,2][2] == 0){
    FrecuenciasImputados.FS22.S[[l]][,2][2] <- 0.000000000000000000000000000000000000001
    m.FS22.S <- FrecuenciasImputados.FS22.S[[l]]
  }
  if(FrecuenciasImputados.FS22.S[[l]][,1][3] == 0 && FrecuenciasImputados.FS22.S[[l]][,2][3] == 0){
    FrecuenciasImputados.FS22.S[[l]][,2][3] <- 0.000000000000000000000000000000000000001
    m.FS22.S <- FrecuenciasImputados.FS22.S[[l]]
  }
  
  v.FS22.S <- chisq.test(x = m.FS22.S[,1],p = m.FS22.S[,2], rescale.p = TRUE)
  #We add the scored Chi2 parameter and the corresponding p-value to the data.frame in columns 27 and 28  
  
  bialelicosimput.FS22.S[g,27] <- v.FS22.S$statistic
  bialelicosimput.FS22.S[g,28] <- v.FS22.S$p.value
  bialelicosimput.FS22.S[g+1,27] <- v.FS22.S$statistic
  bialelicosimput.FS22.S[g+1,28] <- v.FS22.S$p.value
  
  
  
  g <- g+2
  l <- l+1
}

write.table(bialelicosimput.FS22.S , "BialelicosImputados-FS22-S.txt")
#We output the results to the working directory in Tab separated text format
write.table(bialelicosimput.FS22.S , 
            "BialelicosImputados-FS22-S.txt", 
            sep="\t",
            dec=",")

#####################################################################################



##################
#                #
#       END      #
#                #
##################




