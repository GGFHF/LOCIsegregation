###########################################
#                                         #
#       5-Crossing with RNAseq data       #
#                                         #
###########################################
#This software has been developed by:

#  Dpto. Sistemas y Recursos Naturales
#  ETSI Montes, Forestal y del Medio Natural
#  Universidad Politecnica de Madrid
# https://github.com/ggfhf/


#We run the config file
#source('0-Config.R')

#We execute the necessary scripts (hosted in the same working directory)
#The WD in the scripts must be the same that in this script
#We set the working directory to run the scripts. Be sure all the scripts are within this directory
setwd(WD.input.0)
#source('1-Chi2-test-adults.R')
#source('2-Chi2-test-suber-ilex-progenies.R')
#source('3-Chi2-test-hybrid-progenies.R')
source('4-Significant-SNP-Venn-Diagrams.R')

#We set the working directory. Be sure all the input files are within this directory
setwd(WD.input)

#Once the scripts are executed, we upload the DEGs table
#Important columns are cluster and Ortholog with arabidopsis to use them in ShinyGO
RNASEQ.DEGS <- read.table("RNASEQ-DEGSa.txt", 
                          header=TRUE,
                          sep="\t",
                          dec=",")
#Full correspondence table between Quercus suber LOC and Arabidopsis ATG
full.correspondence.LOC.ATG <- read.table("CorrespondenceLOC-ATG.txt", 
                          header=TRUE,
                          sep="\t",
                          dec=",")

#We set the working directory to produce the output.
setwd(WD.output.2)

#We'll keep only those significant LOC significativos that were differentially expressed in the RNASEQ study
inner.sig.res.bialelicosNOimput.EFS <- inner_join(sig.res.bialelicosNOimput.EFS,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.EFS ,
            "inner.sig.res.bialelicosNOimput-EFS.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.AL <- inner_join(sig.res.bialelicosNOimput.AL,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.AL ,
            "inner.sig.res.bialelicosNOimput-AL.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.EN <- inner_join(sig.res.bialelicosNOimput.EN,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.EN ,
            "inner.sig.res.bialelicosNOimput-EN.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.A05 <- inner_join(sig.res.bialelicosNOimput.A05,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.A05 ,
            "inner.sig.res.bialelicosNOimput-A05.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.A07 <- inner_join(sig.res.bialelicosNOimput.A07,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.A07 ,
            "inner.sig.res.bialelicosNOimput-A07.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.A09 <- inner_join(sig.res.bialelicosNOimput.A09,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.A09 ,
            "inner.sig.res.bialelicosNOimput-A09.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.A10 <- inner_join(sig.res.bialelicosNOimput.A10,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.A10 ,
            "inner.sig.res.bialelicosNOimput-A10.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.E28 <- inner_join(sig.res.bialelicosNOimput.E28,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.E28 ,
            "inner.sig.res.bialelicosNOimput-E28.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.E31 <- inner_join(sig.res.bialelicosNOimput.E31,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.E31 ,
            "inner.sig.res.bialelicosNOimput-E31.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.E41 <- inner_join(sig.res.bialelicosNOimput.E41,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.E41 ,
            "inner.sig.res.bialelicosNOimput-E41.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.E96 <- inner_join(sig.res.bialelicosNOimput.E96,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.E96 ,
            "inner.sig.res.bialelicosNOimput-E96.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.FS16.S <- inner_join(sig.res.bialelicosNOimput.FS16.S,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.FS16.S ,
            "inner.sig.res.bialelicosNOimput-FS16.-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.FS19.S <- inner_join(sig.res.bialelicosNOimput.FS19.S,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.FS19.S ,
            "inner.sig.res.bialelicosNOimput-FS19-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.FS20.S <- inner_join(sig.res.bialelicosNOimput.FS20.S,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.FS20.S ,
            "inner.sig.res.bialelicosNOimput-FS20-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosNOimput.FS22.S <- inner_join(sig.res.bialelicosNOimput.FS22.S,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.FS22.S ,
            "inner.sig.res.bialelicosNOimput-FS22-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.EFS <- inner_join(sig.res.bialelicosimput.EFS,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.EFS ,
            "inner.sig.res.bialelicosimput-EFS.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.AL <- inner_join(sig.res.bialelicosimput.AL,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.AL ,
            "inner.sig.res.bialelicosimput-AL.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.EN <- inner_join(sig.res.bialelicosimput.EN,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.EN ,
            "inner.sig.res.bialelicosimput-EN.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.A05 <- inner_join(sig.res.bialelicosimput.A05,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.A05 ,
            "inner.sig.res.bialelicosimput-A05.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.A07 <- inner_join(sig.res.bialelicosimput.A07,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.A07 ,
            "inner.sig.res.bialelicosimput-A07.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.A09 <- inner_join(sig.res.bialelicosimput.A09,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.A09 ,
            "inner.sig.res.bialelicosimput-A09.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.A10 <- inner_join(sig.res.bialelicosimput.A10,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.A10 ,
            "inner.sig.res.bialelicosimput-A10.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.E28 <- inner_join(sig.res.bialelicosimput.E28,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.E28 ,
            "inner.sig.res.bialelicosimput-E28.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.E31 <- inner_join(sig.res.bialelicosimput.E31,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.E31 ,
            "inner.sig.res.bialelicosimput-E31.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.E41 <- inner_join(sig.res.bialelicosimput.E41,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.E41 ,
            "inner.sig.res.bialelicosimput-E41.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.E96 <- inner_join(sig.res.bialelicosimput.E96,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.E96 ,
            "inner.sig.res.bialelicosimput-E96.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.FS16.S <- inner_join(sig.res.bialelicosimput.FS16.S,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.FS16.S ,
            "inner.sig.res.bialelicosimput-FS16.-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.FS19.S <- inner_join(sig.res.bialelicosimput.FS19.S,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.FS19.S ,
            "inner.sig.res.bialelicosimput-FS19-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.FS20.S <- inner_join(sig.res.bialelicosimput.FS20.S,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.FS20.S ,
            "inner.sig.res.bialelicosimput-FS20-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
inner.sig.res.bialelicosimput.FS22.S <- inner_join(sig.res.bialelicosimput.FS22.S,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.FS22.S ,
            "inner.sig.res.bialelicosimput-FS22-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)


###########################################################################################################
#Now we repeat this for the intersected SNP for the adult hybrids

inner.sig.res.bialelicosNOimput.EFS.SNP <- inner_join(sig.res.bialelicosNOimput.EFS,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.EFS.SNP ,
            "inner-sig-res-bialelicosNOimput-EFS-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

inner.sig.res.bialelicosimput.EFS.SNP <- inner_join(sig.res.bialelicosimput.EFS,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.EFS.SNP ,
            "inner-sig-res-bialelicosimput-EFS-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)


#Now we repeat this for the intersected SNP within the families of each species

##hybrid
inner.sig.res.bialelicosNOimput.FS16.FS19.FS20.FS22.SNP <- inner_join(sig.res.bialelicosNOimput.FS16.FS19.FS20.FS22.SNP,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.FS16.FS19.FS20.FS22.SNP ,
            "inner-sig-res-bialelicosNOimput-FS16-FS19-FS20-FS22-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

##suber
inner.sig.res.bialelicosNOimput.A05.A07.A09.A10.SNP <- inner_join(sig.res.bialelicosNOimput.A05.A07.A09.A10.SNP,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.A05.A07.A09.A10.SNP ,
            "inner-sig-res-bialelicosNOimput-A05-A07-A09-A10-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

##ilex
inner.sig.res.bialelicosNOimput.E28.E31.E41.E96.SNP <- inner_join(sig.res.bialelicosNOimput.E28.E31.E41.E96.SNP,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosNOimput.E28.E31.E41.E96.SNP ,
            "inner-sig-res-bialelicosNOimput-E28-E31-E41-E96-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

##hybrid
inner.sig.res.bialelicosimput.FS16.FS19.FS20.FS22.SNP <- inner_join(sig.res.bialelicosimput.FS16.FS19.FS20.FS22.SNP,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.FS16.FS19.FS20.FS22.SNP ,
            "inner-sig-res-bialelicosimput-FS16-FS19-FS20-FS22-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

##suber
inner.sig.res.bialelicosimput.A05.A07.A09.A10.SNP <- inner_join(sig.res.bialelicosimput.A05.A07.A09.A10.SNP,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.A05.A07.A09.A10.SNP ,
            "inner-sig-res-bialelicosimput-A05-A07-A09-A10-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

##ilex
inner.sig.res.bialelicosimput.E28.E31.E41.E96.SNP <- inner_join(sig.res.bialelicosimput.E28.E31.E41.E96.SNP,RNASEQ.DEGS, by="gene" ) 
write.table(inner.sig.res.bialelicosimput.E28.E31.E41.E96.SNP ,
            "inner-sig-res-bialelicosimput-E28-E31-E41-E96-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)


#We create unique data.frames combining imputed and non-imputed data with rbind
inner.sig.res.bialelicos.FS16.FS19.FS20.FS22.SNP <- rbind(inner.sig.res.bialelicosNOimput.FS16.FS19.FS20.FS22.SNP,
                                                                 inner.sig.res.bialelicosimput.FS16.FS19.FS20.FS22.SNP)

inner.sig.res.bialelicos.A05.A07.A09.A10.SNP <- rbind(inner.sig.res.bialelicosNOimput.A05.A07.A09.A10.SNP,
                                                          inner.sig.res.bialelicosimput.A05.A07.A09.A10.SNP)

inner.sig.res.bialelicos.E28.E31.E41.E96.SNP <- rbind(inner.sig.res.bialelicosNOimput.E28.E31.E41.E96.SNP,
                                                      inner.sig.res.bialelicosimput.E28.E31.E41.E96.SNP)

###############################################################################################



#To obtain the Arabidopsis ID for shinygo, we do the same, but using the full correspondence table

##hybrid adults
full.inner.sig.res.bialelicosNOimput.EFS.SNP <- inner_join(sig.res.bialelicosNOimput.EFS,full.correspondence.LOC.ATG, by="gene" ) 
write.table(full.inner.sig.res.bialelicosNOimput.EFS.SNP ,
            "full.inner-sig-res-bialelicosNOimput-EFS-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)


##hybrid families
full.inner.sig.res.bialelicosNOimput.FS16.FS19.FS20.FS22.SNP <- inner_join(sig.res.bialelicosNOimput.FS16.FS19.FS20.FS22.SNP,full.correspondence.LOC.ATG, by="gene" ) 
write.table(full.inner.sig.res.bialelicosNOimput.FS16.FS19.FS20.FS22.SNP ,
            "full.inner-sig-res-bialelicosNOimput-FS16-FS19-FS20-FS22-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

##suber families
full.inner.sig.res.bialelicosNOimput.A05.A07.A09.A10.SNP <- inner_join(sig.res.bialelicosNOimput.A05.A07.A09.A10.SNP,full.correspondence.LOC.ATG, by="gene" ) 
write.table(full.inner.sig.res.bialelicosNOimput.A05.A07.A09.A10.SNP ,
            "full.inner-sig-res-bialelicosNOimput-A05-A07-A09-A10-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

##ilex families
full.inner.sig.res.bialelicosNOimput.E28.E31.E41.E96.SNP <- inner_join(sig.res.bialelicosNOimput.E28.E31.E41.E96.SNP,full.correspondence.LOC.ATG, by="gene" ) 
write.table(full.inner.sig.res.bialelicosNOimput.E28.E31.E41.E96.SNP ,
            "full.inner-sig-res-bialelicosNOimput-E28-E31-E41-E96-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

#hybrid adults
full.inner.sig.res.bialelicosimput.EFS.SNP <- inner_join(sig.res.bialelicosimput.EFS,full.correspondence.LOC.ATG, by="gene" ) 
write.table(full.inner.sig.res.bialelicosimput.EFS.SNP ,
            "full.inner-sig-res-bialelicosimput-EFS-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)


##hybrid families
full.inner.sig.res.bialelicosimput.FS16.FS19.FS20.FS22.SNP <- inner_join(sig.res.bialelicosimput.FS16.FS19.FS20.FS22.SNP,full.correspondence.LOC.ATG, by="gene" ) 
write.table(full.inner.sig.res.bialelicosimput.FS16.FS19.FS20.FS22.SNP ,
            "full.inner-sig-res-bialelicosimput-FS16-FS19-FS20-FS22-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

##suber families
full.inner.sig.res.bialelicosimput.A05.A07.A09.A10.SNP <- inner_join(sig.res.bialelicosimput.A05.A07.A09.A10.SNP,full.correspondence.LOC.ATG, by="gene" ) 
write.table(full.inner.sig.res.bialelicosimput.A05.A07.A09.A10.SNP ,
            "full.inner-sig-res-bialelicosimput-A05-A07-A09-A10-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

##ilex families
full.inner.sig.res.bialelicosimput.E28.E31.E41.E96.SNP <- inner_join(sig.res.bialelicosimput.E28.E31.E41.E96.SNP,full.correspondence.LOC.ATG, by="gene" ) 
write.table(full.inner.sig.res.bialelicosimput.E28.E31.E41.E96.SNP ,
            "full.inner-sig-res-bialelicosimput-E28-E31-E41-E96-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)


#We create unique data.frames combining imputed and non-imputed data with rbind

full.inner.sig.res.bialelicos.EFS.SNP <- rbind(full.inner.sig.res.bialelicosNOimput.EFS.SNP,
                                              full.inner.sig.res.bialelicosimput.EFS.SNP)


full.inner.sig.res.bialelicos.FS16.FS19.FS20.FS22.SNP <- rbind(full.inner.sig.res.bialelicosNOimput.FS16.FS19.FS20.FS22.SNP,
                                                               full.inner.sig.res.bialelicosimput.FS16.FS19.FS20.FS22.SNP)

full.inner.sig.res.bialelicos.A05.A07.A09.A10.SNP <- rbind(full.inner.sig.res.bialelicosNOimput.A05.A07.A09.A10.SNP,
                                                           full.inner.sig.res.bialelicosimput.A05.A07.A09.A10.SNP)

full.inner.sig.res.bialelicos.E28.E31.E41.E96.SNP <- rbind(full.inner.sig.res.bialelicosNOimput.E28.E31.E41.E96.SNP,
                                                           full.inner.sig.res.bialelicosimput.E28.E31.E41.E96.SNP)

#We extract a vector with unique ATG identifiers to run ShinyGO
hybrid.adults.ATG <- sort(unique(as.vector(full.inner.sig.res.bialelicos.EFS.SNP$Arabidopsis.ID)))
length(hybrid.adults.ATG)
write.table(hybrid.adults.ATG, 
            "hybrid-adults-ATG.txt",
            row.names =FALSE)

hybrid.families.ATG <- sort(unique(as.vector(full.inner.sig.res.bialelicos.FS16.FS19.FS20.FS22.SNP$Arabidopsis.ID)))
length(hybrid.families.ATG)
write.table(hybrid.families.ATG, 
            "hybrid-families-ATG.txt",
            row.names =FALSE)
suber.families.ATG <- sort(unique(as.vector(full.inner.sig.res.bialelicos.A05.A07.A09.A10.SNP$Arabidopsis.ID)))
length(suber.families.ATG)
write.table(suber.families.ATG, 
            "suber-families-ATG.txt",
            row.names =FALSE)
ilex.families.ATG <- sort(unique(as.vector(full.inner.sig.res.bialelicos.E28.E31.E41.E96.SNP$Arabidopsis.ID)))
length(ilex.families.ATG)
write.table(ilex.families.ATG, 
            "ilex-families-ATG.txt",
            row.names =FALSE)

#We create a Venn diagram to classify all the Arabidopsis genes found
x <- list(HY.AD = hybrid.adults.ATG, 
          HY.PG = hybrid.families.ATG, 
          SU.PG = suber.families.ATG,
          IL.PG = ilex.families.ATG)

pdf("ArabidopsisIDs-Common.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
 )
)
dev.off()

