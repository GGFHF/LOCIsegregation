##########################################################
#                                                        #
#       4-Subset significant SNP and Venn diagrams       #
#                                                        #
##########################################################
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
source('1-Chi2-test-adults.R')
setwd(WD.input.0)
source('2-Chi2-test-suber-ilex-progenies.R')
setwd(WD.input.0)
source('3-Chi2-test-hybrid-progenies.R')


#We set the working directory to produce the output (don't need an WD with input files).
setwd(WD.output.2)

#Once the scripts are executed, we construct global data.frames
#to include only the relevant results in each Adult population/family
#We use Distinct(), which is a function of dplyr that removes duplicated rows

#Construimos ahora unos data.frame globales con los resultados 
#de bialelicosImputados y bialelicosNOimputados para cada conjunto de datos
#Distinct() es una función del paquete dplyr que elimina las filas duplicadas
#We also change the column names to ease indexing
colnames.df<- c("variant_id",
                "imputation",
                "genomic_zone",
                "gene",
                "description",
                "robur_chr",
                "Chi",
                "pvalue")

res.bialelicosNOimput.EFS <- distinct(data.frame(bialelicosNOimput.EFS$FrecuenciasAlelicas.variant_id,
                                        bialelicosNOimput.EFS$FrecuenciasAlelicas.imputations,
                                        bialelicosNOimput.EFS$FrecuenciasAlelicas.genomic_zone,
                                        bialelicosNOimput.EFS$FrecuenciasAlelicas.gene.fragment,
                                        bialelicosNOimput.EFS$FrecuenciasAlelicas.description,
                                        bialelicosNOimput.EFS$FrecuenciasAlelicas.chromosome_id,
                                        bialelicosNOimput.EFS$V26,
                                        bialelicosNOimput.EFS$V27))

res.bialelicosNOimput.EN <- distinct(data.frame(bialelicosNOimput.EN$FrecuenciasAlelicas.variant_id,
                                                bialelicosNOimput.EN$FrecuenciasAlelicas.imputations,
                                                bialelicosNOimput.EN$FrecuenciasAlelicas.genomic_zone,
                                                bialelicosNOimput.EN$FrecuenciasAlelicas.gene.fragment,
                                                bialelicosNOimput.EN$FrecuenciasAlelicas.description,
                                                bialelicosNOimput.EN$FrecuenciasAlelicas.chromosome_id,
                                                bialelicosNOimput.EN$V26,
                                                bialelicosNOimput.EN$V27))

res.bialelicosNOimput.AL <- distinct(data.frame(bialelicosNOimput.AL$FrecuenciasAlelicas.variant_id,
                                                bialelicosNOimput.AL$FrecuenciasAlelicas.imputations,
                                                bialelicosNOimput.AL$FrecuenciasAlelicas.genomic_zone,
                                                bialelicosNOimput.AL$FrecuenciasAlelicas.gene.fragment,
                                                bialelicosNOimput.AL$FrecuenciasAlelicas.description,
                                                bialelicosNOimput.AL$FrecuenciasAlelicas.chromosome_id,
                                                bialelicosNOimput.AL$V26,
                                                bialelicosNOimput.AL$V27))

res.bialelicosimput.EFS <- distinct(data.frame(bialelicosimput.EFS$FrecuenciasAlelicas.variant_id,
                                               bialelicosimput.EFS$FrecuenciasAlelicas.imputations,
                                               bialelicosimput.EFS$FrecuenciasAlelicas.genomic_zone,
                                               bialelicosimput.EFS$FrecuenciasAlelicas.gene.fragment,
                                               bialelicosimput.EFS$FrecuenciasAlelicas.description,
                                               bialelicosimput.EFS$FrecuenciasAlelicas.chromosome_id,
                                               bialelicosimput.EFS$V26,
                                               bialelicosimput.EFS$V27))

res.bialelicosimput.EN <- distinct(data.frame(bialelicosimput.EN$FrecuenciasAlelicas.variant_id,
                                              bialelicosimput.EN$FrecuenciasAlelicas.imputations,
                                              bialelicosimput.EN$FrecuenciasAlelicas.genomic_zone,
                                              bialelicosimput.EN$FrecuenciasAlelicas.gene.fragment,
                                              bialelicosimput.EN$FrecuenciasAlelicas.description,
                                              bialelicosimput.EN$FrecuenciasAlelicas.chromosome_id,
                                              bialelicosimput.EN$V26,
                                              bialelicosimput.EN$V27))

res.bialelicosimput.AL <- distinct(data.frame(bialelicosimput.AL$FrecuenciasAlelicas.variant_id,
                                              bialelicosimput.AL$FrecuenciasAlelicas.imputations,
                                              bialelicosimput.AL$FrecuenciasAlelicas.genomic_zone,
                                              bialelicosimput.AL$FrecuenciasAlelicas.gene.fragment,
                                              bialelicosimput.AL$FrecuenciasAlelicas.description,
                                              bialelicosimput.AL$FrecuenciasAlelicas.chromosome_id,
                                              bialelicosimput.AL$V26,
                                              bialelicosimput.AL$V27))

res.bialelicosNOimput.A05 <- distinct(data.frame(bialelicosNOimput.A05$FrecuenciasAlelicas.variant_id,
                                                 bialelicosNOimput.A05$FrecuenciasAlelicas.imputations,
                                                 bialelicosNOimput.A05$FrecuenciasAlelicas.genomic_zone,
                                                 bialelicosNOimput.A05$FrecuenciasAlelicas.gene.fragment,
                                                 bialelicosNOimput.A05$FrecuenciasAlelicas.description,
                                                 bialelicosNOimput.A05$FrecuenciasAlelicas.chromosome_id,
                                                 bialelicosNOimput.A05$V27,
                                                 bialelicosNOimput.A05$V28))

res.bialelicosNOimput.A07 <- distinct(data.frame(bialelicosNOimput.A07$FrecuenciasAlelicas.variant_id,
                                                 bialelicosNOimput.A07$FrecuenciasAlelicas.imputations,
                                                 bialelicosNOimput.A07$FrecuenciasAlelicas.genomic_zone,
                                                 bialelicosNOimput.A07$FrecuenciasAlelicas.gene.fragment,
                                                 bialelicosNOimput.A07$FrecuenciasAlelicas.description,
                                                 bialelicosNOimput.A07$FrecuenciasAlelicas.chromosome_id,
                                                 bialelicosNOimput.A07$V27,
                                                 bialelicosNOimput.A07$V28))
res.bialelicosNOimput.A09 <- distinct(data.frame(bialelicosNOimput.A09$FrecuenciasAlelicas.variant_id,
                                                 bialelicosNOimput.A09$FrecuenciasAlelicas.imputations,
                                                 bialelicosNOimput.A09$FrecuenciasAlelicas.genomic_zone,
                                                 bialelicosNOimput.A09$FrecuenciasAlelicas.gene.fragment,
                                                 bialelicosNOimput.A09$FrecuenciasAlelicas.description,
                                                 bialelicosNOimput.A09$FrecuenciasAlelicas.chromosome_id,
                                                 bialelicosNOimput.A09$V27,
                                                 bialelicosNOimput.A09$V28))

res.bialelicosNOimput.A10 <- distinct(data.frame(bialelicosNOimput.A10$FrecuenciasAlelicas.variant_id,
                                                 bialelicosNOimput.A10$FrecuenciasAlelicas.imputations,
                                                 bialelicosNOimput.A10$FrecuenciasAlelicas.genomic_zone,
                                                 bialelicosNOimput.A10$FrecuenciasAlelicas.gene.fragment,
                                                 bialelicosNOimput.A10$FrecuenciasAlelicas.description,
                                                 bialelicosNOimput.A10$FrecuenciasAlelicas.chromosome_id,
                                                 bialelicosNOimput.A10$V27,
                                                 bialelicosNOimput.A10$V28))

res.bialelicosimput.A05 <- distinct(data.frame(bialelicosimput.A05$FrecuenciasAlelicas.variant_id,
                                               bialelicosimput.A05$FrecuenciasAlelicas.imputations,
                                               bialelicosimput.A05$FrecuenciasAlelicas.genomic_zone,
                                               bialelicosimput.A05$FrecuenciasAlelicas.gene.fragment,
                                               bialelicosimput.A05$FrecuenciasAlelicas.description,
                                               bialelicosimput.A05$FrecuenciasAlelicas.chromosome_id,
                                               bialelicosimput.A05$V27,
                                               bialelicosimput.A05$V28))

res.bialelicosimput.A07 <- distinct(data.frame(bialelicosimput.A07$FrecuenciasAlelicas.variant_id,
                                               bialelicosimput.A07$FrecuenciasAlelicas.imputations,
                                               bialelicosimput.A07$FrecuenciasAlelicas.genomic_zone,
                                               bialelicosimput.A07$FrecuenciasAlelicas.gene.fragment,
                                               bialelicosimput.A07$FrecuenciasAlelicas.description,
                                               bialelicosimput.A07$FrecuenciasAlelicas.chromosome_id,
                                               bialelicosimput.A07$V27,
                                               bialelicosimput.A07$V28))

res.bialelicosimput.A09 <- distinct(data.frame(bialelicosimput.A09$FrecuenciasAlelicas.variant_id,
                                               bialelicosimput.A09$FrecuenciasAlelicas.imputations,
                                               bialelicosimput.A09$FrecuenciasAlelicas.genomic_zone,
                                               bialelicosimput.A09$FrecuenciasAlelicas.gene.fragment,
                                               bialelicosimput.A09$FrecuenciasAlelicas.description,
                                               bialelicosimput.A09$FrecuenciasAlelicas.chromosome_id,
                                               bialelicosimput.A09$V27,
                                               bialelicosimput.A09$V28))

res.bialelicosimput.A10 <- distinct(data.frame(bialelicosimput.A10$FrecuenciasAlelicas.variant_id,
                                               bialelicosimput.A10$FrecuenciasAlelicas.imputations,
                                               bialelicosimput.A10$FrecuenciasAlelicas.genomic_zone,
                                               bialelicosimput.A10$FrecuenciasAlelicas.gene.fragment,
                                               bialelicosimput.A10$FrecuenciasAlelicas.description,
                                               bialelicosimput.A10$FrecuenciasAlelicas.chromosome_id,
                                               bialelicosimput.A10$V27,
                                               bialelicosimput.A10$V28))

res.bialelicosNOimput.E28 <- distinct(data.frame(bialelicosNOimput.E28$FrecuenciasAlelicas.variant_id,
                                                 bialelicosNOimput.E28$FrecuenciasAlelicas.imputations,
                                                 bialelicosNOimput.E28$FrecuenciasAlelicas.genomic_zone,
                                                 bialelicosNOimput.E28$FrecuenciasAlelicas.gene.fragment,
                                                 bialelicosNOimput.E28$FrecuenciasAlelicas.description,
                                                 bialelicosNOimput.E28$FrecuenciasAlelicas.chromosome_id,
                                                 bialelicosNOimput.E28$V27,
                                                 bialelicosNOimput.E28$V28))

res.bialelicosNOimput.E31 <- distinct(data.frame(bialelicosNOimput.E31$FrecuenciasAlelicas.variant_id,
                                                 bialelicosNOimput.E31$FrecuenciasAlelicas.imputations,
                                                 bialelicosNOimput.E31$FrecuenciasAlelicas.genomic_zone,
                                                 bialelicosNOimput.E31$FrecuenciasAlelicas.gene.fragment,
                                                 bialelicosNOimput.E31$FrecuenciasAlelicas.description,
                                                 bialelicosNOimput.E31$FrecuenciasAlelicas.chromosome_id,
                                                 bialelicosNOimput.E31$V27,
                                                 bialelicosNOimput.E31$V28))

res.bialelicosNOimput.E41 <- distinct(data.frame(bialelicosNOimput.E41$FrecuenciasAlelicas.variant_id,
                                                 bialelicosNOimput.E41$FrecuenciasAlelicas.imputations,
                                                 bialelicosNOimput.E41$FrecuenciasAlelicas.genomic_zone,
                                                 bialelicosNOimput.E41$FrecuenciasAlelicas.gene.fragment,
                                                 bialelicosNOimput.E41$FrecuenciasAlelicas.description,
                                                 bialelicosNOimput.E41$FrecuenciasAlelicas.chromosome_id,
                                                 bialelicosNOimput.E41$V27,
                                                 bialelicosNOimput.E41$V28))

res.bialelicosNOimput.E96 <- distinct(data.frame(bialelicosNOimput.E96$FrecuenciasAlelicas.variant_id,
                                                 bialelicosNOimput.E96$FrecuenciasAlelicas.imputations,
                                                 bialelicosNOimput.E96$FrecuenciasAlelicas.genomic_zone,
                                                 bialelicosNOimput.E96$FrecuenciasAlelicas.gene.fragment,
                                                 bialelicosNOimput.E96$FrecuenciasAlelicas.description,
                                                 bialelicosNOimput.E96$FrecuenciasAlelicas.chromosome_id,
                                                 bialelicosNOimput.E96$V27,
                                                 bialelicosNOimput.E96$V28))

res.bialelicosimput.E28 <- distinct(data.frame(bialelicosimput.E28$FrecuenciasAlelicas.variant_id,
                                               bialelicosimput.E28$FrecuenciasAlelicas.imputations,
                                               bialelicosimput.E28$FrecuenciasAlelicas.genomic_zone,
                                               bialelicosimput.E28$FrecuenciasAlelicas.gene.fragment,
                                               bialelicosimput.E28$FrecuenciasAlelicas.description,
                                               bialelicosimput.E28$FrecuenciasAlelicas.chromosome_id,
                                               bialelicosimput.E28$V27,
                                               bialelicosimput.E28$V28))

res.bialelicosimput.E31 <- distinct(data.frame(bialelicosimput.E31$FrecuenciasAlelicas.variant_id,
                                               bialelicosimput.E31$FrecuenciasAlelicas.imputations,
                                               bialelicosimput.E31$FrecuenciasAlelicas.genomic_zone,
                                               bialelicosimput.E31$FrecuenciasAlelicas.gene.fragment,
                                               bialelicosimput.E31$FrecuenciasAlelicas.description,
                                               bialelicosimput.E31$FrecuenciasAlelicas.chromosome_id,
                                               bialelicosimput.E31$V27,
                                               bialelicosimput.E31$V28))

res.bialelicosimput.E41 <- distinct(data.frame(bialelicosimput.E41$FrecuenciasAlelicas.variant_id,
                                               bialelicosimput.E41$FrecuenciasAlelicas.imputations,
                                               bialelicosimput.E41$FrecuenciasAlelicas.genomic_zone,
                                               bialelicosimput.E41$FrecuenciasAlelicas.gene.fragment,
                                               bialelicosimput.E41$FrecuenciasAlelicas.description,
                                               bialelicosimput.E41$FrecuenciasAlelicas.chromosome_id,
                                               bialelicosimput.E41$V27,
                                               bialelicosimput.E41$V28))

res.bialelicosimput.E96 <- distinct(data.frame(bialelicosimput.E96$FrecuenciasAlelicas.variant_id,
                                               bialelicosimput.E96$FrecuenciasAlelicas.imputations,
                                               bialelicosimput.E96$FrecuenciasAlelicas.genomic_zone,
                                               bialelicosimput.E96$FrecuenciasAlelicas.gene.fragment,
                                               bialelicosimput.E96$FrecuenciasAlelicas.description,
                                               bialelicosimput.E96$FrecuenciasAlelicas.chromosome_id,
                                               bialelicosimput.E96$V27,
                                               bialelicosimput.E96$V28))

res.bialelicosNOimput.FS16.S <- distinct(data.frame(bialelicosNOimput.FS16.S$FrecuenciasAlelicas.variant_id,
                                                    bialelicosNOimput.FS16.S$FrecuenciasAlelicas.imputations,
                                                    bialelicosNOimput.FS16.S$FrecuenciasAlelicas.genomic_zone,
                                                    bialelicosNOimput.FS16.S$FrecuenciasAlelicas.gene.fragment,
                                                    bialelicosNOimput.FS16.S$FrecuenciasAlelicas.description,
                                                    bialelicosNOimput.FS16.S$FrecuenciasAlelicas.chromosome_id,
                                                    bialelicosNOimput.FS16.S$V27,
                                                    bialelicosNOimput.FS16.S$V28))

res.bialelicosNOimput.FS19.S <- distinct(data.frame(bialelicosNOimput.FS19.S$FrecuenciasAlelicas.variant_id,
                                                    bialelicosNOimput.FS19.S$FrecuenciasAlelicas.imputations,
                                                    bialelicosNOimput.FS19.S$FrecuenciasAlelicas.genomic_zone,
                                                    bialelicosNOimput.FS19.S$FrecuenciasAlelicas.gene.fragment,
                                                    bialelicosNOimput.FS19.S$FrecuenciasAlelicas.description,
                                                    bialelicosNOimput.FS19.S$FrecuenciasAlelicas.chromosome_id,
                                                    bialelicosNOimput.FS19.S$V27,
                                                    bialelicosNOimput.FS19.S$V28))

res.bialelicosNOimput.FS20.S <- distinct(data.frame(bialelicosNOimput.FS20.S$FrecuenciasAlelicas.variant_id,
                                                    bialelicosNOimput.FS20.S$FrecuenciasAlelicas.imputations,
                                                    bialelicosNOimput.FS20.S$FrecuenciasAlelicas.genomic_zone,
                                                    bialelicosNOimput.FS20.S$FrecuenciasAlelicas.gene.fragment,
                                                    bialelicosNOimput.FS20.S$FrecuenciasAlelicas.description,
                                                    bialelicosNOimput.FS20.S$FrecuenciasAlelicas.chromosome_id,
                                                    bialelicosNOimput.FS20.S$V27,
                                                    bialelicosNOimput.FS20.S$V28))

res.bialelicosNOimput.FS22.S <- distinct(data.frame(bialelicosNOimput.FS22.S$FrecuenciasAlelicas.variant_id,
                                                    bialelicosNOimput.FS22.S$FrecuenciasAlelicas.imputations,
                                                    bialelicosNOimput.FS22.S$FrecuenciasAlelicas.genomic_zone,
                                                    bialelicosNOimput.FS22.S$FrecuenciasAlelicas.gene.fragment,
                                                    bialelicosNOimput.FS22.S$FrecuenciasAlelicas.description,
                                                    bialelicosNOimput.FS22.S$FrecuenciasAlelicas.chromosome_id,
                                                    bialelicosNOimput.FS22.S$V27,
                                                    bialelicosNOimput.FS22.S$V28))

res.bialelicosimput.FS16.S <- distinct(data.frame(bialelicosimput.FS16.S$FrecuenciasAlelicas.variant_id,
                                                  bialelicosimput.FS16.S$FrecuenciasAlelicas.imputations,
                                                  bialelicosimput.FS16.S$FrecuenciasAlelicas.genomic_zone,
                                                  bialelicosimput.FS16.S$FrecuenciasAlelicas.gene.fragment,
                                                  bialelicosimput.FS16.S$FrecuenciasAlelicas.description,
                                                  bialelicosimput.FS16.S$FrecuenciasAlelicas.chromosome_id,
                                                  bialelicosimput.FS16.S$V27,
                                                  bialelicosimput.FS16.S$V28))

res.bialelicosimput.FS19.S <- distinct(data.frame(bialelicosimput.FS19.S$FrecuenciasAlelicas.variant_id,
                                                  bialelicosimput.FS19.S$FrecuenciasAlelicas.imputations,
                                                  bialelicosimput.FS19.S$FrecuenciasAlelicas.genomic_zone,
                                                  bialelicosimput.FS19.S$FrecuenciasAlelicas.gene.fragment,
                                                  bialelicosimput.FS19.S$FrecuenciasAlelicas.description,
                                                  bialelicosimput.FS19.S$FrecuenciasAlelicas.chromosome_id,
                                                  bialelicosimput.FS19.S$V27,
                                                  bialelicosimput.FS19.S$V28))

res.bialelicosimput.FS20.S <- distinct(data.frame(bialelicosimput.FS20.S$FrecuenciasAlelicas.variant_id,
                                                  bialelicosimput.FS20.S$FrecuenciasAlelicas.imputations,
                                                  bialelicosimput.FS20.S$FrecuenciasAlelicas.genomic_zone,
                                                  bialelicosimput.FS20.S$FrecuenciasAlelicas.gene.fragment,
                                                  bialelicosimput.FS20.S$FrecuenciasAlelicas.description,
                                                  bialelicosimput.FS20.S$FrecuenciasAlelicas.chromosome_id,
                                                  bialelicosimput.FS20.S$V27,
                                                  bialelicosimput.FS20.S$V28))

res.bialelicosimput.FS22.S <- distinct(data.frame(bialelicosimput.FS22.S$FrecuenciasAlelicas.variant_id,
                                                  bialelicosimput.FS22.S$FrecuenciasAlelicas.imputations,
                                                  bialelicosimput.FS22.S$FrecuenciasAlelicas.genomic_zone,
                                                  bialelicosimput.FS22.S$FrecuenciasAlelicas.gene.fragment,
                                                  bialelicosimput.FS22.S$FrecuenciasAlelicas.description,
                                                  bialelicosimput.FS22.S$FrecuenciasAlelicas.chromosome_id,
                                                  bialelicosimput.FS22.S$V27,
                                                  bialelicosimput.FS22.S$V28))

colnames(res.bialelicosNOimput.EFS) <- colnames.df 
colnames(res.bialelicosNOimput.EN) <- colnames.df 
colnames(res.bialelicosNOimput.AL) <- colnames.df 
colnames(res.bialelicosimput.EFS) <- colnames.df 
colnames(res.bialelicosimput.EN) <- colnames.df 
colnames(res.bialelicosimput.AL) <- colnames.df 
colnames(res.bialelicosNOimput.A05) <- colnames.df 
colnames(res.bialelicosNOimput.A07) <- colnames.df 
colnames(res.bialelicosNOimput.A09) <- colnames.df 
colnames(res.bialelicosNOimput.A10) <- colnames.df 
colnames(res.bialelicosimput.A05) <- colnames.df 
colnames(res.bialelicosimput.A07) <- colnames.df 
colnames(res.bialelicosimput.A09) <- colnames.df 
colnames(res.bialelicosimput.A10) <- colnames.df 
colnames(res.bialelicosNOimput.E28) <- colnames.df 
colnames(res.bialelicosNOimput.E31) <- colnames.df 
colnames(res.bialelicosNOimput.E41) <- colnames.df 
colnames(res.bialelicosNOimput.E96) <- colnames.df
colnames(res.bialelicosimput.E28) <- colnames.df 
colnames(res.bialelicosimput.E31) <- colnames.df 
colnames(res.bialelicosimput.E41) <- colnames.df 
colnames(res.bialelicosimput.E96) <- colnames.df 
colnames(res.bialelicosNOimput.FS16.S) <- colnames.df 
colnames(res.bialelicosNOimput.FS19.S) <- colnames.df 
colnames(res.bialelicosNOimput.FS20.S) <- colnames.df 
colnames(res.bialelicosNOimput.FS22.S) <- colnames.df 
colnames(res.bialelicosimput.FS16.S) <- colnames.df 
colnames(res.bialelicosimput.FS19.S) <- colnames.df 
colnames(res.bialelicosimput.FS20.S) <- colnames.df 
colnames(res.bialelicosimput.FS22.S) <- colnames.df 

#Now, we subest the above data.frames to extract
#only those SNP variants with a p-value < 0.05
sig.res.bialelicosNOimput.EFS <- subset(res.bialelicosNOimput.EFS,
                                        res.bialelicosNOimput.EFS$pvalue <0.05)
sig.res.bialelicosNOimput.EN <- subset(res.bialelicosNOimput.EN,
                                        res.bialelicosNOimput.EN$pvalue <0.05)
sig.res.bialelicosNOimput.AL <- subset(res.bialelicosNOimput.AL,
                                       res.bialelicosNOimput.AL$pvalue <0.05)
sig.res.bialelicosNOimput.A05 <- subset(res.bialelicosNOimput.A05,
                                        res.bialelicosNOimput.A05$pvalue <0.05)
sig.res.bialelicosNOimput.A07 <- subset(res.bialelicosNOimput.A07,
                                        res.bialelicosNOimput.A07$pvalue <0.05)
sig.res.bialelicosNOimput.A09 <- subset(res.bialelicosNOimput.A09,
                                        res.bialelicosNOimput.A09$pvalue <0.05)
sig.res.bialelicosNOimput.A10 <- subset(res.bialelicosNOimput.A10,
                                        res.bialelicosNOimput.A10$pvalue <0.05)
sig.res.bialelicosNOimput.E28 <- subset(res.bialelicosNOimput.E28,
                                        res.bialelicosNOimput.E28$pvalue <0.05)
sig.res.bialelicosNOimput.E31 <- subset(res.bialelicosNOimput.E31,
                                        res.bialelicosNOimput.E31$pvalue <0.05)
sig.res.bialelicosNOimput.E41 <- subset(res.bialelicosNOimput.E41,
                                        res.bialelicosNOimput.E41$pvalue <0.05)
sig.res.bialelicosNOimput.E96 <- subset(res.bialelicosNOimput.E96,
                                        res.bialelicosNOimput.E96$pvalue <0.05)
sig.res.bialelicosNOimput.FS16.S <- subset(res.bialelicosNOimput.FS16.S,
                                           res.bialelicosNOimput.FS16.S$pvalue <0.05)
sig.res.bialelicosNOimput.FS19.S <- subset(res.bialelicosNOimput.FS19.S,
                                           res.bialelicosNOimput.FS19.S$pvalue <0.05)
sig.res.bialelicosNOimput.FS20.S <- subset(res.bialelicosNOimput.FS20.S,
                                           res.bialelicosNOimput.FS20.S$pvalue <0.05)
sig.res.bialelicosNOimput.FS22.S <- subset(res.bialelicosNOimput.FS22.S,
                                           res.bialelicosNOimput.FS22.S$pvalue <0.05)
sig.res.bialelicosimput.EFS <- subset(res.bialelicosimput.EFS,
                                      res.bialelicosimput.EFS$pvalue <0.05)
sig.res.bialelicosimput.EN <- subset(res.bialelicosimput.EN,
                                     res.bialelicosimput.EN$pvalue <0.05)
sig.res.bialelicosimput.AL <- subset(res.bialelicosimput.AL,
                                     res.bialelicosimput.AL$pvalue <0.05)
sig.res.bialelicosimput.A05 <- subset(res.bialelicosimput.A05,
                                      res.bialelicosimput.A05$pvalue <0.05)
sig.res.bialelicosimput.A07 <- subset(res.bialelicosimput.A07,
                                      res.bialelicosimput.A07$pvalue <0.05)
sig.res.bialelicosimput.A09 <- subset(res.bialelicosimput.A09,
                                      res.bialelicosimput.A09$pvalue <0.05)
sig.res.bialelicosimput.A10 <- subset(res.bialelicosimput.A10,
                                      res.bialelicosimput.A10$pvalue <0.05)
sig.res.bialelicosimput.E28 <- subset(res.bialelicosimput.E28,
                                      res.bialelicosimput.E28$pvalue <0.05)
sig.res.bialelicosimput.E31 <- subset(res.bialelicosimput.E31,
                                      res.bialelicosimput.E31$pvalue <0.05)
sig.res.bialelicosimput.E41 <- subset(res.bialelicosimput.E41,
                                      res.bialelicosimput.E41$pvalue <0.05)
sig.res.bialelicosimput.E96 <- subset(res.bialelicosimput.E96,
                                      res.bialelicosimput.E96$pvalue <0.05)
sig.res.bialelicosimput.FS16.S <- subset(res.bialelicosimput.FS16.S,
                                         res.bialelicosimput.FS16.S$pvalue <0.05)
sig.res.bialelicosimput.FS19.S <- subset(res.bialelicosimput.FS19.S,
                                         res.bialelicosimput.FS19.S$pvalue <0.05)
sig.res.bialelicosimput.FS20.S <- subset(res.bialelicosimput.FS20.S,
                                         res.bialelicosimput.FS20.S$pvalue <0.05)
sig.res.bialelicosimput.FS22.S <- subset(res.bialelicosimput.FS22.S,
                                         res.bialelicosimput.FS22.S$pvalue <0.05)

##We save the results in tab separated files
write.table(sig.res.bialelicosNOimput.EFS ,
            "sig.res.bialelicosNOimput-EFS.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.EN ,
            "sig.res.bialelicosNOimput-EN.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.AL ,
            "sig.res.bialelicosNOimput-AL.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.A05 ,
            "sig.res.bialelicosNOimput-A05.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.A07 ,
            "sig.res.bialelicosNOimput-A07.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.A09 ,
            "sig.res.bialelicosNOimput-A09.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.A10 ,
            "sig.res.bialelicosNOimput-A10.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.E28 ,
            "sig.res.bialelicosNOimput-E28.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.E31 ,
            "sig.res.bialelicosNOimput-E31.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.E41 ,
            "sig.res.bialelicosNOimput-E41.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.E96 ,
            "sig.res.bialelicosNOimput-E96.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.FS16.S ,
            "sig.res.bialelicosNOimput-FS16-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.FS19.S ,
            "sig.res.bialelicosNOimput-FS19-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.FS20.S ,
            "sig.res.bialelicosNOimput-FS20-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosNOimput.FS22.S ,
            "sig.res.bialelicosNOimput-FS22-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.EFS ,
            "sig.res.bialelicosimput-EFS.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.EN ,
            "sig.res.bialelicosimput-EN.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.AL ,
            "sig.res.bialelicosimput-AL.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.A05 ,
            "sig.res.bialelicosimput-A05.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.A07 ,
            "sig.res.bialelicosimput-A07.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.A09 ,
            "sig.res.bialelicosimput-A09.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.A10 ,
            "sig.res.bialelicosimput-A10.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.E28 ,
            "sig.res.bialelicosimput-E28.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.E31 ,
            "sig.res.bialelicosimput-E31.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.E41 ,
            "sig.res.bialelicosimput-E41.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.E96 ,
            "sig.res.bialelicosimput-E96.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.FS16.S ,
            "sig.res.bialelicosimput-FS16-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.FS19.S ,
            "sig.res.bialelicosimput-FS19-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.FS20.S ,
            "sig.res.bialelicosimput-FS20-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

write.table(sig.res.bialelicosimput.FS22.S ,
            "sig.res.bialelicosimput-FS22-S.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)


###########################################################################################################


###############################
#                             #
#        Venn diagrams        #
#                             #
###############################

#We will use the ggvenn library to produce Venn Diagrams in pdf
#for SNP and for genes. We will intersect the SNP for the four
#hybrids, the four suber and the four ilex (imputed and Non-imputed separately)
library(ggvenn)
#SNP-hybrid families
x <- list(FS16 = sig.res.bialelicosNOimput.FS16.S$variant_id, 
          FS19 = sig.res.bialelicosNOimput.FS19.S$variant_id, 
          FS20 = sig.res.bialelicosNOimput.FS20.S$variant_id,
          FS22 = sig.res.bialelicosNOimput.FS22.S$variant_id)

pdf("NOimput-FS16-19-20-22-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()

x <- list(FS16 = sig.res.bialelicosimput.FS16.S$variant_id, 
          FS19 = sig.res.bialelicosimput.FS19.S$variant_id, 
          FS20 = sig.res.bialelicosimput.FS20.S$variant_id,
          FS22 = sig.res.bialelicosimput.FS22.S$variant_id)

pdf("Imput-FS16-19-20-22-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()

#Gene-hybrid families

x <- list(FS16 = sig.res.bialelicosNOimput.FS16.S$gene, 
          FS19 = sig.res.bialelicosNOimput.FS19.S$gene, 
          FS20 = sig.res.bialelicosNOimput.FS20.S$gene,
          FS22 = sig.res.bialelicosNOimput.FS22.S$gene)

pdf("NOimput-fragments-FS16-19-20-22-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()

x <- list(FS16 = sig.res.bialelicosimput.FS16.S$gene, 
          FS19 = sig.res.bialelicosimput.FS19.S$gene, 
          FS20 = sig.res.bialelicosimput.FS20.S$gene,
          FS22 = sig.res.bialelicosimput.FS22.S$gene)

pdf("Imput-fragments-FS16-19-20-22-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()

#SNP-suber families
x <- list(A05 = sig.res.bialelicosNOimput.A05$variant_id, 
          A07 = sig.res.bialelicosNOimput.A07$variant_id, 
          A09 = sig.res.bialelicosNOimput.A09$variant_id,
          A10 = sig.res.bialelicosNOimput.A10$variant_id)

pdf("NOimput-A05-A07-A09-A10-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()


x <- list(A05 = sig.res.bialelicosimput.A05$variant_id, 
          A07 = sig.res.bialelicosimput.A07$variant_id, 
          A09 = sig.res.bialelicosimput.A09$variant_id,
          A10 = sig.res.bialelicosimput.A10$variant_id)

pdf("Imput-A05-A07-A09-A10-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()

#Gene-Suber families
x <- list(A05 = sig.res.bialelicosNOimput.A05$gene, 
          A07 = sig.res.bialelicosNOimput.A07$gene, 
          A09 = sig.res.bialelicosNOimput.A09$gene,
          A10 = sig.res.bialelicosNOimput.A10$gene)

pdf("NOimput-fragments-A05-A07-A09-A10-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()

x <- list(A05 = sig.res.bialelicosimput.A05$gene, 
          A07 = sig.res.bialelicosimput.A07$gene, 
          A09 = sig.res.bialelicosimput.A09$gene,
          A10 = sig.res.bialelicosimput.A10$gene)

pdf("Imput-fragments-A05-A07-A09-A10-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()

#SNP-Ilex families
x <- list(E28 = sig.res.bialelicosNOimput.E28$variant_id, 
          E31 = sig.res.bialelicosNOimput.E31$variant_id, 
          E41 = sig.res.bialelicosNOimput.E41$variant_id,
          E96 = sig.res.bialelicosNOimput.E96$variant_id)

pdf("NOimput-E28-E31-E41-E96-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()


x <- list(E28 = sig.res.bialelicosimput.E28$variant_id, 
          E31 = sig.res.bialelicosimput.E31$variant_id, 
          E41 = sig.res.bialelicosimput.E41$variant_id,
          E96 = sig.res.bialelicosimput.E96$variant_id)

pdf("Imput-E28-E31-E41-E96-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()

#Gene-Ilex families
x <- list(E28 = sig.res.bialelicosNOimput.E28$gene, 
          E31 = sig.res.bialelicosNOimput.E31$gene, 
          E41 = sig.res.bialelicosNOimput.E41$gene,
          E96 = sig.res.bialelicosNOimput.E96$gene)

pdf("NOimput-fragments-E28-E31-E41-E96-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()


x <- list(E28 = sig.res.bialelicosimput.E28$gene, 
          E31 = sig.res.bialelicosimput.E31$gene, 
          E41 = sig.res.bialelicosimput.E41$gene,
          E96 = sig.res.bialelicosimput.E96$gene)

pdf("Imput-fragments-E28-E31-E41-E96-comunes-significativos.pdf",width=28, height = 28)
print(
  ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 18,text_size = 20
)
)
dev.off()

######################################
#                                    #
#     Stats of significant data      #  
#                                    #
######################################

#We will produce a numerical tab separated file with the number of significant SNP
#and fragment in each case

rows<- c("EFS",
          "AL",
          "EN",
          "A05",
          "A07",
          "A09",
          "A10",
          "E28",
          "E31",
          "E41",
          "E96",
          "FS16",
          "FS19",
          "FS20",
          "FS22",
          "Hybrid progenies",
          "Suber progenies",
          "Ilex progenies", 
          "All progenies")

cols <- c("Sig SNP non-imp",
          "Sig SNP imp")

######################################################################################
#We already have the individual families data, we need the intersected SNP
#Hybrid families
SNP.NOimput.FS16.FS19.FS20.FS22 <- Reduce(intersect, list(sig.res.bialelicosNOimput.FS16.S$variant_id,
                                                          sig.res.bialelicosNOimput.FS19.S$variant_id,
                                                          sig.res.bialelicosNOimput.FS20.S$variant_id,
                                                          sig.res.bialelicosNOimput.FS22.S$variant_id))

sig.res.bialelicosNOimput.FS16.FS19.FS20.FS22.SNP <-sig.res.bialelicosNOimput.FS22.S[sig.res.bialelicosNOimput.FS22.S$variant_id %in% SNP.NOimput.FS16.FS19.FS20.FS22,]


write.table(sig.res.bialelicosNOimput.FS16.FS19.FS20.FS22.SNP,
            "sig-res-bialelicosNOimput-FS16-FS19-FS20-FS22-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
SNP.imput.FS16.FS19.FS20.FS22 <- Reduce(intersect, list(sig.res.bialelicosimput.FS16.S$variant_id,
                                                        sig.res.bialelicosimput.FS19.S$variant_id,
                                                        sig.res.bialelicosimput.FS20.S$variant_id,
                                                        sig.res.bialelicosimput.FS22.S$variant_id))

sig.res.bialelicosimput.FS16.FS19.FS20.FS22.SNP <-sig.res.bialelicosimput.FS22.S[sig.res.bialelicosimput.FS22.S$variant_id %in% SNP.imput.FS16.FS19.FS20.FS22,]


write.table(sig.res.bialelicosimput.FS16.FS19.FS20.FS22.SNP,
            "sig-res-bialelicosimput-FS16-FS19-FS20-FS22-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

#suber families
SNP.NOimput.A05.A07.A09.A10 <- Reduce(intersect, list(sig.res.bialelicosNOimput.A05$variant_id,
                                                      sig.res.bialelicosNOimput.A07$variant_id,
                                                      sig.res.bialelicosNOimput.A09$variant_id,
                                                      sig.res.bialelicosNOimput.A10$variant_id))

sig.res.bialelicosNOimput.A05.A07.A09.A10.SNP <-sig.res.bialelicosNOimput.A10[sig.res.bialelicosNOimput.A10$variant_id %in% SNP.NOimput.A05.A07.A09.A10,]
write.table(sig.res.bialelicosNOimput.A05.A07.A09.A10.SNP,
            "sig-res-bialelicosNOimput-A05-A07-A09-A10-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

SNP.imput.A05.A07.A09.A10 <- Reduce(intersect, list(sig.res.bialelicosimput.A05$variant_id,
                                                    sig.res.bialelicosimput.A07$variant_id,
                                                    sig.res.bialelicosimput.A09$variant_id,
                                                    sig.res.bialelicosimput.A10$variant_id))


sig.res.bialelicosimput.A05.A07.A09.A10.SNP <-sig.res.bialelicosimput.A10[sig.res.bialelicosimput.A10$variant_id %in% SNP.imput.A05.A07.A09.A10,]
write.table(sig.res.bialelicosimput.A05.A07.A09.A10.SNP,
            "sig-res-bialelicosimput-A05-A07-A09-A10-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

#ilex families
SNP.NOimput.E28.E31.E41.E96 <- Reduce(intersect, list(sig.res.bialelicosNOimput.E28$variant_id,
                                                      sig.res.bialelicosNOimput.E31$variant_id,
                                                      sig.res.bialelicosNOimput.E41$variant_id,
                                                      sig.res.bialelicosNOimput.E96$variant_id))

sig.res.bialelicosNOimput.E28.E31.E41.E96.SNP <-sig.res.bialelicosNOimput.E96[sig.res.bialelicosNOimput.E96$variant_id %in% SNP.NOimput.E28.E31.E41.E96,]

write.table(sig.res.bialelicosNOimput.E28.E31.E41.E96.SNP,
            "sig-res-bialelicosNOimput-E28-E31-E41-E96-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)
SNP.imput.E28.E31.E41.E96 <- Reduce(intersect, list(sig.res.bialelicosimput.E28$variant_id,
                                                    sig.res.bialelicosimput.E31$variant_id,
                                                    sig.res.bialelicosimput.E41$variant_id,
                                                    sig.res.bialelicosimput.E96$variant_id))

sig.res.bialelicosimput.E28.E31.E41.E96.SNP <-sig.res.bialelicosimput.E96[sig.res.bialelicosimput.E96$variant_id %in% SNP.imput.E28.E31.E41.E96,]

write.table(sig.res.bialelicosimput.E28.E31.E41.E96.SNP,
            "sig-res-bialelicosimput-E28-E31-E41-E96-SNP.txt",
            sep="\t",
            dec=",",
            row.names = FALSE)

#All families
##Progenies alcornoques con mayor número de muestras (NO IMPUTADOS):
SNP.NOimput.A05.A07.A09.A10.E28.E31.E41.E96.FS16.FS19.FS20.FS22 <- Reduce(intersect, list(SNP.NOimput.FS16.FS19.FS20.FS22,
                                                                                          SNP.NOimput.A05.A07.A09.A10,
                                                                                          SNP.NOimput.E28.E31.E41.E96))
length(SNP.NOimput.A05.A07.A09.A10.E28.E31.E41.E96.FS16.FS19.FS20.FS22)
SNP.imput.A05.A07.A09.A10.E28.E31.E41.E96.FS16.FS19.FS20.FS22 <- Reduce(intersect, list(SNP.imput.FS16.FS19.FS20.FS22,
                                                                                        SNP.imput.A05.A07.A09.A10,
                                                                                        SNP.imput.E28.E31.E41.E96))




################################################################################################

#We fill in the data.frame creating vector with the number of SNP
sig.SNP.NOimput <- c(length(sig.res.bialelicosNOimput.EFS$variant_id),
                     length(sig.res.bialelicosNOimput.AL$variant_id),
                     length(sig.res.bialelicosNOimput.EN$variant_id),
                     length(sig.res.bialelicosNOimput.A05$variant_id),
                     length(sig.res.bialelicosNOimput.A07$variant_id),
                     length(sig.res.bialelicosNOimput.A09$variant_id),
                     length(sig.res.bialelicosNOimput.A10$variant_id),
                     length(sig.res.bialelicosNOimput.E28$variant_id),
                     length(sig.res.bialelicosNOimput.E31$variant_id),
                     length(sig.res.bialelicosNOimput.E41$variant_id),
                     length(sig.res.bialelicosNOimput.E96$variant_id),
                     length(sig.res.bialelicosNOimput.FS16.S$variant_id),
                     length(sig.res.bialelicosNOimput.FS19.S$variant_id),
                     length(sig.res.bialelicosNOimput.FS20.S$variant_id),
                     length(sig.res.bialelicosNOimput.FS22.S$variant_id),
                     length(SNP.NOimput.FS16.FS19.FS20.FS22),
                     length(SNP.NOimput.A05.A07.A09.A10),
                     length(SNP.NOimput.E28.E31.E41.E96),
                     length(SNP.NOimput.A05.A07.A09.A10.E28.E31.E41.E96.FS16.FS19.FS20.FS22))
sig.SNP.imput <- c(length(sig.res.bialelicosimput.EFS$variant_id),
                   length(sig.res.bialelicosimput.AL$variant_id),
                   length(sig.res.bialelicosimput.EN$variant_id),
                   length(sig.res.bialelicosimput.A05$variant_id),
                   length(sig.res.bialelicosimput.A07$variant_id),
                   length(sig.res.bialelicosimput.A09$variant_id),
                   length(sig.res.bialelicosimput.A10$variant_id),
                   length(sig.res.bialelicosimput.E28$variant_id),
                   length(sig.res.bialelicosimput.E31$variant_id),
                   length(sig.res.bialelicosimput.E41$variant_id),
                   length(sig.res.bialelicosimput.E96$variant_id),
                   length(sig.res.bialelicosimput.FS16.S$variant_id),
                   length(sig.res.bialelicosimput.FS19.S$variant_id),
                   length(sig.res.bialelicosimput.FS20.S$variant_id),
                   length(sig.res.bialelicosimput.FS22.S$variant_id),
                   length(SNP.imput.FS16.FS19.FS20.FS22),
                   length(SNP.imput.A05.A07.A09.A10),
                   length(SNP.imput.E28.E31.E41.E96),
                   length(SNP.imput.A05.A07.A09.A10.E28.E31.E41.E96.FS16.FS19.FS20.FS22))

significant.stats <- data.frame(sig.SNP.NOimput,
                                sig.SNP.imput)

colnames(significant.stats) <- cols

rownames(significant.stats) <- rows

write.table(significant.stats,
            "significant.stats.txt",
            sep="\t",
            dec=",",
            row.names = TRUE,
            col.names = TRUE)

##################
#                #
#       END      #
#                #
##################




