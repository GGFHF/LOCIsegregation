###############################################
#                                             #
#    R script to find correspondence          #
#   between Linkage and physical maps         #
#                                             #
###############################################
#This software has been developed by:

#  Dpto. Sistemas y Recursos Naturales
#  ETSI Montes, Forestal y del Medio Natural
#  Universidad Politecnica de Madrid
# https://github.com/ggfhf/

#We set the working directory to run the scripts. Be sure all the scripts are within this directory
setwd(WD.input.0)

#We upload the necessary libraries
#source('0-Config.R')
#We execute the necessary scripts (hosted in the same working directory)
#The WD in the scripts must be the same that in this script

#source('1-Chi2-test-adults.R')
#source('2-Chi2-test-suber-ilex-progenies.R')
#source('3-Chi2-test-hybrid-progenies.R')
#source('4-Significant-SNP-Venn-Diagrams.R')
source('5-Crossing-RNAseq-data.R')

#We set the working directory. Be sure all the input files are within this directory
setwd(WD.input)


###We scan the correspondence between linkage groups and physical positions
#in the Q. robur genome assembly for each gene found in ScnI using input_LMAP_viewer.txt"

input.data<-read.table("input_LMAP_viewer4.txt",
                       header=TRUE,
                       sep="\t",
                       dec=",")

#We set the working directory to produce the output.
setwd(WD.output.3)

#Significant genes common to FS16, FS19, FS20 and FS22
#We create the list of the markers that will appear 
flist <- list()
locus <- as.vector(full.inner.sig.res.bialelicos.FS16.FS19.FS20.FS22.SNP$gene)
font  <- c(4)   #italic
col <- c("blue")
flist[[1]] <- list(locus = locus, font = font, col = col)

#For Chr1-Chr3
lmv.linkage.plot(input.data,
                 "ScnI-HYBRID-FS16-FS19-FS20-FS22-Chr01-03.pdf",
                 mapthese = c("LG1","LG2","LG3"),
                 at.axis = NULL,
                 autoconnadj = TRUE,
                 cex.axis = 1,
                 cex.lgtitle = par("cex.main"),
                 cex.main = par("cex.main"),
                 col.axis = par("col.axis"),
                 col.lgtitle = par("col.main"),
                 col.main = par("col.main"),
                 conndf = NULL,
                 denmap = FALSE,
                 dupnbr = FALSE,
                 font.axis = par("font.axis"),
                 font.lgtitle = par("font.main"),
                 font.main = par("font.main"),
                 header = TRUE,
                 labdist = 0.3,
                 labels.axis = TRUE,
                 lcex = par("cex"),
                 lcol = par("col"),
                 lfont = par("font"),
                 lgperrow = NULL,
                 lgtitles = NULL,
                 lgw = 0.25,
                 lg.col =  "lightblue1",
                 lg.lwd = par("lwd"),
                 lty.axis = "solid",
                 lwd.axis = 1,
                 lwd.ticks.axis = lwd.axis,
                 main = NULL,
                 markerformatlist = flist,
                 maxnbrcolsfordups = 2,
                 pdf.bg = "transparent",
                 pdf.family = "Helvetica",
                 pdf.fg = "black",
                 pdf.width = NULL,
                 pdf.height = NULL,
                 pdf.pointsize = 12,
                 pdf.title = "LinkageMapView R output",
                 posonleft = NULL,
                 prtlgtitles = TRUE,
                 qtldf = NULL,
                 revthese = NULL,
                 rcex = par("cex"),
                 rcol = par("col"),
                 rfont = par("font"),
                 roundpos = 1,
                 rsegcol = TRUE,
                 ruler = FALSE,
                 sectcoldf = NULL,
                 segcol = NULL,
                 qtlscanone = NULL,
                 showonly = NULL,
                 units = "cM",
                 ylab = units
)


#For Chr4-Chr6
lmv.linkage.plot(input.data,
                 "ScnI-HYBRID-FS16-FS19-FS20-FS22-Chr04-06.pdf",
                 mapthese = c("LG4","LG5","LG6"),
                 at.axis = NULL,
                 autoconnadj = TRUE,
                 cex.axis = 1,
                 cex.lgtitle = par("cex.main"),
                 cex.main = par("cex.main"),
                 col.axis = par("col.axis"),
                 col.lgtitle = par("col.main"),
                 col.main = par("col.main"),
                 conndf = NULL,
                 denmap = FALSE,
                 dupnbr = FALSE,
                 font.axis = par("font.axis"),
                 font.lgtitle = par("font.main"),
                 font.main = par("font.main"),
                 header = TRUE,
                 labdist = 0.3,
                 labels.axis = TRUE,
                 lcex = par("cex"),
                 lcol = par("col"),
                 lfont = par("font"),
                 lgperrow = NULL,
                 lgtitles = NULL,
                 lgw = 0.25,
                 lg.col =  "lightblue1",
                 lg.lwd = par("lwd"),
                 lty.axis = "solid",
                 lwd.axis = 1,
                 lwd.ticks.axis = lwd.axis,
                 main = NULL,
                 markerformatlist = flist,
                 maxnbrcolsfordups = 2,
                 pdf.bg = "transparent",
                 pdf.family = "Helvetica",
                 pdf.fg = "black",
                 pdf.width = NULL,
                 pdf.height = NULL,
                 pdf.pointsize = 12,
                 pdf.title = "LinkageMapView R output",
                 posonleft = NULL,
                 prtlgtitles = TRUE,
                 qtldf = NULL,
                 revthese = NULL,
                 rcex = par("cex"),
                 rcol = par("col"),
                 rfont = par("font"),
                 roundpos = 1,
                 rsegcol = TRUE,
                 ruler = FALSE,
                 sectcoldf = NULL,
                 segcol = NULL,
                 qtlscanone = NULL,
                 showonly = NULL,
                 units = "cM",
                 ylab = units
)

#For Chr7-Chr9
lmv.linkage.plot(input.data,
                 "ScnI-HYBRID-FS16-FS19-FS20-FS22-Chr07-09.pdf",
                 mapthese = c("LG7","LG8","LG9"),
                 at.axis = NULL,
                 autoconnadj = TRUE,
                 cex.axis = 1,
                 cex.lgtitle = par("cex.main"),
                 cex.main = par("cex.main"),
                 col.axis = par("col.axis"),
                 col.lgtitle = par("col.main"),
                 col.main = par("col.main"),
                 conndf = NULL,
                 denmap = FALSE,
                 dupnbr = FALSE,
                 font.axis = par("font.axis"),
                 font.lgtitle = par("font.main"),
                 font.main = par("font.main"),
                 header = TRUE,
                 labdist = 0.3,
                 labels.axis = TRUE,
                 lcex = par("cex"),
                 lcol = par("col"),
                 lfont = par("font"),
                 lgperrow = NULL,
                 lgtitles = NULL,
                 lgw = 0.25,
                 lg.col =  "lightblue1",
                 lg.lwd = par("lwd"),
                 lty.axis = "solid",
                 lwd.axis = 1,
                 lwd.ticks.axis = lwd.axis,
                 main = NULL,
                 markerformatlist = flist,
                 maxnbrcolsfordups = 2,
                 pdf.bg = "transparent",
                 pdf.family = "Helvetica",
                 pdf.fg = "black",
                 pdf.width = NULL,
                 pdf.height = NULL,
                 pdf.pointsize = 12,
                 pdf.title = "LinkageMapView R output",
                 posonleft = NULL,
                 prtlgtitles = TRUE,
                 qtldf = NULL,
                 revthese = NULL,
                 rcex = par("cex"),
                 rcol = par("col"),
                 rfont = par("font"),
                 roundpos = 1,
                 rsegcol = TRUE,
                 ruler = FALSE,
                 sectcoldf = NULL,
                 segcol = NULL,
                 qtlscanone = NULL,
                 showonly = NULL,
                 units = "cM",
                 ylab = units
)

#For Chr10-Chr12
lmv.linkage.plot(input.data,
                 "ScnI-HYBRID-FS16-FS19-FS20-FS22-Chr10-12.pdf",
                 mapthese = c("LG10","LG11","LG12"),
                 at.axis = NULL,
                 autoconnadj = TRUE,
                 cex.axis = 1,
                 cex.lgtitle = par("cex.main"),
                 cex.main = par("cex.main"),
                 col.axis = par("col.axis"),
                 col.lgtitle = par("col.main"),
                 col.main = par("col.main"),
                 conndf = NULL,
                 denmap = FALSE,
                 dupnbr = FALSE,
                 font.axis = par("font.axis"),
                 font.lgtitle = par("font.main"),
                 font.main = par("font.main"),
                 header = TRUE,
                 labdist = 0.3,
                 labels.axis = TRUE,
                 lcex = par("cex"),
                 lcol = par("col"),
                 lfont = par("font"),
                 lgperrow = NULL,
                 lgtitles = NULL,
                 lgw = 0.25,
                 lg.col =  "lightblue1",
                 lg.lwd = par("lwd"),
                 lty.axis = "solid",
                 lwd.axis = 1,
                 lwd.ticks.axis = lwd.axis,
                 main = NULL,
                 markerformatlist = flist,
                 maxnbrcolsfordups = 2,
                 pdf.bg = "transparent",
                 pdf.family = "Helvetica",
                 pdf.fg = "black",
                 pdf.width = NULL,
                 pdf.height = NULL,
                 pdf.pointsize = 12,
                 pdf.title = "LinkageMapView R output",
                 posonleft = NULL,
                 prtlgtitles = TRUE,
                 qtldf = NULL,
                 revthese = NULL,
                 rcex = par("cex"),
                 rcol = par("col"),
                 rfont = par("font"),
                 roundpos = 1,
                 rsegcol = TRUE,
                 ruler = FALSE,
                 sectcoldf = NULL,
                 segcol = NULL,
                 qtlscanone = NULL,
                 showonly = NULL,
                 units = "cM",
                 ylab = units
)




#Significant genes common to EFS
#We create the list of the markers that will appear 
flist <- list()
locus <- as.vector(full.inner.sig.res.bialelicos.EFS.SNP$gene)
font  <- c(4)   #italic
col <- c("blue")
flist[[1]] <- list(locus = locus, font = font, col = col)

#For Chr1-Chr3
lmv.linkage.plot(input.data,
                 "ScnI-HYBRID-EFS-FS22-Chr01-03.pdf",
                 mapthese = c("LG1","LG2","LG3"),
                 at.axis = NULL,
                 autoconnadj = TRUE,
                 cex.axis = 1,
                 cex.lgtitle = par("cex.main"),
                 cex.main = par("cex.main"),
                 col.axis = par("col.axis"),
                 col.lgtitle = par("col.main"),
                 col.main = par("col.main"),
                 conndf = NULL,
                 denmap = FALSE,
                 dupnbr = FALSE,
                 font.axis = par("font.axis"),
                 font.lgtitle = par("font.main"),
                 font.main = par("font.main"),
                 header = TRUE,
                 labdist = 0.3,
                 labels.axis = TRUE,
                 lcex = par("cex"),
                 lcol = par("col"),
                 lfont = par("font"),
                 lgperrow = NULL,
                 lgtitles = NULL,
                 lgw = 0.25,
                 lg.col =  "lightblue1",
                 lg.lwd = par("lwd"),
                 lty.axis = "solid",
                 lwd.axis = 1,
                 lwd.ticks.axis = lwd.axis,
                 main = NULL,
                 markerformatlist = flist,
                 maxnbrcolsfordups = 2,
                 pdf.bg = "transparent",
                 pdf.family = "Helvetica",
                 pdf.fg = "black",
                 pdf.width = NULL,
                 pdf.height = NULL,
                 pdf.pointsize = 12,
                 pdf.title = "LinkageMapView R output",
                 posonleft = NULL,
                 prtlgtitles = TRUE,
                 qtldf = NULL,
                 revthese = NULL,
                 rcex = par("cex"),
                 rcol = par("col"),
                 rfont = par("font"),
                 roundpos = 1,
                 rsegcol = TRUE,
                 ruler = FALSE,
                 sectcoldf = NULL,
                 segcol = NULL,
                 qtlscanone = NULL,
                 showonly = NULL,
                 units = "cM",
                 ylab = units
)


#For Chr4-Chr6
lmv.linkage.plot(input.data,
                 "ScnI-HYBRID-EFS-FS22-Chr04-06.pdf",
                 mapthese = c("LG4","LG5","LG6"),
                 at.axis = NULL,
                 autoconnadj = TRUE,
                 cex.axis = 1,
                 cex.lgtitle = par("cex.main"),
                 cex.main = par("cex.main"),
                 col.axis = par("col.axis"),
                 col.lgtitle = par("col.main"),
                 col.main = par("col.main"),
                 conndf = NULL,
                 denmap = FALSE,
                 dupnbr = FALSE,
                 font.axis = par("font.axis"),
                 font.lgtitle = par("font.main"),
                 font.main = par("font.main"),
                 header = TRUE,
                 labdist = 0.3,
                 labels.axis = TRUE,
                 lcex = par("cex"),
                 lcol = par("col"),
                 lfont = par("font"),
                 lgperrow = NULL,
                 lgtitles = NULL,
                 lgw = 0.25,
                 lg.col =  "lightblue1",
                 lg.lwd = par("lwd"),
                 lty.axis = "solid",
                 lwd.axis = 1,
                 lwd.ticks.axis = lwd.axis,
                 main = NULL,
                 markerformatlist = flist,
                 maxnbrcolsfordups = 2,
                 pdf.bg = "transparent",
                 pdf.family = "Helvetica",
                 pdf.fg = "black",
                 pdf.width = NULL,
                 pdf.height = NULL,
                 pdf.pointsize = 12,
                 pdf.title = "LinkageMapView R output",
                 posonleft = NULL,
                 prtlgtitles = TRUE,
                 qtldf = NULL,
                 revthese = NULL,
                 rcex = par("cex"),
                 rcol = par("col"),
                 rfont = par("font"),
                 roundpos = 1,
                 rsegcol = TRUE,
                 ruler = FALSE,
                 sectcoldf = NULL,
                 segcol = NULL,
                 qtlscanone = NULL,
                 showonly = NULL,
                 units = "cM",
                 ylab = units
)

#For Chr7-Chr9
lmv.linkage.plot(input.data,
                 "ScnI-HYBRID-EFS-FS22-Chr07-09.pdf",
                 mapthese = c("LG7","LG8","LG9"),
                 at.axis = NULL,
                 autoconnadj = TRUE,
                 cex.axis = 1,
                 cex.lgtitle = par("cex.main"),
                 cex.main = par("cex.main"),
                 col.axis = par("col.axis"),
                 col.lgtitle = par("col.main"),
                 col.main = par("col.main"),
                 conndf = NULL,
                 denmap = FALSE,
                 dupnbr = FALSE,
                 font.axis = par("font.axis"),
                 font.lgtitle = par("font.main"),
                 font.main = par("font.main"),
                 header = TRUE,
                 labdist = 0.3,
                 labels.axis = TRUE,
                 lcex = par("cex"),
                 lcol = par("col"),
                 lfont = par("font"),
                 lgperrow = NULL,
                 lgtitles = NULL,
                 lgw = 0.25,
                 lg.col =  "lightblue1",
                 lg.lwd = par("lwd"),
                 lty.axis = "solid",
                 lwd.axis = 1,
                 lwd.ticks.axis = lwd.axis,
                 main = NULL,
                 markerformatlist = flist,
                 maxnbrcolsfordups = 2,
                 pdf.bg = "transparent",
                 pdf.family = "Helvetica",
                 pdf.fg = "black",
                 pdf.width = NULL,
                 pdf.height = NULL,
                 pdf.pointsize = 12,
                 pdf.title = "LinkageMapView R output",
                 posonleft = NULL,
                 prtlgtitles = TRUE,
                 qtldf = NULL,
                 revthese = NULL,
                 rcex = par("cex"),
                 rcol = par("col"),
                 rfont = par("font"),
                 roundpos = 1,
                 rsegcol = TRUE,
                 ruler = FALSE,
                 sectcoldf = NULL,
                 segcol = NULL,
                 qtlscanone = NULL,
                 showonly = NULL,
                 units = "cM",
                 ylab = units
)

#For Chr10-Chr12
lmv.linkage.plot(input.data,
                 "ScnI-HYBRID-EFS-FS22-Chr10-12.pdf",
                 mapthese = c("LG10","LG11","LG12"),
                 at.axis = NULL,
                 autoconnadj = TRUE,
                 cex.axis = 1,
                 cex.lgtitle = par("cex.main"),
                 cex.main = par("cex.main"),
                 col.axis = par("col.axis"),
                 col.lgtitle = par("col.main"),
                 col.main = par("col.main"),
                 conndf = NULL,
                 denmap = FALSE,
                 dupnbr = FALSE,
                 font.axis = par("font.axis"),
                 font.lgtitle = par("font.main"),
                 font.main = par("font.main"),
                 header = TRUE,
                 labdist = 0.3,
                 labels.axis = TRUE,
                 lcex = par("cex"),
                 lcol = par("col"),
                 lfont = par("font"),
                 lgperrow = NULL,
                 lgtitles = NULL,
                 lgw = 0.25,
                 lg.col =  "lightblue1",
                 lg.lwd = par("lwd"),
                 lty.axis = "solid",
                 lwd.axis = 1,
                 lwd.ticks.axis = lwd.axis,
                 main = NULL,
                 markerformatlist = flist,
                 maxnbrcolsfordups = 2,
                 pdf.bg = "transparent",
                 pdf.family = "Helvetica",
                 pdf.fg = "black",
                 pdf.width = NULL,
                 pdf.height = NULL,
                 pdf.pointsize = 12,
                 pdf.title = "LinkageMapView R output",
                 posonleft = NULL,
                 prtlgtitles = TRUE,
                 qtldf = NULL,
                 revthese = NULL,
                 rcex = par("cex"),
                 rcol = par("col"),
                 rfont = par("font"),
                 roundpos = 1,
                 rsegcol = TRUE,
                 ruler = FALSE,
                 sectcoldf = NULL,
                 segcol = NULL,
                 qtlscanone = NULL,
                 showonly = NULL,
                 units = "cM",
                 ylab = units
)

