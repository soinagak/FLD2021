###########make scatter plots in Fig 1a, and Ext. Data Fig 1a,c ##############

Col_H3_1_coverage <- read.delim("path_to_folder/Col_H3_1_coverage.bed", header=FALSE, col.names = c("chr","start","end","ID","annotation","direction","WT_H3_reads","a","b","c"))

Col_H3_1_coverage$a <- NULL
Col_H3_1_coverage$b <- NULL
Col_H3_1_coverage$c <- NULL

fld_H3_1_coverage <- read.delim("path_to_folder/fld_H3_1_coverage.bed", header=FALSE, col.names = c("chr","start","end","ID","annotation","direction","fld_H3_1_reads","a","b","c"))
fld_H3_2_coverage <- read.delim("path_to_folder/fld_H3_2_coverage.bed", header=FALSE, col.names = c("chr","start","end","ID","annotation","direction","fld_H3_2_reads","a","b","c"))

Col_K4me1_1_coverage <- read.delim("path_to_folder/Col_K4me1_1_coverage.bed", header=FALSE, col.names = c("chr","start","end","ID","annotation","direction","WT_K4me1_reads","a","b","c"))
fld_K4me1_1_coverage <- read.delim("path_to_folder/fld_K4me1_1_coverage.bed", header=FALSE, col.names = c("chr","start","end","ID","annotation","direction","fld_K4me1_1_reads","a","b","c"))
fld_K4me1_2_coverage <- read.delim("path_to_folder/fld_K4me1_2_coverage.bed", header=FALSE, col.names = c("chr","start","end","ID","annotation","direction","fld_K4me1_2_reads","a","b","c"))

Col_K4me2_1_coverage <- read.delim("path_to_folder/Col_K4me2_1_coverage.bed", header=FALSE, col.names = c("chr","start","end","ID","annotation","direction","WT_K4me2_reads","a","b","c"))
fld_K4me2_1_coverage <- read.delim("path_to_folder/fld_K4me2_1_coverage.bed", header=FALSE, col.names = c("chr","start","end","ID","annotation","direction","fld_K4me2_1_reads","a","b","c"))
fld_K4me2_2_coverage <- read.delim("path_to_folder/fld_K4me2_2_coverage.bed", header=FALSE, col.names = c("chr","start","end","ID","annotation","direction","fld_K4me2_2_reads","a","b","c"))

All_reads <- Col_H3_1_coverage
All_reads$fld_H3_1_reads <- fld_H3_1_coverage$fld_H3_1_reads
All_reads$fld_H3_2_reads <- fld_H3_2_coverage$fld_H3_2_reads

All_reads$WT_K4me1_reads <- Col_K4me1_1_coverage$WT_K4me1_reads
All_reads$fld_K4me1_1_reads <- fld_K4me1_1_coverage$fld_K4me1_1_reads
All_reads$fld_K4me1_2_reads <- fld_K4me1_2_coverage$fld_K4me1_2_reads

All_reads$WT_K4me2_reads <- Col_K4me2_1_coverage$WT_K4me2_reads
All_reads$fld_K4me2_1_reads <- fld_K4me2_1_coverage$fld_K4me2_1_reads
All_reads$fld_K4me2_2_reads <- fld_K4me2_2_coverage$fld_K4me2_2_reads

All_reads$WT_H3_RPM <- All_reads$WT_H3_reads/18.825640
All_reads$fld_H3_1_RPM <- All_reads$fld_H3_1_reads/15.818767
All_reads$fld_H3_2_RPM <- All_reads$fld_H3_2_reads/19.743565

All_reads$WT_K4me1_RPM <- All_reads$WT_K4me1_reads/21.387638
All_reads$fld_K4me1_1_RPM <- All_reads$fld_K4me1_1_reads/22.678456
All_reads$fld_K4me1_2_RPM <- All_reads$fld_K4me1_2_reads/24.943972

All_reads$WT_K4me2_RPM <- All_reads$WT_K4me2_reads/20.684162
All_reads$fld_K4me2_1_RPM <- All_reads$fld_K4me2_1_reads/16.943065
All_reads$fld_K4me2_2_RPM <- All_reads$fld_K4me2_2_reads/15.393192

Macs2_1659 <- read.delim("path_to_finder/macs2_fld_K4me1_increase_replicates_1659.txt", header=FALSE, col.names = c("ID"))

keywords <- Macs2_1659$ID
param <- 4
obj <- is.element(as.character(All_reads[,param]), keywords) 
All_reads$Macs2_1659 <- obj                                                

All_reads$WT_fld_K4me1_diff <- (All_reads$fld_K4me1_1_RPM + All_reads$fld_K4me1_2_RPM)/2 - All_reads$WT_K4me1_RPM


##############Scatter plots#####################
plot(sqrt(All_reads$WT_H3_RPM),(sqrt(All_reads$fld_H3_1_RPM)+sqrt(All_reads$fld_H3_2_RPM))/2, pch=16, cex=1, col=rgb(0,0,0,0.1),xlim=c(0,15),ylim=c(0,15),xlab="",ylab="")
par(new=T)
plot(sqrt(All_reads$WT_H3_RPM)[All_reads$Macs2_1659 == "TRUE" & All_reads$WT_K4me1_reads>0],(sqrt(All_reads$fld_H3_1_RPM)[All_reads$Macs2_1659 == "TRUE" & All_reads$WT_K4me1_reads>0]+sqrt(All_reads$fld_H3_2_RPM)[All_reads$Macs2_1659 == "TRUE" & All_reads$WT_K4me1_reads>0])/2, pch=16, cex=1, col=rgb(1,0,1,0.2),xlim=c(0,15),ylim=c(0,15),xlab="",ylab="")
par(new=T)
plot(sqrt(All_reads$WT_H3_RPM)[All_reads$ID == "AT5G10140"],(sqrt(All_reads$fld_H3_1_RPM)[All_reads$ID == "AT5G10140"] + sqrt(All_reads$fld_H3_2_RPM)[All_reads$ID == "AT5G10140"])/2, pch=16, cex=1, col=rgb(1,0,0,1),xlim=c(0,15),ylim=c(0,15),xlab="", ylab="")
par(new=T)
plot(sqrt(All_reads$WT_H3_RPM)[All_reads$ID == "AT3G10390"],(sqrt(All_reads$fld_H3_1_RPM)[All_reads$ID == "AT3G10390"] + sqrt(All_reads$fld_H3_2_RPM)[All_reads$ID == "AT3G10390"])/2, pch=16, cex=1, col=rgb(0,0,1,1),xlim=c(0,15),ylim=c(0,15),xlab="H3H3 in WT", ylab="H3H3 in fld")

plot(sqrt(All_reads$WT_K4me1_RPM[All_reads$Macs2_1659 == "FALSE"]),(sqrt(All_reads$fld_K4me1_1_RPM[All_reads$Macs2_1659 == "FALSE"])+sqrt(All_reads$fld_K4me1_2_RPM[All_reads$Macs2_1659 == "FALSE"]))/2, pch=16, cex=1, col=rgb(0,0,0,0.1),xlim=c(0,17),ylim=c(0,17),xlab="",ylab="")
par(new=T)
plot(sqrt(All_reads$WT_K4me1_RPM)[All_reads$Macs2_1659 == "TRUE" & All_reads$WT_K4me1_reads>0],(sqrt(All_reads$fld_K4me1_1_RPM)[All_reads$Macs2_1659 == "TRUE" & All_reads$WT_K4me1_reads>0]+sqrt(All_reads$fld_K4me1_2_RPM)[All_reads$Macs2_1659 == "TRUE" & All_reads$WT_K4me1_reads>0])/2, pch=16, cex=1, col=rgb(1,0,1,0.2),xlim=c(0,17),ylim=c(0,17),xlab="",ylab="")
par(new=T)
plot(sqrt(All_reads$WT_K4me1_RPM)[All_reads$ID == "AT5G10140"],(sqrt(All_reads$fld_K4me1_1_RPM)[All_reads$ID == "AT5G10140"] + sqrt(All_reads$fld_K4me1_2_RPM)[All_reads$ID == "AT5G10140"])/2, pch=16, cex=1, col=rgb(1,0,0,1),xlim=c(0,17),ylim=c(0,17),xlab="", ylab="")
par(new=T)
plot(sqrt(All_reads$WT_K4me1_RPM)[All_reads$ID == "AT3G10390"],(sqrt(All_reads$fld_K4me1_1_RPM)[All_reads$ID == "AT3G10390"] + sqrt(All_reads$fld_K4me1_2_RPM)[All_reads$ID == "AT3G10390"])/2, pch=16, cex=1, col=rgb(0,0,1,1),xlim=c(0,17),ylim=c(0,17),xlab="H3K4me1 in WT", ylab="H3K4me1 in fld")

plot(sqrt(All_reads$WT_K4me2_RPM),(sqrt(All_reads$fld_K4me2_1_RPM)+sqrt(All_reads$fld_K4me2_2_RPM))/2, pch=16, cex=1, col=rgb(0,0,0,0.1),xlim=c(0,15),ylim=c(0,15),xlab="",ylab="")
par(new=T)
plot(sqrt(All_reads$WT_K4me2_RPM)[All_reads$Macs2_1659 == "TRUE" & All_reads$WT_K4me1_reads>0],(sqrt(All_reads$fld_K4me2_1_RPM)[All_reads$Macs2_1659 == "TRUE" & All_reads$WT_K4me1_reads>0]+sqrt(All_reads$fld_K4me2_2_RPM)[All_reads$Macs2_1659 == "TRUE" & All_reads$WT_K4me1_reads>0])/2, pch=16, cex=1, col=rgb(1,0,1,0.2),xlim=c(0,15),ylim=c(0,15),xlab="",ylab="")
par(new=T)
plot(sqrt(All_reads$WT_K4me2_RPM)[All_reads$ID == "AT5G10140"],(sqrt(All_reads$fld_K4me2_1_RPM)[All_reads$ID == "AT5G10140"] + sqrt(All_reads$fld_K4me2_2_RPM)[All_reads$ID == "AT5G10140"])/2, pch=16, cex=1, col=rgb(1,0,0,1),xlim=c(0,15),ylim=c(0,15),xlab="", ylab="")
par(new=T)
plot(sqrt(All_reads$WT_K4me2_RPM)[All_reads$ID == "AT3G10390"],(sqrt(All_reads$fld_K4me2_1_RPM)[All_reads$ID == "AT3G10390"] + sqrt(All_reads$fld_K4me2_2_RPM)[All_reads$ID == "AT3G10390"])/2, pch=16, cex=1, col=rgb(0,0,1,1),xlim=c(0,15),ylim=c(0,15),xlab="K4me2K4me2 in WT", ylab="K4me2K4me2 in fld")

###############Comparison of biological replicates##################

plot(sqrt(All_reads$fld_K4me1_1_RPM[All_reads$Macs2_1659 == "FALSE"])-sqrt(All_reads$WT_K4me1_RPM[All_reads$Macs2_1659 == "FALSE"]),sqrt(All_reads$fld_K4me1_2_RPM[All_reads$Macs2_1659 == "FALSE"])-sqrt(All_reads$WT_K4me1_RPM[All_reads$Macs2_1659 == "FALSE"]), pch=16, cex=1, col=rgb(0,0,0,0.1),xlim=c(-1.5,3.5),ylim=c(-1.5,3.5),xlab="",ylab="")
par(new=T)
plot(sqrt(All_reads$fld_K4me1_1_RPM[All_reads$Macs2_1659 == "TRUE"])-sqrt(All_reads$WT_K4me1_RPM[All_reads$Macs2_1659 == "TRUE"]),sqrt(All_reads$fld_K4me1_2_RPM[All_reads$Macs2_1659 == "TRUE"])-sqrt(All_reads$WT_K4me1_RPM[All_reads$Macs2_1659 == "TRUE"]), pch=16, cex=1, col=rgb(1,0,1,0.2),xlim=c(-1.5,3.5),ylim=c(-1.5,3.5),xlab="",ylab="")
par(new=T)
plot(sqrt(All_reads$fld_K4me1_1_RPM[All_reads$ID == "AT5G10140"])-sqrt(All_reads$WT_K4me1_RPM[All_reads$ID == "AT5G10140"]),sqrt(All_reads$fld_K4me1_2_RPM[All_reads$ID == "AT5G10140"])-sqrt(All_reads$WT_K4me1_RPM[All_reads$ID == "AT5G10140"]), pch=16, cex=1, col=rgb(1,0,0,1),xlim=c(-1.5,3.5),ylim=c(-1.5,3.5),xlab="", ylab="")

cor(x=sqrt(All_reads$fld_K4me1_1_RPM)-sqrt(All_reads$WT_K4me1_RPM),y=sqrt(All_reads$fld_K4me1_2_RPM)-sqrt(All_reads$WT_K4me1_RPM),method="pearson")
#0.8580761