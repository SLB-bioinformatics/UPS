# UPS_Plotting.R
library("ggplot2")
Path_to_CNV <- "Copy_Number_Variation/CNVkit/"
Path_to_SV <- "Structural_Variant/Lumpy/"
Path_to_SNV <- "Single_Nucleotide_Variant/Mutect/"
Path_to_Support_Data <- "/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/Support_Data/"

GenomicPositions <- "/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/Support_Data/HG38_Genome_Location_For_Plotting.csv"
HumanChromosomePositions <- read.csv(GenomicPositions)


CNV_files <- list.files(Path_to_CNV)
print(CNV_files )
CNV_df      <- read.delim(paste(Path_to_CNV,CNV_files[grep(".call.cns",CNV_files )], sep=""))


CNV_df$chromosome <- gsub("chr", "Chromosome ",CNV_df$chromosome)

print(CNV_df[1,])
CNV_df$GenomicPositionStart <- 0
CNV_df$GenomicPositionEnd <- 0

for(curChr in HumanChromosomePositions$ChromosomeName){
  curChromosomeStart <- as.numeric(HumanChromosomePositions$Start[curChr == HumanChromosomePositions$ChromosomeName])
  CNV_df$GenomicPositionStart[curChr == CNV_df$chromosome] <- as.numeric(CNV_df$start[curChr == CNV_df$chromosome]) + curChromosomeStart
  CNV_df$GenomicPositionEnd[curChr == CNV_df$chromosome] <- as.numeric(CNV_df$end[curChr == CNV_df$chromosome]) + curChromosomeStart
}

pdf("Reports/Plots/Genome_And_Chromosomes.pdf", width=20)

copynumber_Plot <- ggplot(data=CNV_df, x= end) +
                     geom_rect(data = HumanChromosomePositions, aes(xmin=Start,xmax=End, ymin=-1, ymax=0 )) +
                     geom_rect(data = CNV_df, aes(xmin=GenomicPositionStart,xmax=GenomicPositionEnd, ymin=0, ymax=cn ), fill="#002147") 

plot(copynumber_Plot)

cytobandDf <- read.csv("/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/Support_Data/cytoBand.csv")
print(cytobandDf[1:10,])
print(unique(cytobandDf$Colour))
for(curChr in unique(CNV_df$chromosome)){    
    curcytobandDf <- cytobandDf[curChr == cytobandDf$Chromosome,]
    curcytobandDf$row_num <- seq_len(nrow(curcytobandDf))
    curCNV_df <- CNV_df[curChr == CNV_df$chromosome,]
    copynumber_Plot <- ggplot(data=curCNV_df, x= end) +
                         geom_rect(data = curcytobandDf, aes(xmin=Start,xmax=End, ymin=-1, ymax=0), fill=curcytobandDf$Colour) +
                         geom_text(data = curcytobandDf, aes(x = (Start + End) / 2, y = -0.5, 
                                   label = Cytoband, angle = 45),
                                   color = ifelse(curcytobandDf$Colour %in% c("#000000", "#333333", "#666666"), "white", "black"),
                                   size = 3) +
                         geom_rect(data = curCNV_df, aes(xmin=start,xmax=end, ymin=0, ymax=cn ), fill="#002147") +
                         scale_fill_identity(guide = "none") +
                         scale_y_continuous(breaks = 0:max(curCNV_df$cn)) +
                         ggtitle(curChr)
    plot(copynumber_Plot)
}

