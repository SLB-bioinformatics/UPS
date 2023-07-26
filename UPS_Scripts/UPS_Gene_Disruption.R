#UPS_Gene_Disruption.R
library("GenomicRanges")
args <- commandArgs(trailingOnly = TRUE)

# Define the input gene names
CNV_df    <- args[1]
Gene_path <- "/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/Support_Data/Ensembl_Gene_Chr_Start_End.txt"
Gene_df   <- read.csv(Gene_path)
SL_path <- "/mnt/beegfs/hassan/Stuart/UniformProcessingOfSamples/Support_Data/Human_SL_with_Ensembl_Ids.csv"
SL_df   <- read.csv(SL_path)


print("Genes")
print(Gene_df[1:10,])
print("CNV")
print(CNV_df[1,])

CNV_df$chromosome <- as.character(CNV_df$chromosome)
CNV_df$start      <- as.numeric(as.character(CNV_df$start))
CNV_df$end        <- as.numeric(as.character(CNV_df$end))

Gene_df$Chromosome.scaffold.name     <- as.character(Gene_df$Chromosome.scaffold.name)
Gene_df$Gene.start..bp.              <- as.numeric(as.character(Gene_df$Gene.start..bp.))
Gene_df$Gene.end..bp.                <- as.numeric(as.character(Gene_df$Gene.end..bp.))

CNV_Range <- GRanges(seqnames = CNV_df$chromosome,
                  ranges = IRanges(start = CNV_df$start, end = CNV_df$end))

Gene_Range <- GRanges(seqnames = paste("chr", Gene_df$Chromosome.scaffold.name, sep=""),
                  ranges = IRanges(start = Gene_df$Gene.start..bp., end = Gene_df$Gene.end..bp.))


print(CNV_Range)
print(Gene_Range)

CNV_Gene_overlaps <- findOverlaps(CNV_Range, Gene_Range)

# Print the overlaps
print(CNV_Gene_overlaps)

if(! file.exists("Gene_Disruption_By_CNV.csv")){
    CNV_Gene_overlaps_df <- as.data.frame(CNV_Gene_overlaps)
    Gene_Disruptuon_Df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 8))
    print("Assessing gene disruption by copy number")
    i <- 0
    for(curGenes in unique(CNV_Gene_overlaps_df$subjectHits)){
        i <- i + 1
        print(length(unique(CNV_Gene_overlaps_df$subjectHits)) - i)
        cur_OverLaps <- CNV_Gene_overlaps_df[curGenes == CNV_Gene_overlaps_df$subjectHits,]
        cur_Gene_df <- Gene_df[curGenes ,]
        cur_CNV_df  <- CNV_df[cur_OverLaps$queryHits ,]

        if(nrow(cur_OverLaps) == 1){
            newrow <- as.data.frame(matrix(data=c(cur_Gene_df$Gene.stable.ID,
                                                  cur_Gene_df$Gene.name,
                                                  cur_Gene_df$Chromosome.scaffold.name,
                                                  cur_Gene_df$Gene.start..bp.,
                                                  cur_Gene_df$Gene.end..bp.,
                                                  cur_CNV_df$cn,
                                                  cur_CNV_df$cn,
                                                  cur_CNV_df$cn), nrow=1))
            colnames(newrow) <- c("Gene_Id","Gene_name",
                                  "Chr", "Start", "End",
                                  "Average_CN","Max_CN","Min_CN")    
            Gene_Disruptuon_Df <- rbind(Gene_Disruptuon_Df,newrow)                            
        }else{
            cur_CNV_df$start[cur_CNV_df$start < cur_Gene_df$Gene.start..bp.] <- cur_Gene_df$Gene.start..bp.
            cur_CNV_df$end[cur_CNV_df$end > cur_Gene_df$Gene.end..bp.] <- cur_Gene_df$Gene.end..bp.
            cur_CNV_df$width <- cur_CNV_df$end - cur_CNV_df$start
            Gene_Size <- cur_Gene_df$Gene.end..bp. - cur_Gene_df$Gene.start..bp.
            cur_CNV_df$CN_waited_width <- cur_CNV_df$width  * cur_CNV_df$cn
            Average_CN  <- sum(cur_CNV_df$CN_waited_width) / Gene_Size 
            Max_CN      <- min(cur_CNV_df$cn)
            Min_CN      <- max(cur_CNV_df$cn)
            newrow <- as.data.frame(matrix(data=c(cur_Gene_df$Gene.stable.ID,
                                                  cur_Gene_df$Gene.name,
                                                  cur_Gene_df$Chromosome.scaffold.name,
                                                  cur_Gene_df$Gene.start..bp.,
                                                  cur_Gene_df$Gene.end..bp.,
                                                  Average_CN,
                                                  Max_CN,
                                                  Min_CN), nrow=1))
            colnames(newrow) <- c("Gene_Id","Gene_name",
                                  "Chr", "Start", "End",
                                  "Average_CN","Max_CN","Min_CN")    
            Gene_Disruptuon_Df <- rbind(Gene_Disruptuon_Df,newrow)     
        }  
    }

    print(Gene_Disruptuon_Df)
    write.csv(Gene_Disruptuon_Df, "Gene_Disruption_By_CNV.csv")
}

Gene_Disruptuon_Df <- read.csv("Gene_Disruption_By_CNV.csv")

Gene_Disruptuon_Df[1:20,]

SL_df[1:10,]

Gene_Disruptuon_Df$Average_CN <- as.numeric(as.character(Gene_Disruptuon_Df$Average_CN))
Genes_With_CN0 <- Gene_Disruptuon_Df[Gene_Disruptuon_Df$Average_CN == 0,]

i <- 0
Target_Df <- as.data.frame(matrix(data=NA, nrow=0, ncol=6))
colnames(Target_Df) <- c("Disrupted_Gene_Id","Synthetic_Lethal_Id",
                             "Disrupted_Gene_Name","Synthetic_Lethal_Name",
                             "r.source", "r.statistic_score")
for(curGene in unique(Genes_With_CN0$Gene_Id)){
    i <- i + 1
    print(length(unique(Genes_With_CN0$Gene_Id)) - i)
    SL_Gene1 <- SL_df[grep(curGene, SL_df$Gene_ID_1),]
    SL_Gene2 <- SL_df[grep(curGene, SL_df$Gene_ID_2),]
    if(nrow(SL_Gene1) > 0){
        newrow <- as.data.frame(matrix(data=c(SL_Gene1$Gene_ID_1, SL_Gene1$Gene_ID_2,
                                              SL_Gene1$n1.name, SL_Gene1$n2.name,
                                              SL_Gene1$r.source, SL_Gene1$r.statistic_score), 
                                              ncol=6))
        colnames(newrow) <- c("Disrupted_Gene_Id","Synthetic_Lethal_Id",
                             "Disrupted_Gene_Name","Synthetic_Lethal_Name",
                             "r.source", "r.statistic_score")
        Target_Df <- rbind(Target_Df,newrow)
    }

    if(nrow(SL_Gene2) > 0){
        newrow <- as.data.frame(matrix(data=c(SL_Gene2$Gene_ID_2, SL_Gene2$Gene_ID_1,
                                              SL_Gene2$n2.name,  SL_Gene2$n1.name,
                                              SL_Gene2$r.source, SL_Gene2$r.statistic_score), 
                                              ncol=6))
        colnames(newrow) <- c("Disrupted_Gene_Id","Synthetic_Lethal_Id",
                             "Disrupted_Gene_Name","Synthetic_Lethal_Name",
                             "r.source", "r.statistic_score")
        Target_Df <- rbind(Target_Df, newrow)
    }
}

write.csv(Target_Df, "CN_0_SytheticPredctions_SLDB.csv")