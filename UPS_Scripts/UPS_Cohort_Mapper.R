library("ggplot2")
library("reshape2")
library("dplyr")
library("stringr")


pdf("CohortMap.pdf")
SampleMap <- "../UPS_Demo_Data/SampleMapOSShort.csv"
SampleMapDf <- read.csv(SampleMap)
SampleMapDf$Cohort <- gsub("_"," ",SampleMapDf$Cohort)
SampleMapDf$Cancer_Type <- gsub("_"," ",SampleMapDf$Cancer_Type)

NumberOfPatients <- length(unique(SampleMapDf$Patient_ID))
NumberOfPatients

NumberOfCancerTypes <- length(unique(SampleMapDf$Cancer_Type))
NumberOfCancerTypes

NumberOfCohorts <- length(unique(SampleMapDf$Cohort))
NumberOfCohorts 

counts <- SampleMapDf %>%
  group_by(Cancer_Type) %>%
  summarise(Unique_Samples = n_distinct(Sample_ID))
counts

ggplot(counts, aes(x = "", y = Unique_Samples, fill = str_wrap(Cancer_Type, width = 25))) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  labs(title = "Samples in Each Cancer Type") +
  scale_fill_discrete(name = "Cohort") +
  theme_void()+
  theme(
    legend.spacing.x = unit(0.5, "cm"),  
    legend.spacing.y = unit(0.5, "cm")   
  )

counts <- SampleMapDf %>%
  group_by(Cohort) %>%
  summarise(Unique_Samples = n_distinct(Sample_ID))

ggplot(counts, aes(x = "", y = Unique_Samples, fill = str_wrap(Cohort, width = 25))) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  labs(title = "Samples in Each Cohort") +
  scale_fill_discrete(name = "Cohort") +
  theme_void()+
  theme(
    legend.spacing.x = unit(0.5, "cm"),  
    legend.spacing.y = unit(0.5, "cm")   
  )


outputdf <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
for(curSample in SampleMapDf$Patient_ID){
    curdf <- SampleMapDf[grep(curSample, SampleMapDf$Patient_ID),]
    WGSdf <- curdf[grep("WGS", curdf$Sequencing_Type ),]
    WESdf <- curdf[grep("WES", curdf$Sequencing_Type ),]
    RNAdf <- curdf[grep("RNA", curdf$Sequencing_Type ),]
    Tdf   <- curdf[grep("Yes", curdf$Is_Tumor),]
    Ndf   <- curdf[grep("No", curdf$Is_Tumor),]
    if(nrow(WGSdf) > 0 && nrow(RNAdf) > 0 ){
        seqeunincing <- "WGS_And_RNA"
    } else if(nrow(WESdf) > 0 && nrow(RNAdf) > 0 ){
        seqeunincing <- "WES_And_RNA"  
    } else if(nrow(WGSdf) > 0 ){
        seqeunincing <- "WGS"
    } else if(nrow(WESdf) > 0 ){
        seqeunincing <- "WES"
    } else if(nrow(RNAdf) > 0 ){
        seqeunincing <- "RNA"
    }

    if(nrow(Tdf) > 0 && nrow(Ndf) > 0 ){
        Samples <- "TumourNormal"
    } else if(nrow(Tdf) > 0 ){
        Samples <- "Tumour"
    } else if(nrow(Ndf) > 0 ){
        Samples <- "Normal"
    }
    newrow <- as.data.frame(matrix(data = c(curSample, seqeunincing, Samples), nrow = 1, ncol = 3))
    colnames(newrow) <- c("Sample_Name","Seqeuncing_Type","Sample_Type")
    outputdf <- rbind(outputdf, newrow)
}

outputdf[c(1,100,50,220,300),]

outputdf$Seqeuncing_Type <- as.character(outputdf$Seqeuncing_Type)
outputdf$Sample_Type <- as.character(outputdf$Sample_Type)

outputdf$Seqeuncing_Type <- as.character(outputdf$Seqeuncing_Type)
outputdf$Sample_Type <- as.character(outputdf$Sample_Type)

table_data <- table(outputdf$Seqeuncing_Type, outputdf$Sample_Type)
df_heatmap <- as.data.frame.matrix(table_data)
print(df_heatmap)

df <- data.frame(Method = rownames(df_heatmap))
for (i in 1:ncol(df_heatmap)) {
  df[[colnames(df_heatmap)[i]]] <- df_heatmap[, i]
}

df_long <- melt(df, id.vars = "Method")

ggplot(df_long, aes(x = Method, y = variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), vjust = 1) +
  scale_fill_gradient(low = "#ADD8E6", high = "#228B22") +
  theme_minimal() +
  labs(
    title = "Sequencing Method and Sample Type",
    x = "Method",
    y = "Sample Type",
    fill = "Value")+
  guides(fill = "none")

