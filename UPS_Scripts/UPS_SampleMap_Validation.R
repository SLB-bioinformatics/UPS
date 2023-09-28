
args <- commandArgs(trailingOnly = TRUE)

suppressWarnings(df <- read.csv(args[1]))

SMcols <- colnames(df)

required_Cols <- c("Sample_ID","Cohort","Patient_ID",
                   "Cancer_Type","Sequencing_Type","Data_Format",
                   "Source","Sex","Location_Of_Data")

# Validate ColNames
OutputDf <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 2))
for(CurCol in required_Cols){
    CurColPos <- grep(CurCol,SMcols)
    if(length(CurColPos)>0){
        newrow <- as.data.frame(matrix(data = c(CurCol,CurColPos), nrow = 1, ncol = 2))
    }else{
        CurColPos <- "Column Missing"
        newrow <- as.data.frame(matrix(data = c(CurCol,CurColPos), nrow = 1, ncol = 2))

    }
    colnames(newrow) <- c("Column","Found_in_Col")
    OutputDf <- rbind(OutputDf, newrow)
}

# Validate Data types
OutputDf$Validate <- NA
for(i in 1:nrow(OutputDf)){
    curCol <- OutputDf$Column[i]
    CurLoc <- OutputDf$Found_in_Col[i]
    if(CurLoc != "Column Missing"){

        if(curCol == "Sample_ID"){
            nSamples    <- nrow(df)
            N_Sample_ID <- length(unique(df$Sample_ID))
            if(nSamples == N_Sample_ID){
                report <- "Valid"
            } else {
                report <- "Invalid: All Sample Id Must Be Unique"
            }
        } else if(curCol == "Sequencing_Type"){
            valid_values <- c("WGS", "WES")
            valid <- all(grepl(paste(valid_values, collapse = "|"), df$Sequencing_Type))
            if(valid == TRUE){
                report <- "Valid"
            } else {
                report <- "Invalid: Sequencing_Type only accepts WGS or WES"
            }
        } else if(curCol == "Source"){
            valid_values <- c("Cell_Line", "Patient")
            valid <- all(grepl(paste(valid_values, collapse = "|"), df$Source))
            if(valid == TRUE){
                report <- "Valid"
            } else {
                report <- "Invalid: Source only accepts Cell_Line or Patient"
            }
        } else if(curCol == "Location_Of_Data"){
            valid <- all(file.exists(df$Location_Of_Data))
            if(valid == TRUE){
                report <- "Valid"
            }else{
                report <- "Invalid: Not All Files found for Location_Of_Data"
            }
        } else{
            report <- "Free Text"
        }
    OutputDf$Validate[i] <- report

    } else {
        report <- "Invalid: Column Missing"
    }
}
write.csv(OutputDf, "Sample_Map_Validation_Report.csv")
