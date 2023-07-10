

library(stringr)
library(utils)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Input file path is missing.")
}

# Input and output file paths
input_file <- args[1]
output_file <- args[2]

# Function to convert read identifiers
convert_read_identifier <- function(read_identifier) {
  split_id <- strsplit(read_identifier, ".", fixed = TRUE)[[1]]
  converted_id <- paste0(split_id[1], ".", split_id[2], "/",split_id[3])
  return(converted_id)
}

# Open input and output files
con_in <- gzfile(input_file, "r")
con_out <- file(output_file, "w")

# Process input file line by line
while (length(line <- readLines(con_in, n = 1, warn = FALSE)) > 0) {
  if (startsWith(line, "@")) {
    # Read identifier line
    converted_id <- convert_read_identifier(line)
    writeLines(c(converted_id), con_out)
  } else {
    # Write other lines as they are
    writeLines(c(line), con_out)
  }
}

# Close input and output files
close(con_in)
close(con_out)
