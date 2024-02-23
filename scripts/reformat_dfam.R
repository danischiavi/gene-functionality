args <- commandArgs(trailingOnly=TRUE)

# Read command line arguments
file1 <- args[1]
file2 <- args[2]
file3 <- args[3]

# Load data
df1 <- read.csv(file1, stringsAsFactors=FALSE)
df2 <- read.csv(file2, stringsAsFactors=FALSE)

# Process data
for(i in 1:nrow(df2)) {
    index1 <- grep(df2[i, 2], df1[, 4])                  
    index2 <- grep(df2[i, 1], df1[, 3])
      
    if (length(index1) > 0 & length(index2) > 0) {
        index <- intersect(index1, index2)
        
        if (length(index) > 0) {
        df1[index, 'Dfam_min'] <- df2[i, 4]
        df1[index, 'Dfam_sum'] <- df2[i, 5]
        }
    }
}

# Write processed data
write.csv(df1, file3, quote=FALSE, row.names=FALSE)