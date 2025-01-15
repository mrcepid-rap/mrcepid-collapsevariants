#!/usr/bin/Rscript --vanilla

# Load required libraries
library(data.table)
library(Matrix)

# Read in inputs. They are:
# 1. [1] Samples file – the rows of our final matrix
# 2. [2] Our list of variants – the columns of our final matrix
# 3. [3] matrix file – The formated matrix file for readMM
# 4. [4] Name of the .rds matrix file to save to
args <-  commandArgs(trailingOnly = T)
samples_file <- args[1]
variants_file <- args[2]
matrix_file <- args[3]
out_file <- args[4]

# Process the samples file:
samples <- fread(samples_file)
samples[,sampID:=as.character(sampID)] # Make sure to set as a character - typically gets read as a integer

# Process the variants file:
variants <- fread(variants_file)

# Read in the formatted matrix file:
gt_matrix <- readMM(matrix_file)

# Set row and column names according the the samples/variants files:
rownames(gt_matrix) <- samples[,sampID]
colnames(gt_matrix) <- variants[,varID]

# And save the final file to our AWS instance.
saveRDS(gt_matrix, out_file)

