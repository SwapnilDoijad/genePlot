library(optparse)
library(ape)
library(phangorn)

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input matrix file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", help = "Output Newick file", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if input and output files are provided
if (is.null(opt$input) || is.null(opt$output)) {
  stop("Arguments -i (input) and -o (output) are required.")
}

# Read the input matrix
all_distance <- read.table(opt$input)

# Generate the UPGMA tree
treeUPGMA <- upgma(all_distance)
mytree <- as.phylo(treeUPGMA)

# Write the tree to the output file
write.tree(mytree, file = opt$output)

# Plot the tree
plot(treeUPGMA, main = "UPGMA")