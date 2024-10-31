library(ape)

# reduce_tree.R
args <- commandArgs(trailingOnly = TRUE)

# load your tree
tree <- read.tree(args[1])

# we are removing the refseqs and the outgroup from the tree 
# so that parnas only selects reps from the rest of the tree 
# since the refseqs and outgroup are automatically included in the tree
tree_modified <- drop.tip(tree, c(args[3], args[4], args[5]))

# save the tree to a new file 
write.tree(tree_modified, file=args[2])