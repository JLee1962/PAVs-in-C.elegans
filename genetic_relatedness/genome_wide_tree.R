### This scrirpt was made by Daehan Lee.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(phangorn)
library(ape)
library(tibble)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}"))


## tree generation

# Generate phylogeny files

system(glue::glue("python vcf2phylip-master/vcf2phylip.py -i Ce328_complete_sites_16lrstrains.vcf.gz"))
WI_phy <- read.phyDat("Ce328_complete_sites_16lrstrains.min4.phy", format = "interleaved")
save(WI_phy, file = "WI_phy.RData")

load("Processed_Data/WI_phy.RData")

dm <- dist.ml(WI_phy)
treeNJ <- NJ(dm)

tree_NJ_rooted <- root(treeNJ, "XZ1516")

plot_tree <- ggtree(tree_NJ_rooted, branch.length="rate")
plot_tree + geom_tiplab() + geom_treescale(x=0.02, y=-1)

df_tree <- na.omit(plot_tree[[1]])

ggsave("16strain_tree.pdf", width=7.5, height=6)
