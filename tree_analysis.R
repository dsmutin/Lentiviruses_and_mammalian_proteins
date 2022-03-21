library(ape)

compare_trees <- function (tree1, tree2) {
  nwk1 <- read.tree(tree1)
  nwk2 <- read.tree(tree2)
  comparePhylo(tree1, tree2, plot = T, type = "tidy", cex = 1.2, location = NA,)
}

setwd ("~/CCR5_speedrun")
nwks <- c ("CXCR4.nwk", "CD4.nwk", "CCR5.nwk")

CCR5 <- read.tree("CCR5.nwk")
CD4 <- read.tree("CD4.nwk")
CXCR4 <- read.tree("CXCR4.nwk")

plot (chronoMPL(CCR5))

comparePhylo(CCR5, CD4, plot = T, type = "tidy", cex = 0.6, location = NA)
comparePhylo(CCR5, CXCR4, plot = T, type = "tidy", cex = 0.6, location = NA)
comparePhylo(CXCR4, CD4, plot = T, type = "tidy", cex = 0.6, location = NA)

r_CCR5 <- root(chronos(CCR5), "Chrysochloris_asiatica")
r_CD4 <- root(chronos(CD4), "Chrysochloris_asiatica")
r_CXCR4 <- root(chronoMPL(CXCR4), "Chrysochloris_asiatica")
dndlist <- dendextend::dendlist(r_CCR5, r_CD4)

