library(ggplot2)
library(ape)
library(phangorn)
library(ggtree)

# This script calculates boostraps from the outputs of the python script "hammingDistanceTree.py"
# https://www.biostars.org/p/105010/
# Repurposed from boot.phylo: https://github.com/cran/ape/blob/0c004ab30b16263e339dd3e4715b8c14bb631fca/R/dist.topo.R

calc.boot.phylo <-
  function(phy, boot, trees = FALSE, rooted = is.rooted(phy))
  {
      boot.tree <- .compressTipLabel(boot, ref = phy$tip.label)
      boot.tree <- .uncompressTipLabel(boot.tree)
      boot.tree <- unclass(boot.tree) # otherwise countBipartitions crashes
    if (rooted) {
      pp <- prop.part(boot.tree)
      ans <- prop.clades(phy, part = pp, rooted = rooted)
    } else {
      phy <- reorder(phy, "postorder")
      ints <- phy$edge[, 2] > Ntip(phy)
      ans <- countBipartitions(phy, boot.tree)
      ans <- c(B, ans[order(phy$edge[ints, 2])])
    }
    if (trees) {
      class(boot.tree) <- "multiPhylo"
      ans <- list(BP = ans, trees = boot.tree)
    }
    ans
  }

# read in neighbor-joining tree (single phylo)
njtree <- read.tree("")

# read in permuted neighbor-joining trees (multi phylo)
boots <- read.tree("")

# calculate bootstraps from neighbor-joining trees
bootstraps <- calc.boot.phylo(phy = njtree, boot = boots)

# add boostraps to neighbor joining tree
njtree$node.labels <- bootstraps

mpt <- midpoint(njtree)
ggtree(mpt, ladderize = T, alpha = .5) + geom_nodelab(hjust = .5)
