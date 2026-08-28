#!/usr/bin/env Rscript
# export_empirical_networks.R
#
# Export the six empirical networks used in the paper (Bat, Highschool boy,
# House finch, Surfer, Hall, Jazz; main.tex Appendix A) from the `networkdata`
# R package (https://github.com/schochastics/networkdata) to edge-list CSVs
# consumed by core/kuramoto_core.jl's `load_empirical`.
#
# Usage:
#   Rscript export_empirical_networks.R
#
# Requires: install.packages("igraph");
#           remotes::install_github("schochastics/networkdata")

suppressPackageStartupMessages({
  library(networkdata)
  library(igraph)
})

outdir <- "nets"  # run this script from the networks/ directory
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# outname (-> nets/<outname>_edges.csv), dataset name in `networkdata`,
# 1-based index into the dataset if it is a list of graphs (NA otherwise).
spec <- list(
  c("animal_18_1",     "animal_18",       1),   # Bat            (43 nodes, 546 edges)
  c("highschool_boys",  "highschool_boys", NA),  # Highschool boy (70, 274)
  c("animal_6",         "animal_6",        NA),  # House finch    (108, 1026)
  c("surfersb",         "surfersb",        NA),  # Surfer         (43, 336)
  c("hall",             "hall",            NA),  # Hall           (217, 1839)
  c("jazz",             "jazz",            NA)   # Jazz           (198, 2742)
)

get_graph <- function(dataset, index) {
  e <- new.env()
  data(list = dataset, package = "networkdata", envir = e)
  g <- get(dataset, envir = e)
  if (!igraph::is_igraph(g) && is.list(g)) {
    idx <- if (is.na(index)) 1 else as.integer(index)
    g <- g[[idx]]
  }
  g
}

for (s in spec) {
  outname <- s[1]; dataset <- s[2]; index <- s[3]
  g <- get_graph(dataset, index)
  g <- igraph::simplify(igraph::as_undirected(g, mode = "collapse"))
  el <- igraph::as_edgelist(g, names = TRUE)
  write.csv(data.frame(from = el[, 1], to = el[, 2]),
            file.path(outdir, paste0(outname, "_edges.csv")), row.names = FALSE)
  cat(sprintf("%-16s <- %-16s[%s]  %d nodes, %d edges (largest connected component is taken later, in Julia)\n",
              outname, dataset, ifelse(is.na(index), "-", index), gorder(g), gsize(g)))
}
