# Analysis of multiple transmission clusters 
# Final version: search email for "identical seqs in valencia clusters"

library(ape)
library(tidyverse)
library(treedater)
library(TransPhylo)

# Function to read ML trees corresponding to the clusters
# Here these are built from raxml
# Return --- list of trees
read_bestTree <- function(path = ".", pattern){
  f <- list.files(path, pattern, full.names = TRUE, recursive = TRUE)
  purrr::map(f, ape::read.tree)
}

cls_mltree <- read_bestTree(path = "~/Biomath/raxml_ng/output_v2", pattern = "bestTree")
cls_name <- list.dirs("~/Biomath/raxml_ng/output_v2", full.names = FALSE, recursive = FALSE)

# Unroot trees so that treedater runs
cls_mltree_unrooted <- map(cls_mltree, unroot)

# Scale edge lengths as per CC's email
cls_mltree_unrooted <- map(cls_mltree_unrooted, function(x) {x$edge.length <- x$edge.length * 4400; x})

# Copy and modify get_sts2Treedater from analysis4.Rmd so it works with the new version of 
# the (single) diagnosis date file containing dates of all cluster tips
# tre --- ML tree
# clname --- cluster name of tre
get_sts2Treedater <- function(tre, clname, epi_file = "Data/new_v2/valencia_transphylo_clusters_diagnostic_dates.xlsx"){
  epi <- readxl::read_excel(epi_file)
  sts <- epi %>%
    filter(TR_CL == clname) %>%
    select(ID_old, Diagnostic_date) %>% 
    mutate(Diagnostic_date = lubridate::decimal_date(Diagnostic_date))
  
  tips_no_match <- dplyr::setdiff(tre$tip.label, sts$ID_old)
  if(length(tips_no_match) > 0) warning("Found tips on the tree that are missing in epi data")
  
  num_tips <- length(tre$tip.label)
  sts2 <- set_names(rep(NA, num_tips), tre$tip.label)
  sts2[sts$ID_old] <- sts$Diagnostic_date
  
  if(length(tips_no_match) != sum(is.na(sts2)))
    warning("There is record in epi data but the date field is missing.")
  
  if(length(tips_no_match) > 0)
    return(list(sts=sts2, tips_no_match=tips_no_match))
  else
    return(list(sts=sts2))
}

cls_sts <- map2(cls_mltree_unrooted, cls_name, get_sts2Treedater)

# Make time-trees from treedater
# Source tdaterAnalysis from analysis4.Rmd
cls_timetree <- tdaterAnalysis(cls_mltree_unrooted, cls_sts, 1, strictClock = TRUE, meanRateLimits = c(0.3, 0.7))

# Fix zero branch lengths of timed trees by shifting *all* branch lengths by 0.1
fix_edge_length <- function(timed_tree){
  if(any(dplyr::near(timed_tree$edge.length,0))){
    timed_tree$edge.length <- timed_tree$edge.length + 0.1
    return(timed_tree)
  }
  else # do nothing
    return(timed_tree)
}
cls_timetree <- map(cls_timetree, fix_edge_length)

# Run TransPhylo
cls_lastdate <- map(cls_timetree, ~ max(.$sts))
cls_ptree <- map2(cls_timetree, cls_lastdate, ptreeFromPhylo)
source("~/Biomath/TransPhylo/R/proposal.R")
source("~/Biomath/TransPhylo/R/computeHost.R")
source("~/Biomath/TransPhylo/R/infer_multiTTree_shareParam.R")

ws.shape <- 1.1; ws.scale <- 1/0.4
w.shape <- 1.3; w.scale <- 1/0.3
iters <- 2e4; thin <- 10
set.seed(1)
cls_record <- infer_multiTTree_shareParam(cls_ptree, w.shape, w.scale, ws.shape, ws.scale, 
                                         mcmcIterations = iters, thinning = thin, 
                                         share = c("neg","off.r","off.p","pi"))
# Simple trace plot of shared parameters
plot_trace <- function(record){
  par(mfrow = c(1,3))
  plot(map_dbl(record, "neg"), type = "l", ylab = "neg", xlab = "MCMC step")
  plot(map_dbl(record, "off.r"), type = "l", ylab = "off.r", xlab = "MCMC step")
  plot(map_dbl(record, "pi"), type = "l", ylab = "pi", xlab = "MCMC step")  
}
plot_trace(cls_record[[1]])

discard_burnin <- function(record_lst, burnin = 0.2){
  map(record_lst, function(x) x[round(length(x)*burnin):length(x)])
}
cls_record_small <- discard_burnin(cls_record)
rm(cls_record)

# Share params within two groups defined by edge length scale
g1 <- c(2, 5, 10, 16, 17, 21) # scale = 0.5
g2 <- setdiff(1:20, g1) # scale >= 1
cls_record_g1 <- infer_multiTTree_shareParam(cls_ptree[g1], w.shape, w.scale, ws.shape, ws.scale, 
                                          mcmcIterations = iters, thinning = thin, 
                                          share = c("neg","off.r","off.p","pi"))
cls_record_g2 <- infer_multiTTree_shareParam(cls_ptree[g2], w.shape, w.scale, ws.shape, ws.scale, 
                                          mcmcIterations = iters, thinning = thin, 
                                          share = c("neg","off.r","off.p","pi"))
names(cls_record_g1) <- cls_name[g1]
names(cls_record_g2) <- cls_name[g2] # for easy referencing
plot_trace(cls_record_g1[[1]])
plot_trace(cls_record_g2[[1]])

## Observation: g1 shows higher sampling rate pi than g2, and lower neg; whereas sharing params for all clusters 
#               averages out these effects

# Check transmission tree for cluster 77
source("~/Biomath/TransPhylo/R/selectTTree/selectTTree.R")
whichtree <- selectTTree(cls_record_g2[["CL077"]], burnin = 0.2)
plotCTree(cls_record_g2[["CL077"]][[whichtree]]$ctree)
plotTTree2(extractTTree(cls_record_g2[["CL077"]][[whichtree]]$ctree))
