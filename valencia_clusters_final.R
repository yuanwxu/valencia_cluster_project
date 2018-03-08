# Analysis of multiple transmission clusters 
# Final version: search email for "identical seqs in valencia clusters"

library(ape)
library(tidyverse)
library(treedater)

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
