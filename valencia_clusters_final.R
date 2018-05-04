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

# Plot the representative combined tree for all clusters and store in a file
pdf("cls_ctree_g1.pdf")
for(nm in names(cls_record_g1)){
  whichtree <- selectTTree(cls_record_g1[[nm]], burnin = 0.2)
  plotCTree(cls_record_g1[[nm]][[whichtree]]$ctree)
  title(main = nm)
}
dev.off()
pdf("cls_ctree_g2.pdf")
for(nm in names(cls_record_g2)){
  whichtree <- selectTTree(cls_record_g2[[nm]], burnin = 0.2)
  plotCTree(cls_record_g2[[nm]][[whichtree]]$ctree)
  title(main = nm)
}
dev.off()

# Plot ctree of CL002 with visNetwork
# source networkPlot() from PHE_TB/transphylo.R
whichtree <- selectTTree(cls_record_g1[["CL002"]], burnin = 0.2)
best_tree_cl002 <- cls_record_g1[["CL002"]][[whichtree]]$ctree
library(visNetwork)
library(RColorBrewer)
vninfo <- networkPlot(best_tree_cl002)
vn <- visNetwork(vninfo$nodes,vninfo$edges, width = "100%") %>% 
  visLegend(width=0.2,addNodes=vninfo$lnodes,useGroups=F)
vn %>%
  visSave("PHE_TB_transmission.html")


# Plot generation times and TimestoSampling distribution
# source getGenerationTimes and getTimesToSampling from transphylo_extras.R
# A wrapper for plotting histogram generation times for each cluster
# ... --- parameters to be passed to grid.arrange
plot_gen_times <- function(cls_record, ...){
  gs <- list()
  for(nm in names(cls_record)){
    gen_times <- map(cls_record[[nm]], ~ getGenerationTimes(.$ctree))
    gs[[nm]] <- ggplot(data.frame(gen_times = unlist(gen_times)), aes(gen_times)) +
      geom_histogram(bins = 20) + ggtitle(nm)  
    
    # The use of data.frame() in ggplot creates a new plot env for each ggplot object 
    # See here https://stackoverflow.com/questions/39799886/r-assigning-ggplot-objects-to-list-in-loop/39800861#39800861
    # for explanation of why the following doesn't work properly 
    # https://stackoverflow.com/questions/49978925/saving-ggplot-object)
    #gs[[nm]] <- qplot(unlist(gen_times), geom = "histogram", main = nm)
  }
  gridExtra::grid.arrange(grobs = gs, ...)
}

# A wrapper for plotting histogram times-to-sampling for each cluster
plot_times_to_samp <- function(cls_record, ...){
  gs <- list()
  for(nm in names(cls_record)){
    times_to_samp <- map(cls_record[[nm]], ~ getTimesToSampling(.$ctree))
    gs[[nm]] <- ggplot(data.frame(times_to_samp = unlist(times_to_samp)), aes(times_to_samp)) +
      geom_histogram(bins = 20) + ggtitle(nm)  
  }
  gridExtra::grid.arrange(grobs = gs, ...)
}
plot_gen_times(cls_record_g1, ncol = 3)
plot_gen_times(cls_record_g2, ncol = 3)  
plot_times_to_samp(cls_record_g1, ncol = 3)
plot_times_to_samp(cls_record_g2, ncol = 3)


# Get posterior summary table for generation and sampling times for all clusters
get_table_gtst <- function(cls_record){
  tibble(cluster_id = names(cls_record), 
         gen_time = map(cls_record, ~ sapply(lapply(., function(x) getGenerationTimes(x$ctree)), mean)),
         samp_time = map(cls_record, ~ sapply(lapply(., function(x) getTimesToSampling(x$ctree)), mean))) %>%
    mutate(mean_gt = map_dbl(gen_time, mean), sd_gt = map_dbl(gen_time, sd)) %>%
    mutate(mean_st = map_dbl(samp_time, mean), sd_st = map_dbl(samp_time, sd)) %>%
    select(-gen_time, -samp_time)  
}
get_table_gtst(cls_record_g1)
get_table_gtst(cls_record_g2)

# Plot number of unsampled cases for each cluster
plot_num_unsamp <- function(cls_record, ...){
  gs <- list()
  for(nm in names(cls_record)){
    num_unsamp <- get_num_unsampled(cls_record[[nm]])
    gs[[nm]] <- ggplot(data.frame(num_unsamp = num_unsamp), aes(num_unsamp)) +
      geom_bar() + ggtitle(nm)
  }
  gridExtra::grid.arrange(grobs = gs, ...)
}
plot_num_unsamp(cls_record_g1)
plot_num_unsamp(cls_record_g2)

# Get posterior summary table for number of unsampled cases
# source get_num_unsampled() from analysis4.Rmd
get_table_num_unsamp <- function(cls_record){
  tibble(cluster_id = names(cls_record),
         num_unsamp = map(cls_record, get_num_unsampled)) %>%
    mutate(mean_num_unsamp = map_dbl(num_unsamp, mean), sd_num_unsamp = map_dbl(num_unsamp, sd)) %>%
    select(-num_unsamp)
}
get_table_num_unsamp(cls_record_g1)
get_table_num_unsamp(cls_record_g2)


# Function to get max root-to-tip distance from a ctree, used in plotting against number of unsampled cases
get_max_rtt <- function(ctree){
  ttree <- extractTTree(ctree)
  ns <- length(ttree$nam) # number of sampled cases
  ttree <- ttree$ttree
  src <- which(ttree[, 3] == 0)
  max(ttree[1:ns,2] - ttree[src,1])
}

# Plot rtt distance with mean number of unsampled cases for two groups of clusters
plot_rtt_and_mean_unsamp <- function(cls_record1, cls_record2){
  tibble1 <- tibble(cluster_id = names(cls_record1), 
                    group_id = "G1",
                    max_rtt_dist = map(cls_record1, ~ sapply(., function(x) get_max_rtt(x$ctree))),
                    mean_max_rtt_dist = map_dbl(max_rtt_dist, mean),
                    mean_num_unsamp = get_table_num_unsamp(cls_record1)$mean_num_unsamp) %>%
    select(-max_rtt_dist)
  tibble2 <- tibble(cluster_id = names(cls_record2), 
                    group_id = "G2",
                    max_rtt_dist = map(cls_record2, ~ sapply(., function(x) get_max_rtt(x$ctree))),
                    mean_max_rtt_dist = map_dbl(max_rtt_dist, mean),
                    mean_num_unsamp = get_table_num_unsamp(cls_record2)$mean_num_unsamp) %>%
    select(-max_rtt_dist)
  ggplot(data = rbind(tibble1, tibble2), aes(mean_max_rtt_dist, mean_num_unsamp)) +
    geom_point(aes(color = group_id))
}
plot_rtt_and_mean_unsamp(cls_record_g1, cls_record_g2)

