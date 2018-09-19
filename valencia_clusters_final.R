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
names(cls_sts) <- cls_name

# Generate time-trees from a variety of clock rates, sampled from log-normal distribution with suitable 
# mean and variance.  Source tdaterAnalysis from analysis4.Rmd
# cls_mltree, cls_sts, sqlen --- mandatory parameters of of treedater
# n --- number of random variates of rates to sample from
# meanlog, sdlog --- same as in the description of rlnorm
# delta --- determine the margin of rate interval for any sampled rate, used in specifying meanRateLimits in treedater
# ... --- other arguments passed to treedater
# Return a list (of length n) of list (of length #cluster) of time-trees
get_timetrees_with_rates <- function(cls_mltree, cls_sts, sqlen, n, meanlog, sdlog, delta, ...){
  r <- rlnorm(n, meanlog, sdlog)
  out <- vector("list", length = n)
  for(i in seq_len(n)){
    rate_left <- ifelse(r[i]-delta < 0, r[i], r[i]-delta) # ensure the lower bound of rate is positive
    out[[i]] <- tdaterAnalysis(cls_mltree, cls_sts, sqlen, meanRateLimits = c(rate_left, r[i]+delta), ...)
  }
  names(out) <- format(round(r, 3), nsmall = 3) # use the rates to name the result, for easy reference
  out
}

set.seed(1)
lst_cls_timetree <- get_timetrees_with_rates(cls_mltree_unrooted, cls_sts, 1, n = 5, meanlog = -0.7, sdlog = 0.5, 
                                             delta = 0.2, strictClock = TRUE)

# Fix zero branch lengths of timed trees by shifting *all* branch lengths by 0.1
fix_edge_length <- function(timed_tree){
  if(any(dplyr::near(timed_tree$edge.length,0))){
    timed_tree$edge.length <- timed_tree$edge.length + 0.1
    return(timed_tree)
  }
  else # do nothing
    return(timed_tree)
}
for(i in seq_along(lst_cls_timetree)){
  lst_cls_timetree[[i]] <- map(lst_cls_timetree[[i]], fix_edge_length)
}

# Run TransPhylo
cls_lastdate <- map(lst_cls_timetree[[1]], ~ max(.$sts))
lst_cls_ptree <- vector("list", length(lst_cls_timetree))
for(i in seq_along(lst_cls_timetree)){
  lst_cls_ptree[[i]] <- map2(lst_cls_timetree[[i]], cls_lastdate, ptreeFromPhylo)
}
source("~/Biomath/TransPhylo/R/proposal.R")
source("~/Biomath/TransPhylo/R/computeHost.R")
source("infer_multiTTree_shareParam.R")

ws.shape <- 1.1; ws.scale <- 1/0.4
w.shape <- 1.3; w.scale <- 1/0.3
iters <- 4e4; thin <- 50
lst_cls_record <- vector("list", length(lst_cls_ptree))
for(i in seq_along(lst_cls_timetree)){ # run transphylo on all clusters for each clock rate
  lst_cls_record[[i]] <- infer_multiTTree_shareParam(lst_cls_ptree[[i]], w.shape, w.scale, ws.shape, ws.scale, 
                                                     mcmcIterations = iters, thinning = thin, 
                                                     share = c("neg","off.r","off.p","pi"))   
}


# Simple trace plot of shared parameters
plot_trace <- function(record){
  par(mfrow = c(1,3))
  plot(map_dbl(record, "neg"), type = "l", ylab = "neg", xlab = "MCMC step")
  plot(map_dbl(record, "off.r"), type = "l", ylab = "off.r", xlab = "MCMC step")
  plot(map_dbl(record, "pi"), type = "l", ylab = "pi", xlab = "MCMC step")  
}
plot_trace(lst_cls_record[[1]][[1]])


discard_burnin <- function(record_lst, burnin = 0.2){
  map(record_lst, function(x) x[round(length(x)*burnin):length(x)])
}

lst_cls_record2 <- vector("list", length(lst_cls_record)) # to store record with burnin disgarded
for(i in seq_along(lst_cls_record)){
  lst_cls_record2[[i]] <- discard_burnin(lst_cls_record[[i]], burnin = 0.5)  
}


# A helper function to set names of outer list as rates, and inner list as cluster names
set_lst_cls_names <- function(lst_cls_record, lst_name, cls_name){
  names(lst_cls_record) <- lst_name
  map(lst_cls_record, setNames, nm = cls_name)
}
# Set names for easy reference
# lst_cls_record <- set_lst_cls_names(lst_cls_record, lst_name = names(lst_cls_timetree), cls_name)
lst_cls_record2 <- set_lst_cls_names(lst_cls_record2, lst_name = names(lst_cls_record), cls_name)
# saveRDS(lst_cls_record, "Runs/lst_cls_record.rds") # save to external place
rm(lst_cls_record)

# A helper function to pool the results of lst_cls_record by collapsing "clock rates" so the output
# is a list of length #clusters, where each element of the list stores the transphylo run of length
# equal to the number of simulated rates times the number of transphylo iterations, i.e. for each 
# cluster merge the results for all rates together
# === Experiment with a simpler R object that mimic data structure of lst_cls_record ===
# foo <- list(rate1 = list(cl1 = list(list(a=1,t="t1"),list(a=2,t="t2")), 
#                          cl2 = list(list(a=1,t="t1"),list(a=2,t="t2"))), 
#             rate2 = list(cl1 = list(list(a=2,t="t2"),list(a=3,t="t2")),
#                          cl2 = list(list(a=2,t="t2"),list(a=3,t="t2"))),
#             rate3 = list(cl1 = list(list(a=3,t="t3"),list(a=4,t="t3")), 
#                          cl2 = list(list(a=3,t="t3"),list(a=4,t="t3"))))
# transform_lst_cls_record(foo)
# ======================================================================================
transform_lst_cls_record <- function(lst_cls_record){
  lst_cls_record_T <- purrr::transpose(lst_cls_record)
  map(lst_cls_record_T, purrr::flatten)
}

cls_record2 <- transform_lst_cls_record(lst_cls_record2)


# Posterior analysis
#
# ==== Trace plot ====
# Plot traces of parameters as function of MCMC step, colored by different simulated clock rate
# may also be used as convergence diagnostics
# lst_cls_record --- transphylo run on all clusters, with different clock rates
plot_trace_by_rate <- function(lst_cls_record){
  cls_record <- transform_lst_cls_record(lst_cls_record)
  iters <- length(lst_cls_record[[1]][[1]])
  nrates <- length(lst_cls_record)
  df_trace <- tibble(mcmc_step = rep(1:iters, nrates),
                     neg = map_dbl(cls_record[[1]], "neg"),
                     off.r = map_dbl(cls_record[[1]], "off.r"),
                     pi = map_dbl(cls_record[[1]], "pi"),
                     rate = rep(names(lst_cls_record), each = iters)) %>%
    gather(key = "param", value = "value", neg:pi)
  df_trace %>%
    # filter(mcmc_step > 200) %>% # looks like neg is still in burnin
    ggplot(aes(mcmc_step, value)) +
    geom_line(aes(color = rate)) +
    facet_wrap(~ param, scales = "free_y") 
}
# Add info or decorations
plot_trace_by_rate(lst_cls_record2) + labs(x = "MCMC step * 50 after burnin")



# ==== Transmission before symptom onset? ====
source("analysis_functions2.R")
cls_ttrees <- setNames(vector("list", length(cls_record2)), names(cls_record2))
for(i in seq_along(cls_record2)){ # extract ttree from record to save computations required in analysis functions
  cls_ttrees[[i]] <- map(cls_record2[[i]], ~ extractTTree(.$ctree))
}

cls_transmtime_df <- map2(cls_record2, cls_ttrees, get_transmtime_df, rates = names(lst_cls_record2))

# Read epi data for each cluster
# path --- path to directory
# pattern --- file name regex, same as in list.files
# cls_name --- char vector, names of clusters for which epi is desired
# ... --- other params passed on to list.files
read_epi <- function(path = ".", pattern, cls_name){
  f <- list.files(path, pattern)
  f <- f[f %in% paste0("epi_info_", cls_name, ".xlsx")]

  out <- purrr::map(paste0(path, f), readxl::read_excel, na = c("", "NO DATE"))
  names(out) <- map_chr(out, ~ .$`Cluster number`[1])
  out
}

cls_epi <- read_epi(path = "~/Biomath/MLRA/clusters_alignments_valencia/Data/new/", 
                    pattern = "epi_info_CL", cls_name)

#cl045_transmtime_df <- get_transmtime_df(cls_record2[["CL045"]], cls_ttrees[["CL045"]], names(lst_cls_record2))
#cl045_epi <- readxl::read_excel("~/Biomath/MLRA/clusters_alignments_valencia/Data/new/epi_info_CL045.xlsx")

# Helper function to get reference time points: symptime onset and diagnosis/sampling time for a cluster 
# epi --- epi data of a cluster
# sts --- sampling/diagnosis times for the cluster returned from get_sts2Treedater
# hosts --- character vector of host id's in the cluster, the "nam" component of "ctree"
get_symp_dx_time_df <- function(epi, sts, hosts){
  stime <- lubridate::decimal_date(with(epi, setNames(`Symptomps onset`, paste0("G", `Isolate ID`))))
  tibble(host_id = hosts, symp_time = stime[host_id], dgns_time = sts$sts[host_id])
}

cls_hosts <- map(cls_record2, ~ .x[[1]]$ctree$nam)
cls_times_df <- pmap(list(cls_epi, cls_sts, cls_hosts), get_symp_dx_time_df)
cls_times_df <- map(cls_times_df, gather, key = "key", value = "value", symp_time:dgns_time)


#cl045_times_df <- get_symp_dx_time_df(cl045_epi, cls_sts[["CL045"]], hosts = cls_record2[["CL045"]][[1]]$ctree$nam)
#cl045_times_df <- cl045_times_df %>%
#  gather("key", "value", symp_time:dgns_time)

# Box plot of time of first transmission for one clock rate
# r --- char, which clock rate to show
# show_cls --- char vector, clusters to show
plot_transmtime <- function(cls_transmtime_df, cls_reftimes_df, r, show_cls = names(cls_transmtime_df)){
  stopifnot(all.equal(names(cls_reftimes_df), names(cls_transmtime_df)))
  stopifnot(r %in% unique(cls_transmtime_df[[1]]$rate))

  cls_name <- names(cls_transmtime_df)
  df <- map2(cls_transmtime_df, cls_name, ~ mutate(.x, cluster_id = .y)) %>%
    map_dfr(~.) # add column cluster_id then rbind all data frames
  
  df_ref <- map2(cls_reftimes_df, cls_name, ~ mutate(.x, cluster_id = .y)) %>%
    map_dfr(~.) # add column cluster_id then rbind all data frames
  
  df %>% 
    filter(rate == r, cluster_id %in% show_cls) %>%
    na.omit() %>%
    ggplot(aes(x = host_id, y = transmtime)) +
    geom_boxplot(varwidth = TRUE) +
    facet_wrap(~ cluster_id, scales = "free") +
    geom_point(aes(y = value, shape = key, color = key), 
               data = df_ref %>% filter(cluster_id %in% show_cls), size = 3) +
    scale_shape_manual(values = c(17,15)) +
    labs(x = "Host", y = "Time", shape = "", color = "") + coord_flip()
}

plot_transmtime(cls_transmtime_df, cls_times_df, "0.363", show_cls = cls_name[1:6]) 


#cl045_transmtime_df %>%
#  filter(rate == "0.363") %>%
#  na.omit() %>%
#  ggplot(aes(x = host_id, y = transmtime)) +
#  geom_boxplot(varwidth = TRUE) +
#  geom_point(aes(y = value, shape = key, color = key), data = cl045_times_df, size = 4) +
#  scale_shape_manual(values = c(17,15)) +
#  labs(x = "Host", y = "Time", title = "Box plot of time of first transmission event in CL045", shape = "", color = "")

## THE OLD PLOT SHOWING ALL RATES
#ggplot(cl045_transmtime_df, aes(x = transmtime)) +
#  geom_freqpoly(aes(color = rate), binwidth = 0.1) +
#  facet_wrap(~host_id, scales = "free_y") + # because count can differ a lot across hosts and rates
#  geom_vline(aes(xintercept = symptime), data = cl045_symptime_df, linetype = "dashed", size = 1.5, color = "blue") +
#  labs(x = "Time of first transmission event", title = "CL045")
       

# Function to plot historgrams of first transmission time or return invisibly the posterior fraction of 
# non-transmission events for each host -- this is used to provide reference to the plot because not all
# histograms have the same total counts.
# symptime --- vector of symptom onset times, named by host id
# ... --- parameter to be passed on to grid.arrange
transm_before_symptom <- function(record, symptime, plot = TRUE, ...){
  ttree <- map(record, ~ extractTTree(.$ctree))
  hosts <- ttree[[1]]$nam
  gs <- list()
  fni <- numeric(length(hosts)) # store for each host x fraction of posterior trees where x did not infect anyone
  names(fni) <- hosts
  for(nm in hosts){
    inftime_nm <- get_inftime(record, ttree, nm)  
    gentime_nm <- get_gentime(record, ttree, nm, type = "min")
    transmtime_nm <- inftime_nm + gentime_nm # may contain some NAs if there's no transmission event
    fni[nm] <- sum(is.na(transmtime_nm))
    gs[[nm]] <- ggplot(data.frame(first_transm_time = transmtime_nm), aes(first_transm_time)) +
      geom_histogram() + ggtitle(nm) +
      scale_x_continuous(breaks = scales::pretty_breaks(5)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    if(nm %in% names(symptime)) # add ref line of symptom onset if it's known
      gs[[nm]] <- gs[[nm]] + geom_vline(xintercept = symptime[[nm]], color = "red")
  }
  plt <- gridExtra::grid.arrange(grobs = gs, ...)
  if(plot){
   plt
  }
  invisible(list(plt=plt, fni=fni/length(record)))
}

#tbs_cl2 <- transm_before_symptom(cls_record2[["CL002"]], symptime_cl2, ncol = 4)
#tibble(host = names(tbs_cl2$fni), fni = tbs_cl2$fni) %>%
#  ggplot(aes(reorder(host, -fni), fni)) +
#  geom_col() + ylab("Probability of sampling before infecting others") + xlab("Host")



# ==== Unsampled cases ====
# NOTE: "slice(1:n())" a trick to make the tibble work with unnest
cls_unsamp_df <- tibble(cluster_id = names(cls_record2),
                        unsamp = map2(cls_record2, cls_ttrees, get_num_unsampled),
                        rate = list(rep(names(lst_cls_record2), each = length(lst_cls_record2[[1]][[1]])))) %>%
  slice(1:n()) %>%                  
  unnest()
cls_unsamp_df %>%
  filter(cluster_id %in% cls_name[1:10]) %>%
  ggplot(aes(reorder(cluster_id, unsamp, FUN = median), unsamp)) + geom_boxplot(aes(fill = rate)) +
  xlab("Cluster") + ylab("Number of unsampled cases")
cls_unsamp_df %>%
  filter(cluster_id %in% cls_name[11:21]) %>%
  ggplot(aes(reorder(cluster_id, unsamp, FUN = median), unsamp)) + geom_boxplot(aes(fill = rate)) +
  xlab("Cluster") + ylab("Number of unsampled cases")
       


# ==== Was the index case the first diagnosed case? ====
# sts --- named vector of sampling times for the sampled cases in record
# title --- title of plot, e.g. name of cluster 
# max_to_show --- number of cases to show in plot, default showing all cases.
#                 To make plot more compact, specify this to be some integer n
#                 so that no more than n cases will be shown in the plot.
plot_index_first_dgns <- function(record, sts, title, max_to_show = NULL){
  dtime <- c(sts, Unsampled = NA)
  data <- tibble(index_case = map_chr(record, "source")) %>%
    count(index_case, sort = TRUE) %>%
    mutate(dgns_time = map_dbl(index_case, ~ dtime[[.]])) %>%
    mutate(dgns_time = format(round(dgns_time, 2), nsmall = 2)) # make them characters with 2 signf digits
  
  ns <- ifelse(is.null(max_to_show), nrow(data), min(max_to_show, nrow(data)))
  data %>%
    slice(1:ns) %>%
    ggplot(aes(index_case, n)) + 
    geom_col(aes(fill = dgns_time)) +
    ggtitle(title)
}

# Get index case data for all clusters
get_cls_index_df <- function(cls_record, cls_ttrees = NULL, cls_sts){
  if(is.null(cls_ttrees)){
    cls_ttrees <- setNames(vector("list", length(cls_record)), names(cls_record))
    for(i in seq_along(cls_record)){ 
      cls_ttrees[[i]] <- map(cls_record[[i]], ~ extractTTree(.$ctree))
    }
  }
  
  cls_dtime <- map(cls_sts, ~ c(.$sts, Unsampled = NA))
  lst_df <- vector("list", length(cls_record))
  for(i in seq_along(lst_df)){
    lst_df[[i]] <- tibble(cluster_id = names(cls_record)[i], 
                          index_case = map_chr(cls_record[[i]], "source")) %>%
      count(cluster_id, index_case, sort = TRUE) %>%
      mutate(dgns_time = map_dbl(index_case, ~ cls_dtime[[i]][[.]]))
  }
  map_dfr(lst_df, ~ .)
}

cls_index_df <- get_cls_index_df(cls_record2, cls_ttrees, cls_sts)

# Showing index case distribution of all clusters in stacked bar chart
# max_case --- max #cases to show for each cluster
# min_n --- only cases whose counts are greater than *min_n* will be shown, this is more useful than
#           max_case because it ensures that the host label will be visible within each stacked bar
# cls --- char vector, clusters to show in plot, default all clusters
plot_index_all <- function(df, max_case = NULL, min_n = 1, cls = cls_name){
  if(!is.null(max_case)){
    lst_df <- split(df, df$cluster_id)
    len <- map_int(lst_df, nrow)
    ns <- ifelse(max_case < len, max_case, len)
    df <- map2(lst_df, ns, ~ slice(.x, 1:.y)) %>%
      map_dfr(~.)
  }
  
  df %>%
    mutate(index_case = ifelse(index_case == "Unsampled", "U", index_case)) %>% # for concise symbol in plotting
    filter(cluster_id %in% cls, n > min_n) %>%
    ggplot(aes(cluster_id, n)) + 
    geom_col(aes(fill = dgns_time), color = "black") +
    geom_text(aes(label = index_case), position = position_stack(vjust = 0.5))
}

plot_index_all(cls_index_df, min_n = 60, cls = cls_name[1:10]) # show all cases with greater than 3% posterior counts


# Alternative pie chart
plot_index_first_dgns_pie <- function(record, sts, title){
  dtime <- c(sts, Unsampled = NA)
  data <- tibble(index_case = map_chr(record, "source")) %>%
    count(index_case) %>%
    mutate(freq = n / sum(n))
  data <- data %>% 
    mutate(dgns_time = map_dbl(index_case, ~ dtime[[.]])) %>%
    mutate(dgns_time = format(round(dgns_time, 2), nsmall = 2))
  
  ggplot(data, aes(x = "", y = freq, fill = dgns_time)) + 
    geom_col(width = 1) +
    geom_text(aes(label = index_case), position = position_stack(vjust = 0.5)) +
    coord_polar("y") +
    labs(title = title)  
}

cls_plt_first_dgns <- vector("list", length(cls_record2)) # to store plots
for(i in seq_along(cls_record2)){
  cls_plt_first_dgns[[i]] <- plot_index_first_dgns(cls_record2[[i]], cls_sts[[i]]$sts, title = cls_name[[i]], 
                                                   max_to_show = 5)
  #cls_plt_first_dgns[[i]] <- plot_index_first_dgns_pie(cls_record2[[i]], cls_sts[[i]]$sts, title = cls_name[[i]])
  cls_plt_first_dgns[[i]] <- cls_plt_first_dgns[[i]] + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) # adjust angle of x axis for easier inspection
}
gridExtra::grid.arrange(grobs = cls_plt_first_dgns[1:6], ncol = 3)



#  ==== Visualize trees with diagrammeR ====
library(DiagrammeR)
# ncolors --- how many colors should be used to distinguish different infection times
# Return the graph for rendering and further modifications
visz_transm_tree <- function(ttree, ncolors = 5, alt_label = FALSE){
  ncases <- nrow(ttree$ttree)
  nunsamp <- ncases - length(ttree$nam)
  node_df <- create_node_df(n = ncases, label = c(ttree$nam, rep("Unsampled", nunsamp)),
                            shape = rep(c("circle","rectangle"), c(ncases - nunsamp, nunsamp)),
                            time_inf = ttree$ttree[,1],
                            time_samp = ttree$ttree[,2])
  if(alt_label){
    node_df$label[1:length(ttree$nam)] <- as.character(1:length(ttree$nam))
  }
  
  edge_df <- create_edge_df(from = ttree$ttree[,3], to = 1:ncases)
  
  # Delete from edge_df the row where infector is 0
  edge_df <- edge_df[-which(edge_df$from == 0), ]
  graph <- create_graph(node_df, edge_df)

  # Not working yet, wait to see if the package maintainer fix the bug
  # pal <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(ncolors) # generate 100 blues
  # Source: https://www.r-bloggers.com/r-using-rcolorbrewer-to-colour-your-figures-in-r/

  graph <- graph %>% 
    colorize_node_attrs(node_attr_from = time_inf, node_attr_to = fillcolor, 
                        cut_points = seq(min(node_df$time_inf), max(node_df$time_inf)+0.1, length.out = ncolors),
                        palette = "YlOrRd")
  graph
}

# consensus tree from the combined posterior
visz_transm_tree(consTTree(cls_record2[["CL045"]], burnin = 0)) %>%
  render_graph(title = "CL045")

# My select_weighted_MAP_tree function that returns the weighted maximum a-posteriori transmission tree
# where the weight for each tree is given by the posterior probability of #unsampled present in that tree
# For PHE datasest, MAP tree = weighted MAP tree, so no change even each tree is weighted by the prob #unsampled 
# in that tree. Consider introducing a hyperparam lambda in front of log: lambda * log(ct_unsamp / length(unsamp)),
# where lambda > 1, this results in exponential decay, as the prob. decrease, and this does change the tree
select_weighted_MAP_tree <- function(record, ttrees = NULL, burnin = 0.2){
  log_pos_prob <- map_dbl(record, "pTTree") + map_dbl(record, "pPTree")
  
  unsamp <- get_num_unsampled(record, ttrees)
  df <- data.frame(uns = unsamp) %>% dplyr::count(uns)
  ind <- match(unsamp, df$uns) # know which indices to look for in the count table 
  ct_unsamp <- df[ind, ] %>% dplyr::pull(n) # select as vector
  log_pos_prob <- log_pos_prob + log(ct_unsamp / length(unsamp)) # weight by freq of #unsampled
  which.max(log_pos_prob)
}



# ==== Who infected case X? ====
cls_infector_freq <- map2(cls_record2, cls_ttrees, get_infector_freq)

# Plot diagramm showing strength of infection between cases in a cluster
# df --- the data frame of infector frequences for the cluster, returned from get_infector_freq()
# magnify_edge --- the factor used to magnify the width of edges, increase this if low probability edges are not visible
# alt_label --- alternative label, TRUE if want to use natural sequence (1,2,...) as labels. This is useful if either
#               the labels are too long, or for data privacy and confidential reasons.
# ... --- parameters to be passed on to render_graph()
# If alt_label is TRUE, will return a list containing the graph object and the original labels
plot_infector_freq <- function(df, magnify_edge = 5, alt_label = FALSE, ...){
  nodes <- unique(c(df$host_id, df$infector))
  node_df <- create_node_df(n = length(nodes),
                            label = nodes,
                            shape = "circle",
                            fillcolor = dplyr::case_when(nodes == "0" ~ "grey", 
                                                         nodes == "Unsampled" ~ "lightblue",
                                                         TRUE ~ "lavender"))
  if(alt_label){
    node_df$label <- as.character(seq_along(nodes))
    node_df$label[which(nodes == "Unsampled")] <- "Unsampled"
    node_df$label[which(nodes == "0")] <- "0"
  }
                                  
  from <- match(df$infector, nodes)
  to <- match(df$host_id, nodes)
  edge_df <- create_edge_df(from = from, 
                            to = to,
                            color = "gray40",
                            penwidth = df$freq * magnify_edge)
  g <- create_graph(node_df, edge_df) 
  if(alt_label)
    return(list(graph = g, labels = nodes, from = df$infector, to = df$host_id))
  else
    render_graph(g, ...) # show graph
}
plot_infector_freq(cls_infector_freq[["CL045"]])
plot_infector_freq(cls_infector_freq[["CL016"]])






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


# Get posterior summary table for generation and sampling times for all clusters
get_table_gtst <- function(cls_record){
  tibble(cluster_id = names(cls_record), 
         gen_time = map(cls_record, ~ sapply(lapply(., function(x) getGenerationTimes(x$ctree)), mean)),
         samp_time = map(cls_record, ~ sapply(lapply(., function(x) getTimesToSampling(x$ctree)), mean))) %>%
    mutate(mean_gt = map_dbl(gen_time, mean), sd_gt = map_dbl(gen_time, sd)) %>%
    mutate(mean_st = map_dbl(samp_time, mean), sd_st = map_dbl(samp_time, sd)) %>%
    select(-gen_time, -samp_time)  
}


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


# Get posterior summary table for number of unsampled cases
get_table_num_unsamp <- function(cls_record){
  tibble(cluster_id = names(cls_record),
         num_unsamp = map(cls_record, get_num_unsampled)) %>%
    mutate(mean_num_unsamp = map_dbl(num_unsamp, mean), sd_num_unsamp = map_dbl(num_unsamp, sd)) %>%
    select(-num_unsamp)
}




