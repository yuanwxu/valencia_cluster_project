# Functions to extract interesting results from the posterior
# These are cluster-level functions. 

source("~/Biomath/PHE_TB/analysis_functions.R")

# Get number of unsampled cases
get_num_unsampled <- function(record, ttrees = NULL){
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  purrr::map_int(ttrees, ~ sum(is.na(.$ttree[, 2])))
} 


# Get infector frequencies for each host in a cluster
# Return: a data frame containing 1) host ID 2) infector ("0" being the source) and 
#         3) freq, or probability that the host is infected by the "infector"
get_infector_freq <- function(record, ttrees = NULL){
  if(is.null(ttrees)){
    ttrees <- purrr::map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  host_id <- record[[1]]$ctree$nam
  df <- purrr::map_dfc(host_id, get_infectors, record = record, ttrees = ttrees)
  names(df) <- host_id
  df <- df %>% 
    gather(key = "host_id", value = "infector") %>% 
    group_by(host_id) %>%
    summarise(freq = list(table(infector)/length(record))) %>%
    mutate(infector = map(freq, names)) %>%
    unnest() %>%
    select(host_id, infector, freq) # just changing the order of columns
  df
}


# Get time of first transmission event for each host in a cluster, conditioned on transmission happened
# rates --- if not NULL, character vector of rates used to estimate time-trees, in which case "record"
#          will be "records" corresponding to these rates stacked together
# wt --- weight factor associated with each tree (see get_tree_weight in valencia_cluster_final.R)
# Return: a data frame
get_transmtime_df <- function(record, ttrees = NULL, rates = NULL, wt = NULL){
  if(is.null(ttrees)){
    ttrees <- purrr::map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  host_id <- record[[1]]$ctree$nam
  df <- tibble(host_id = host_id,
               inftime = purrr::map(host_id, get_inftime, record = record, 
                             ttrees = ttrees),  
               gentime = purrr::map(host_id, get_gentime, record = record, 
                             ttrees = ttrees, type = "min"))
  if(!is.null(rates)){
    df <- df %>% mutate(rate = list(rep(rates, each = length(record)/length(rates))))
  }
  df <- df %>%
    unnest() %>%
    mutate(transmtime = inftime + gentime)
  if(!is.null(wt)){
    df <- df %>% mutate(wt = rep(wt, length(!! host_id))) # "!!" used to  bypass the data frame and look directly in the context
  }
  df
}




