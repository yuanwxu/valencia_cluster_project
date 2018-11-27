# Functions to extract interesting results from the posterior
# These are host-specific functions. 

# Get mean infection time for a sampled host 
# record --- transphylo output
# nm --- the name of the host
# ttrees --- list of transmission trees extracted from record, if not provided, will extract from record
#            ttrees <- map(record, ~ extractTTree(.$ctree))
get_mean_inftime <- function(record, ttrees = NULL, nm){
  stopifnot(nm %in% record[[1]]$ctree$nam)
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  out <- purrr::map_dbl(ttrees, ~ .$ttree[which(.$nam == nm),1])
  mean(out)
}

# Add a version that returns the raw posterior samples without taking the mean, useful for plotting histograms
get_inftime <- function(record, ttrees = NULL, nm){
  stopifnot(nm %in% record[[1]]$ctree$nam)
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  out <- purrr::map_dbl(ttrees, ~ .$ttree[which(.$nam == nm),1])
  out
}

# Get mean number of infectees for a sampled host
get_num_infectees <- function(record, ttrees = NULL, nm){
  stopifnot(nm %in% record[[1]]$ctree$nam)
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  out <- purrr::map_int(ttrees, ~ sum(.$ttree[,3] == which(.$nam == nm)))
  mean(out)
}

# Add a version that returns the raw posterior samples without taking the mean, useful for plotting histograms
get_infectees <- function(record, ttrees = NULL, nm){
  stopifnot(nm %in% record[[1]]$ctree$nam)
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  out <- purrr::map_int(ttrees, ~ sum(.$ttree[,3] == which(.$nam == nm)))
  out
}

# Get mean genearation time for a sampled host
# type --- either "mean" (default) which uses the mean of gen times of the given host in each tree,
#          or "min" which uses the minimum gen times, in which case the function returns 
#          the mean time between infection and first transmission event for the host in question.
get_mean_gentime <- function(record, ttrees = NULL, nm, type = "mean"){
  stopifnot(nm %in% record[[1]]$ctree$nam, type %in% c("mean", "min"))
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  if(type == "mean"){
    out <- purrr::map_dbl(ttrees, function(x){
      i <- which(x$nam == nm)
      j <- which(x$ttree[,3] == i)
      # NA incurred if the host didn't infect anyone before sampling
      return(ifelse(length(j) == 0, NA, mean(x$ttree[j,1] - x$ttree[i,1])))})  
  }
  if(type == "min"){
    out <- purrr::map_dbl(ttrees, function(x){
      i <- which(x$nam == nm)
      j <- which(x$ttree[,3] == i)
      # NA incurred if the host didn't infect anyone before sampling
      return(ifelse(length(j) == 0, NA, min(x$ttree[j,1] - x$ttree[i,1])))})
  }
  
  mean(out, na.rm = TRUE)
}

# Add a version that returns the raw posterior samples without taking the mean, useful for plotting histograms
get_gentime <- function(record, ttrees = NULL, nm, type = "mean"){
  stopifnot(nm %in% record[[1]]$ctree$nam, type %in% c("mean", "min"))
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  if(type == "mean"){
    out <- purrr::map_dbl(ttrees, function(x){
      i <- which(x$nam == nm)
      j <- which(x$ttree[,3] == i)
      # NA incurred if the host didn't infect anyone before sampling
      return(ifelse(length(j) == 0, NA, mean(x$ttree[j,1] - x$ttree[i,1])))})  
  }
  if(type == "min"){
    out <- purrr::map_dbl(ttrees, function(x){
      i <- which(x$nam == nm)
      j <- which(x$ttree[,3] == i)
      return(ifelse(length(j) == 0, NA, min(x$ttree[j,1] - x$ttree[i,1])))})
  }
  out # out may contain NA, when the host didn't infect anyone before sampling
}

# Probability of a sampled host having unsampled infector, conditioning on it's not index case
prob_infected_by_unsamp <- function(record, ttrees = NULL, nm){
  stopifnot(nm %in% record[[1]]$ctree$nam)
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  # Get indices of ttrees in which nm is not index case
  mask <- purrr::map_lgl(ttrees, function(x){
    i <- which(x$nam == nm)
    k <- x$ttree[i,3]
    return(k != 0)
    })
  out <- purrr::map_lgl(ttrees[mask], function(x){
    i <- which(x$nam == nm)
    k <- x$ttree[i,3]
    return(is.na(x$ttree[k,2]))
    })
  mean(out)
}

# Probability of a sampled host being the index case
prob_is_index_case <- function(record, ttrees = NULL, nm){
  stopifnot(nm %in% record[[1]]$ctree$nam)
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  out <- purrr::map_lgl(ttrees, function(x){
    i <- which(x$nam == nm)
    k <- x$ttree[i,3]
    return(k == 0)
  })
  mean(out)
}

# Get infectors of a sampled host, return "0" if the host is the source
get_infectors <- function(record, ttrees = NULL, nm){
  stopifnot(nm %in% record[[1]]$ctree$nam)
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  out <- purrr::map_chr(ttrees, function(x){
    i <- which(x$nam == nm)
    k <- x$ttree[i,3]
    return(ifelse(k == 0, "0", x$nam[k]))
  })
  out[is.na(out)] <- "Unsampled"
  out
}

# Get number of unsampled in n backward transmission events of a host
# i.e. for the host in question, look at if its infector, infector's infector, etc
# are sampled or not
get_number_unsamp_backward_n <- function(record, ttrees = NULL, nm, n = 1){
  stopifnot(nm %in% record[[1]]$ctree$nam)
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  
  get_it_one_tree <- function(ttree, nm, n){
    out <- integer(n)
    i <- which(ttree$nam == nm)
    k <- i
    for(ii in seq_along(out)){
      k <- ttree$ttree[k,3] # infector 
      if (k == 0) # found index case, break the loop
        break
      out[ii] <- is.na(ttree$ttree[k,2])
    }
    sum(out)
  }
  
  purrr::map_int(ttrees, get_it_one_tree, nm = nm, n = n)
}

# Get number of unsampled in n forward transmission events of a host
# i.e. for the host in question, look at if its infectee, infectee's infectee, etc
# are sampled or not. Because a host can be infected by only one host but can
# infect many others, this function is more complicated than the backward version,
# it sums all unsampled cases accumulated for each level of the transmission chain,
# for all secondary cases.
get_number_unsamp_forward_n <- function(record, ttrees = NULL, nm, n = 1){
  stopifnot(nm %in% record[[1]]$ctree$nam)
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  
  get_it_one_tree <- function(ttree, nm, n){
    out <- integer(n)
    k <- which(ttree$nam == nm)
    parent <- k
    for(ii in seq_along(out)){
      child <- integer()
      for(p in parent){
        if(!(p %in% ttree$ttree[,3])) # not an infector, continue to next node (if any)
          next
        child <- append(child, which(p == ttree$ttree[,3]))
      }
      out[ii] <- out[ii] + sum(is.na(ttree$ttree[child, 2]))  
      if(length(child) == 0)
        break
      parent <- child
    }
    sum(out)
  }
  
  purrr::map_int(ttrees, get_it_one_tree, nm = nm, n = n)
}


# Get mean number of unsampled infectees for a host, given that it infects other
get_cond_mean_unsamp_infectee <- function(record, ttrees = NULL, nm){
  stopifnot(nm %in% record[[1]]$ctree$nam)
  if(is.null(ttrees)){
    ttrees <- map(record, ~ TransPhylo::extractTTree(.$ctree))
  }
  
  mask <- purrr::map_lgl(ttrees, function(x){
    k <- which(x$nam == nm)
    any(x$ttree[,3] == k)
  })
  out <- purrr::map_int(ttrees[mask], function(x){
    k <- which(x$nam == nm)
    infectee_ind <- which(x$ttree[,3] == k)
    sum(is.na(x$ttree[infectee_ind,2]))
  })
  mean(out)
}


