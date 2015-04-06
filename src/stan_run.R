library(rstan)
library(parallel)

## runs stan with caching of compiled Stan models
stan_run <- function(stanModel, ..., parallel=FALSE) {
  if(class(stanModel) == "stanfit") {
    stanExe <- stanModel
  } else {
    stanModel.rda <- gsub("stan$", "rda", stanModel)
    if(!file.exists(stanModel.rda) || file.info(stanModel.rda)$mtime < file.info(stanModel)$mtime) {
      cat("Model",stanModel,"needs recompilation.\n")
      args <- modifyList(list(...), list(file=stanModel, iter=0, warmup=0, chains=0))
      stanExe <- do.call(stan, args)
      saveRDS(stanExe, file=stanModel.rda)
    } else {
      cat("Loading cached stan model", stanModel, ".\n")
      stanExe = readRDS(stanModel.rda)
    }
  }
  if(parallel)
      return(stanParallel(fit=stanExe, ...))
  else
      return(stan(fit=stanExe, ...))
}

# mclapply backend, package parallel
stanParallel <- function(chains=4, ...) {
    args <- list(...)
    if("fit" %in% names(args)) {
        posterior <- mclapply(1:chains, FUN = function(i){ stan(chains=1, chain_id=i, refresh=-1, ...) })
    } else {
        stanModel <- stan(..., chains=0)
        posterior <- mclapply(1:chains, FUN = function(i){ stan(fit=stanModel, chains=1, chain_id=i, refresh=-1, ...) })
    }
    cat("STAN FIT FINISHED...\n")
    cat("MERGING CHAINS...\n")
    posterior <- sflist2stanfit(posterior) # see below for function definition
    invisible(posterior)
}

