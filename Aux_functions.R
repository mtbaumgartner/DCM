#### Auxiliar functions ####

#### Import all sheets from .xls into one list ####
read_excel_allsheets <- function(filename, tibble = FALSE) {
  library(readxl)
  # if you want tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  y <- lapply(x, data.matrix)
  for(i in 1:length(y)){
    y[[i]] <- y[[i]][,-1]
    row.names(y[[i]]) <- x[[i]][,1]
  }
  y
}

#### Coextinction simulations for each network under DCM ####
coextSims <- function(imatrix, rlow, rup, ext.thres, nsims){
  mods <- list()
  bar <- txtProgressBar(min = 0, max = nsims, style = 3)
  for(sim in 1:nsims){
    rvalue <- runif(1, rlow, rup)
    rplants <- rep(rvalue, ncol(imatrix))
    ranim <- rep(rvalue, nrow(imatrix))
    mods[[sim]] <- stochinteract(imatrix = imatrix, ext.thres = ext.thres,
                                 rplants = rplants, ranim = ranim)
    setTxtProgressBar(bar, sim)
  }
  return(mods)
}


#### Simulatios for all networks together under DCM ####
coextNetworks <- function(networks, rlow, rup, ext.thres, nsims){
  simulations <- list()
  for(net in 1:length(networks)){
    t.begin <- proc.time()
    
    imatrix <- networks[[net]]
    simulations[[net]] <- coextSims(imatrix = imatrix, rlow = rlow, rup = rup,
                                    ext.thres = ext.thres, nsims = nsims)
    
    time.used <- proc.time() - t.begin
    cat("\n","Simulation of web", net, "completed:", time.used[3]/60, "min\n")
  }
  return(simulations)
}


#### Functions for histograms of cascade degree and number of extinctions under DCM ####
coextDegHist <- function(network_sims, nsims){
  for(i in 1:length(network_sims)){
    sims <- c()
    for(j in 1:nsims){
      sims[j] <- max(network_sims[[i]][[j]][[1]]$degree)
    }
    h <- hist(sims, breaks = 50, plot = F)
    h$counts <- h$counts/sum(h$counts)
    plot(h, xlab = i)
  }
}

coextNumHist <- function(network_sims, nsims){
  for(i in 1:length(network_sims)){
    sims <- c()
    for(j in 1:nsims){
      sims[j] <- sum(network_sims[[i]][[j]][[1]]$n_extinctions)
    }
    h <- hist(sims, breaks = 50, plot = F)
    h$counts <- h$counts/sum(h$counts)
    plot(h, xlab = i)
  }
}
