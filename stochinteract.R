####################################################
### Dependent Random-search Coextinction Model (DCM)
### Coded by Matheus T. Baumgartner
### (March 2019)
####################################################

#' Function to run extinction cascade
#' @param imatrix quantitative interaction matrix (plants as columns, animals as rows)
#' @param ext.thres threshold for extinction of potentially isolated species
#' @param rplants intrinsic demographic dependence of plants on animals
#' @param ranim intrinsic demographic dependence of animals on plants
#' @param distPlants dissimilarities among plant species
#' @param distAnim dissimilarities among animal species


stochinteract <- function(imatrix, ext.thres, rplants, ranim,
                          distPlants = NULL, distAnim = NULL){
  
  
  #---- Initial checks ----
  if(class(imatrix)!="matrix" || (nrow(imatrix)+ncol(imatrix))<3){stop("'imatrix' must be of class 'matrix' (animals on rows, plants on columns) with at least three species overall")}
  if(is.null(distPlants) || is.null(distAnim)){warning("Distances were not fully provided and were calculated from 'imatrix' using Bray-Curtis dissimilarities")}
  if(!is.null(distAnim)){warning("Dissimilarities among animals will be standardized")}
  if(!is.null(distPlants)){warning("Dissimilarities among plants will be standardized")}
  
  
  #---- Defining start variables ----
  npla <- ncol(imatrix)
  nanim <- nrow(imatrix)
  degplants <- colSums(imatrix > 0)
  deganim <- rowSums(imatrix > 0)
  plants <- 1:npla
  animals <- 1:nanim
  plantNA <- 1:npla
  animNA <- 1:nanim
  
  degree_when_lost_plants <- c()
  degree_when_lost_animals <- c()
  
  
  #---- Calculating dissimilarity matrices ----
  if(is.null(distPlants)){ # if dissimilarities among plants are not provided
    
    # calculate using Bray-Curtis dissimilarity
    distPlants <- matrix(nrow = nrow(t(imatrix)), ncol = nrow(t(imatrix)))
    
    for(j in 1:nrow(t(imatrix))){
      for(k in 1:nrow(t(imatrix))){
        distPlants[j,k] <- sum(abs(t(imatrix)[j,]-t(imatrix)[k,]))/sum(t(imatrix)[j,]+t(imatrix)[k,])
      }
    }
    
  }else{ # standardize dissimilarities if provided
    distPlants <- as.matrix(distPlants)
    distPlants <- (distPlants-min(distPlants))/(max(distPlants)-min(distPlants))
  }
  
  if(is.null(distAnim)){
    distAnim <- matrix(nrow = nrow(imatrix), ncol = nrow(imatrix))
    
    for(j in 1:nrow(imatrix)){
      for(k in 1:nrow(imatrix)){
        distAnim[j,k] <- sum(abs(imatrix[j,]-imatrix[k,]))/sum(imatrix[j,]+imatrix[k,])
      }
    }
    
  }else{
    distAnim <- as.matrix(distAnim)
    distAnim <- (distAnim-min(distAnim))/(max(distAnim)-min(distAnim))
  }
  
  
  #---- Calculating dependence matrices ----
  M <- array(0, dim=c(nanim, npla, 2))
  for(i in 1:npla){
    M[,i,1] <- imatrix[,i]/sum(imatrix[,i])
  } # matrix of plant dependence on each animal
  
  for(i in 1:nanim){
    M[i,,2] <- imatrix[i,]/sum(imatrix[i,])
  } # matrix of animal dependence on each plant
  
  
  #---- Choosing species for primary extinction ----
  coext_animals <- c()
  coext_plants <- c()
  
  plant_anim <- sample(c("animal","plant"), 1) # choose plant or animal for primary extinction
  if(plant_anim == "animal"){
    alive <- animals[!is.na(animNA)] # from those animals alive
    coext_animals <- sample(alive, 1) # select one animal to extinguish
    degree_when_lost_animals <- 1
  }else{
    alive <- plants[!is.na(plantNA)] # from those plants alive
    coext_plants <- sample(alive, 1) # select one plant to extinguish
    degree_when_lost_plants <- 1
  }
  
  # save original values and extinguish selected plant or animal species
  orimatrix <- imatrix
  
  imatrix[coext_animals,] <- 0
  imatrix[,coext_plants] <- 0
  
  
  # final list of animals/plants which were alive in the original community and were extinct primarily + cascade
  lostanimals <- coext_animals
  lostplants <- coext_plants
  ext_ID <- c(lostanimals, lostplants) # which species began the cascade
  
  
  #---- Defining some outputs (some updated in cascade) ----
  degree_table <- data.frame()
  N_likely <- 0 # number of species likely to extinction
  N_newint <- 0 # number of new interactions
  N_enhance <- 0 # number of interactions enhanced (compensating)
  N_weak <- 0 # number of calculated new with weak potential
  
  
  #---- Coextinction cascade ----
  equilibrium <- FALSE
  degree <- 1
  degree_table <- data.frame(degree, guild=factor(plant_anim, levels=c("animal","plant")), n_extinctions = 1, ext_ID)
  
  while(equilibrium==FALSE){ # CASCADE LOOP
    
    ext_animals <- coext_animals
    ext_plants <- coext_plants
    plantNA[ext_plants] <- NA # those already extinct recieve NAs
    animNA[ext_animals] <- NA
    
    aleft <- animals[!is.na(animNA)] # those left are the ones that do not have NAs
    pleft <- plants[!is.na(plantNA)]
    
    asafe <- c()
    psafe <- c()
    
    # DEPENDENT RANDOM-SEARCH MODEL
    if(length(ext_plants)>0){ #---- if some plant has gone extinct ----
      
      for(i in 1:length(ext_plants)){ # for each extinct plant
        
        Pij <- M[aleft,ext_plants[i],2] * ranim[aleft] * mean(distPlants[ext_plants[i],]) # calculate probability of extinction of i following the extinction of j
        unlucky <- Pij > runif(length(aleft))  # check for animals that are likely to become coextinct
        
        if(length(aleft[unlucky])>0){ # if there is some animal likely to coextinction
          N_likely <- N_likely + length(aleft[unlucky])  # update
          
          for(j in 1:length(aleft[unlucky])){ # for each animal likely to coextinction
            
            if(sum(imatrix[aleft[unlucky][j],]!=0)==0){ # if unlucky species interacted only with j
              
              # CREATE NEW INTERACTIONS
              if(Pij[unlucky][j] < ext.thres){ # if Pij is lower than the extinction threshold
                # create new interactions
                eligible <- colSums(imatrix[,pleft] > 0) <= degplants[pleft] + 1 # check for eligible new partners
                
                if(any(eligible == TRUE)){ # any elegible species for new interactions
                  
                  if(length(eligible) > deganim[aleft[unlucky][j]]){ # if the number of elegible is higher than the wish of j
                    asel <- deganim[aleft[unlucky][j]] # select all wish
                  }else{
                    asel <- sum(eligible) # select only possible
                  }
                  
                  chosen_spp <- as.numeric(sample(as.character(pleft[eligible]), asel, prob = (1 + 10e-7) - distPlants[ext_plants[i],pleft[eligible]])) # dependent random-search for new interaction
                  
                  # establish new interaction
                  imatrix[aleft[unlucky][j], chosen_spp] <- orimatrix[aleft[unlucky][j],ext_plants[i]] ^ (1 - distPlants[ext_plants[i], chosen_spp])
                  
                  N_newint <- N_newint + length(chosen_spp)
                  
                }else{
                  coext_animals <- c(coext_animals, aleft[unlucky][j])
                }
                
              }
              
              if(Pij[unlucky][j] >= ext.thres){ # if Pij is higher than the extinction threshold
                # do not create new interactions and species i is extinct
                coext_animals <- c(coext_animals, aleft[unlucky][j])
                
                N_weak <- N_weak + 1 # update
              }
              
              
            }else{ # if unlucky species interacted with more than j
              
              yesno <- sample(c(0,1), 1) # extinct or not?
              if(yesno==1){ # if extinct == YES
                coext_animals <- c(coext_animals, aleft[unlucky][j])
              }else{
                
                # COMPENSATE = ENHANCE INTERACTIONS
                # enhancement is proportional to its dependence on extinct j
                # (multiply row by (1 + dependence of i on j))
                imatrix[aleft[unlucky][j],] <- imatrix[aleft[unlucky][j],] * (1 + M[aleft[unlucky][j],ext_plants[i],2]) # compensate
                asafe <- c(asafe, aleft[unlucky][j]) # species i is no longer likely to extinction
                
                N_enhance <- N_enhance + 1 # update
              }
            }
          }
        }else{ # if there is no animal likely to coextinction
          coext_animals <- c(coext_animals, aleft[unlucky]) # no animal is coextinct
        }
      }
      
      unlucky[asafe] <- FALSE # species are safe
      
      # Update outputs
      coext_plants <- c()
      coext_animals <- unique(coext_animals)
      animNA[coext_animals] <- NA
      lostanimals <- c(lostanimals, coext_animals)
      imatrix[coext_animals,] <- 0
      orimatrix <- imatrix
      
      # Recalculate dependencies
      for(i in 1:nanim){
        if(sum(imatrix[i,])==0){
          M[i,,2] <- 0
        }else{
          M[i,,2] <- imatrix[i,]/sum(imatrix[i,])
        }
      }
      
      if(length(coext_animals)>0){ # if any animal has gone coextinct
        degree <- degree + 1 # increase degree
        degree_when_lost_animals <- c(degree_when_lost_animals, rep(degree, length(coext_animals)))
        degree_table[degree,] <- data.frame(degree, "animal", length(coext_animals))
      }
      
    }else{ #---- if some animal has gone extinct ----
      
      for(i in 1:length(ext_animals)){ # for each extinct animal
        
        Pij <- M[ext_animals[i],pleft,1] * rplants[pleft] * mean(distAnim[ext_animals[i],]) # calculate probability of extinction of i following the extinction of j
        unlucky <- Pij > runif(length(pleft)) # check for animals that are likely to become coextinct
        
        if(length(pleft[unlucky])>0){ # if there is some plant likely to coextinction
          N_likely <- N_likely + length(pleft[unlucky]) # update
          
          for(j in 1:length(pleft[unlucky])){ # for each plant likely to coextinction
            
            if(sum(imatrix[,pleft[unlucky][j]])==0){ # if unlucky species interacted only with j
              
              # CREATE NEW INTERACTIONS
              if(Pij[unlucky][j] < ext.thres){ # if Pij is lower than the extinction threshold
                # create new interactions
                eligible <- rowSums(imatrix[aleft,] > 0) <= deganim[aleft] + 1
                
                if(any(eligible == TRUE)){ # any elegible species for new interactions
                  
                  if(length(eligible) > degplants[pleft[unlucky][j]]){ # if the number of elegible is higher than the wish of j
                    psel <- degplants[pleft[unlucky][j]] # select all wish
                  }else{
                    psel <- sum(eligible) # select only possible
                  }
                  
                  chosen_spp <- as.numeric(sample(as.character(aleft[eligible]), psel, prob = (1 + 10e-7) - distAnim[ext_animals[i],aleft[eligible]])) # dependent random-search for new interaction
                  
                  # establish new interaction
                  imatrix[chosen_spp, pleft[unlucky][j]] <- orimatrix[ext_animals[i], pleft[unlucky][j]] ^ (1 - distAnim[ext_animals[i], chosen_spp])
                  
                  N_newint <- N_newint + length(chosen_spp)
                  
                }else{
                  coext_plants <- c(coext_plants, pleft[unlucky][j])
                }
                
              }
              
              if(Pij[unlucky][j] >= ext.thres){ # if Pij is higher than the extinction threshold
                # do not create new interactions and species i is extinct
                coext_plants <- c(coext_plants, pleft[unlucky][j])
                
                N_weak <- N_weak + 1 # update
              }
              
              
            }else{ # if unlucky species interacted with more than j
              
              yesno <- sample(c(0,1), 1) # extinct or not?
              if(yesno==1){ # if extinct == YES
                coext_plants <- c(coext_plants, pleft[unlucky][j])
              }else{
                
                # COMPENSATE = ENHANCE INTERACTIONS
                # enhancement is proportional to its dependence on extinct j
                # (multiply row by (1 + dependence of i on j))
                imatrix[,pleft[unlucky][j]] <- imatrix[,pleft[unlucky][j]] * (1 + M[ext_animals[i], pleft[unlucky][j],1]) # compensate
                psafe <- c(psafe, pleft[unlucky][j]) # species i is no longer likely to extinction
                
                N_enhance <- N_enhance + 1 # update
              }
            }
          }
        }else{ # if there is no plant likely to coextinction
          coext_plants <- c(coext_plants, pleft[unlucky]) # no plant is coextinct
        }
      }
      
      unlucky[psafe] <- FALSE # species are safe
      
      # Update outputs
      coext_animals <- c()
      coext_plants <- unique(coext_plants)
      plantNA[coext_plants] <- NA
      lostplants <- c(lostplants, coext_plants)
      imatrix[,coext_plants] <- 0
      orimatrix <- imatrix
      
      # Recalculate dependencies
      for(i in 1:npla){
        if(sum(imatrix[,i])==0){
          M[,i,1] <- 0
        }else{
          M[,i,1] <- imatrix[,i]/sum(imatrix[,i])
        }
      }
      
      if(length(coext_plants)>0){ # if any plant has gone coextinct
        degree <- degree + 1
        degree_when_lost_plants <- c(degree_when_lost_plants, rep(degree, length(coext_plants)))
        degree_table[degree,] <- data.frame(degree, "plant", length(coext_plants))
      }
    }
    
    equilibrium <- equilibrium + (length(coext_plants) + length(coext_animals))==0
  } # END Cascade loop
  
  
  #---- Calculating some outputs ----
  Av_distPlants <- mean(distPlants)
  Av_distAnim <- mean(distAnim)
  
  
  #---- Output ----
  return(list(degree_table=degree_table, N_likely=N_likely, N_newint=N_newint,
              N_enhance=N_enhance, N_weak=N_weak))
  
  
} # END stochinteract
