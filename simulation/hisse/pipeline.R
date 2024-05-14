#######################################
##            Chapter 1              ##
##   HiSSE + FBD Simulation Study    ##
##     Bruno do Rosario Petrucci     ##
#######################################

###
# goal: replicate a portion of Beaulieu & O'Meara 2016 HiSSE
# tests using trees including fossil samples as well (FBD)

### 
# packages required

# devtools 
suppressMessages(library(devtools))

# paleobuddy - need to remember to be in dev_traits branch
suppressMessages(library(paleobuddy))

# APE (for saving trees)
library(ape)

# readr (for write_tsv)
library(readr)

# load paleobuddy
load_all()

###
# simulation pipeline settings:
# Number of traits: 1 real, 3 neutral (q from 0.01 to 1)
# Trait evolution: HiSSE
  # states: c(0A, 0B, 1A, 1B)
  # X0: 0A
  # Transition rates: qjk = q01 = 0.01, and q10 = 0.005
# mu: mujk = mu0A = 0.03; mu0B = 0.06
# lambda: lambdajk = lambda0A = 0.1; lambda1B = 0.2
# psi: 0.01, 0.05, 0.1
# N: 100

###
# set up the baseline simulation settings

# initial number of species
n0 <- 1

# expected number of extant species
N <- 100

# base (null model) rates
# equivalent to BiSSE pipeline
# 0 = 0A, 1 = 1A, 2 = 0B, 3 = 1B
lambda0 <- lambda1 <- lambda2 <- lambda3 <- 0.1

mu0 <- mu1 <- mu2 <- mu3 <- 0.03

q01 <- q10 <- q23 <- q32 <- q02 <- q20 <- q13 <- q31 <- 0.01

q03 <- q30 <- q12 <- q21 <- 0

psi <- 0.05

# create variations we want to test -- also from BiSSE
lambda1Var <- lambda3Var <- 0.2

mu0Var <- mu2Var <- 0.06

q10Var <- q32Var <- 0.005

lambda3Var2 <- 0.2

mu3Var2 <- 0.06

###
# create directories for saving data

# make a smarter dir.create function
smart.dir.create <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir)
}

###
# consider parameter changes and run simulations

# collect baseline parameters in a data
# frame that will later become the key

key <- data.frame(ref = "base",
                  lambda0 = lambda0,
                  lambda1 = lambda1,
                  lambda2 = lambda2,
                  lambda3 = lambda3,
                  mu0 = mu0,
                  mu1 = mu1, 
                  mu2 = mu2,
                  mu3 = mu3,
                  q01 = q01,
                  q10 = q10,
                  q02 = q02,
                  q20 = q20,
                  q03 = q03,
                  q30 = q30,
                  q12 = q12,
                  q21 = q21,
                  q13 = q13,
                  q31 = q31,
                  q23 = q23,
                  q32 = q32,
                  psi = psi)

# make a vector to hold the baseline parameters
base <- key[1, ]

# reference vector for each parameter change
refs <- c("high_lambda1", "high_mu0", "low_q10", "high_lambda1B", "low_mu0B")

# the corresponding column numbers to change in the base
change <- list(c(3, 5), c(6, 8), c(11, 21), 5, 8)

# the new values
new <- list(c(lambda1Var, lambda3Var),  c(mu0Var, mu2Var), c(q10Var, q32Var),
         lambda3Var2, mu3Var2)

# copy the baseline a bunch of times to fill the key
key <- rbind(key, key[rep(1, length(change)), ])
rownames(key) <- 1:nrow(key)

# populate it with the values in aux
for (i in 2:nrow(key)) {
  key$ref[i] <- refs[i - 1]
  key[i, change[[i - 1]]] <- new[[i - 1]]
}

###
# auxiliary functions

# make a make.phylo output usable by FBD in RevBayes
# credit: Joshua Justison
collapseFossils <- function(tree) {
  # number of tips in the original tree
  nTips <- length(tree$tip.label)
  
  # indices to delete - edges with length 0
  delInd <- tree$edge.length == 0
  
  # tips to delete
  deletedTips<- tree$edge[delInd, 2]
  
  # get the labels of the deleted tips
  internalLabels <- tree$tip.label[deletedTips]
  
  # node numbers for the tips we are deleting
  nodeNumbers <- tree$edge[delInd, 1]
  
  # drop tips, but keep them as nodes
  tree <- drop.tip(tree, deletedTips, collapse.singles = FALSE)
  
  # reset node labels
  tree$node.label <- rep("", rep(tree$Nnode))
  
  # add previous deleted tips' labels to the new nodes
  tree$node.label[nodeNumbers - nTips] <- internalLabels
  
  return(tree)
}

# create simulation function for one rep
simulate_rep <- function(lambda0, lambda1, lambda2, lambda3,
                         mu0, mu1, mu2, mu3, 
                         q01, q10, q02, q20, q03, q30,
                         q12, q21, q13, q31, q23, q32, N, psi) {
  ### set up parameters
  # initial number of species
  n0 <- 1
  
  # speciation
  lambda <- c(lambda0, lambda1, lambda2, lambda3)
  
  # extinction
  mu <- c(mu0, mu1, mu2, mu3)
  
  # transition rate matrix
  Q <- list(matrix(c(0, q10, q20, q30, q01, 0, q21, q31, 
                     q02, q12, 0, q32, q03, q13, q23, 0), 4, 4),
            matrix(c(0, 0.01, 0.01, 0, 0.01, 0, 0, 0.01, 
                     0.01, 0, 0, 0.01, 0, 0.01, 0.01, 0), 4, 4),
            matrix(c(0, 0.1, 0.1, 0, 0.1, 0, 0, 0.1, 
                     0.1, 0, 0, 0.1, 0, 0.1, 0.1, 0), 4, 4),
            matrix(c(0, 1, 1, 0, 1, 0, 0, 1, 
                     1, 0, 0, 1, 0, 1, 1, 0), 4, 4))
  
  # initialize check
  bounds <- FALSE
  
  # counter to make sure we don't keep insisting on bad parameter values
  counter <- 1
  ratioSampled <- c()
  ratioTrait <- c()
  
  while (!bounds) {
    # run simulation
    sim <- bd.sim.musse(n0, lambda, mu, "number", N = N, 
                        nTraits = 4, nStates = 4, Q = Q)
    
    # run again until we have some extinct species
    while (sum(!sim$SIM$EXTANT) == 0) {
      sim <- bd.sim.musse(n0, lambda, mu, "number", N = N, 
                          nTraits = 4, nStates = 4, Q = Q)
      
    }
    
    # traits
    real <- lapply(sim$TRAITS, function(x) x[[1]])
    neutral_low <- lapply(sim$TRAITS, function(x) x[[2]])
    neutral_mid <- lapply(sim$TRAITS, function(x) x[[3]])
    neutral_high <- lapply(sim$TRAITS, function(x) x[[4]])
    
    # and get the sim object
    sim <- sim$SIM
  
    # get tMax
    tMax <- max(sim$TS)
    
    # sample data frame
    sample <- data.frame()
    while (nrow(sample) == 0) {
      sample <- suppressMessages(sample.clade(sim, psi, max(sim$TS)))
    }
    
    # get complete tree
    tree <- make.phylo(sim)
    
    # create an ultrametric tree
    treeUltra <- drop.fossil(tree)
    
    # start spec to false
    spec <- FALSE
    
    # while the tMRCA is not a speciation event
    while (!spec) {
      # get a sampled ancestor tree
      treeSA <- make.phylo(sim, sample)
      
      # and make an FBD tree
      # drop fossils
      treeFBD <- drop.tip(treeSA, tip = paste0("t", which(!sim$EXTANT)))
      
      # get a list of taxa for tree and SAtree
      fbdTaxa <- as.numeric(unlist(lapply(treeFBD$tip.label, 
                                          function(x) sub("t", "", x))))
      ultraTaxa <- as.numeric(unlist(lapply(treeUltra$tip.label, 
                                            function(x) sub("t", "", x))))
      
      # and make 0 length branches nodes
      treeFBD <- collapseFossils(treeFBD)
      
      # age of most recent common ancestor (MRCA)
      tMRCA <- max(node.depth.edgelength(treeFBD))
      
      # whether tMRCA is a speciation
      spec <- any(tMRCA == sim$TS)
      
      # if spec is not true, remove sample corresponding to tMRCA
      if (!spec) {
        sample <- sample[sample$SampT != tMRCA, ]
      }
    }
    
    # sampled also including the extant species
    sampled <- unlist(lapply(1:length(sim$TS), function(x) 
      x %in% (unique(floor(as.numeric(sub("t", "", sample$Species))))) ||
        x %in% which(sim$EXTANT)))
    
    # check how many extinct species were sampled
    nSampled <- sum(sampled & !sim$EXTANT)
    
    # get percentage of species sampled
    ratioSampled <- nSampled/sum(!sim$EXTANT)
    
    # times to sample traits from the saTree taxa
    treeFBDTraitT <- 
      unlist(lapply(fbdTaxa, function(x) {
        # part of sample that we care about
        sampleCut <- sample[sample$Species == paste0("t", floor(x)), 
                            "SampT"]
        
        # get the time it was sampled
        ifelse(x == floor(x), 0, 
               sampleCut[round(10^ceiling(log(length(sampleCut) + 1, 10))*
                                 (x - floor(x)))])}))
    
    # name it
    names(treeFBDTraitT) <- paste0("t", fbdTaxa)
    
    # and make the trait lists
    # for the ultrametric tree
    # and make traits hidden
    treeUltraReal <- 
      as.numeric(unlist(lapply(ultraTaxa, function(x)
        tail(real[[x]]$value, 1))))
    treeUltraReal[treeUltraReal == 2] <- 0
    treeUltraReal[treeUltraReal == 3] <- 1
    
    treeUltraNeutralLow <- 
      as.numeric(unlist(lapply(ultraTaxa, function(x)
        tail(neutral_low[[x]]$value, 1))))
    treeUltraNeutralLow[treeUltraNeutralLow == 2] <- 0
    treeUltraNeutralLow[treeUltraNeutralLow == 3] <- 1
    
    treeUltraNeutralMid <- 
      as.numeric(unlist(lapply(ultraTaxa, function(x)
        tail(neutral_mid[[x]]$value, 1))))
    treeUltraNeutralMid[treeUltraNeutralMid == 2] <- 0
    treeUltraNeutralMid[treeUltraNeutralMid == 3] <- 1
    
    treeUltraNeutralHigh <- 
      as.numeric(unlist(lapply(ultraTaxa, function(x)
        tail(neutral_high[[x]]$value, 1))))
    treeUltraNeutralHigh[treeUltraNeutralHigh == 2] <- 0
    treeUltraNeutralHigh[treeUltraNeutralHigh == 3] <- 1
    
    names(treeUltraReal) <- names(treeUltraNeutralLow) <- 
      names(treeUltraNeutralMid) <- names(treeUltraNeutralHigh) <- 
      paste0("t", ultraTaxa)
    
    # function to extract a trait list
    traitList <- function(traits, x) {
      # get traits df
      tr <-  traits[[floor(x)]]
      
      # and time to sample
      trT <- treeFBDTraitT[which(fbdTaxa == x)]
      
      # extract specific time
      tr$value[tr$min <= trT & tr$max > trT]
    }
    
    # and the fbd tree
    treeFBDReal <- 
      as.numeric(unlist(lapply(fbdTaxa, function(x) traitList(real, x))))
    treeFBDReal[treeFBDReal == 2] <- 0
    treeFBDReal[treeFBDReal == 3] <- 1
    
    treeFBDNeutralLow <- 
      as.numeric(unlist(lapply(fbdTaxa, function(x) traitList(neutral_low, x))))
    treeFBDNeutralLow[treeFBDNeutralLow == 2] <- 0
    treeFBDNeutralLow[treeFBDNeutralLow == 3] <- 1
    
    treeFBDNeutralMid <- 
      as.numeric(unlist(lapply(fbdTaxa, function(x) traitList(neutral_mid, x))))
    treeFBDNeutralMid[treeFBDNeutralMid == 2] <- 0
    treeFBDNeutralMid[treeFBDNeutralMid == 3] <- 1
    
    treeFBDNeutralHigh <- 
      as.numeric(unlist(lapply(fbdTaxa, function(x) traitList(neutral_high, x))))
    treeFBDNeutralHigh[treeFBDNeutralHigh == 2] <- 0
    treeFBDNeutralHigh[treeFBDNeutralHigh == 3] <- 1
    
    names(treeFBDReal) <- names(treeFBDNeutralLow) <- 
      names(treeFBDNeutralMid) <- names(treeFBDNeutralHigh) <- 
      paste0("t", fbdTaxa)
    
    # check the minimum percentage of the rare trait
    ratioTrait <- min(sum(treeUltraReal)/length(ultraTaxa),
                      sum(treeFBDReal)/length(fbdTaxa))
    
    # checks - at least 50% of extinct species sampled, and
    # both traits represented at least 10%
    bounds <- (ratioTrait >= 0.1) &&
      (ratioTrait <= 0.9)
    
    # increase counter
    counter <- counter + 1

    # if counter is higher than 10, maybe rethink the parameters
    if (counter > 100) {
      print(ratioSampled)
      print(ratioTrait)
      stop("Hard to find replicate within bounds")
    }
  }
  
  # make a list to return
  res <- list(SIM = sim, SAMPLE = sample,
              TREE = tree, ULTRATREE = treeUltra, FBDTREE = treeFBD,
              REALULTRA = treeUltraReal, REALFBD = treeFBDReal,
              NEUTRALLOWULTRA = treeUltraNeutralLow, 
              NEUTRALLOWFBD = treeFBDNeutralLow,
              NEUTRALMIDULTRA = treeUltraNeutralMid, 
              NEUTRALMIDFBD = treeFBDNeutralMid,
              NEUTRALHIGHULTRA = treeUltraNeutralHigh, 
              NEUTRALHIGHFBD = treeFBDNeutralHigh)
  
  return(res)
}

# function to save simulation lists
save_sims <- function(nReps, simReps, targetDir) {
  # extract simulations
  simList <- lapply(1:nReps, function(x) simReps[[x]]$SIM)
  
  # trait data 
  realUltra <- lapply(1:nReps, function(x) simReps[[x]]$REALULTRA)
  realFBD <- lapply(1:nReps, function(x) simReps[[x]]$REALFBD)
  neutralLowUltra <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALLOWULTRA)
  neutralLowFBD <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALLOWFBD)
  neutralMidUltra <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALMIDULTRA)
  neutralMidFBD <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALMIDFBD)
  neutralHighUltra <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALMIDULTRA)
  neutralHighFBD <- lapply(1:nReps, function(x) simReps[[x]]$NEUTRALMIDFBD)
  
  # fossil records
  fossilsList <- lapply(1:nReps, function(x) simReps[[x]]$SAMPLE)
  
  # trees
  treeList <-  lapply(1:nReps, function(x) simReps[[x]]$TREE)
  treeUltraList <- lapply(1:nReps, function(x) simReps[[x]]$ULTRATREE)
  treeFBDList <- lapply(1:nReps, function(x) simReps[[x]]$FBDTREE)
  
  # save simulation RData
  save(simList, file = paste0(targetDir, "sim_list.RData"))
  
  # simulations directory
  simsDir <- paste0(targetDir, "sims/")
  smart.dir.create(simsDir)
  
  # trait lists directories
  traitsDir <- paste0(targetDir, "traits/")
  smart.dir.create(traitsDir)
  
  traitsUltraDir <- paste0(traitsDir, "ultrametric/")
  smart.dir.create(traitsUltraDir)
  
  neutralUltraDir <- paste0(traitsUltraDir, "neutral/")
  smart.dir.create(neutralUltraDir)
  
  neutralUltraLowDir <- paste0(neutralUltraDir, "low/")
  smart.dir.create(neutralUltraLowDir)
  
  neutralUltraMidDir <- paste0(neutralUltraDir, "mid/")
  smart.dir.create(neutralUltraMidDir)

  neutralUltraHighDir <- paste0(neutralUltraDir, "high/")
  smart.dir.create(neutralUltraHighDir)
  
  traitsFBDDir <- paste0(traitsDir, "fbd/")
  smart.dir.create(traitsFBDDir)
  
  neutralFBDDir <- paste0(traitsFBDDir, "neutral/")
  smart.dir.create(neutralFBDDir)
  
  neutralFBDLowDir <- paste0(neutralFBDDir, "low/")
  smart.dir.create(neutralFBDLowDir)
  
  neutralFBDMidDir <- paste0(neutralFBDDir, "mid/")
  smart.dir.create(neutralFBDMidDir)
  
  neutralFBDHighDir <- paste0(neutralFBDDir, "high/")
  smart.dir.create(neutralFBDHighDir)
  
  # fossils directory
  fossilsDir <- paste0(targetDir, "fossils/")
  smart.dir.create(fossilsDir)
  
  # and tree directories
  treeDir <- paste0(targetDir, "trees/")
  smart.dir.create(treeDir)
  
  treeUltraDir <- paste0(treeDir, "ultrametric/")
  smart.dir.create(treeUltraDir)
  
  treeFBDDir <- paste0(treeDir, "fbd/")
  smart.dir.create(treeFBDDir)
  
  # print simulations to file
  invisible(lapply(1:nReps, function(x)
    capture.output(print(simList[[x]]),
      file = paste0(simsDir, "sim_", x, ".txt")))) 
  invisible(lapply(1:nReps, function(x)
    capture.output(lapply(simList[[x]], function(x) {
      cat("\n")
      x}),
      file = paste0(simsDir, "sim_", x, ".txt"),
      append = TRUE))) 
  
  # save trait data as as .nex
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(realUltra[[x]], 
                     file = paste0(traitsUltraDir, "real_ultra_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralLowUltra[[x]], 
                     file = paste0(neutralUltraLowDir, 
                                   "neutral_low_ultra_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralMidUltra[[x]], 
                     file = paste0(neutralUltraMidDir, 
                                   "neutral_mid_ultra_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralHighUltra[[x]], 
                     file = paste0(neutralUltraHighDir, 
                                   "neutral_high_ultra_", x, ".nex"),
                     format = "standard")))
  
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(realFBD[[x]], 
                     file = paste0(traitsFBDDir, "real_fbd_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralLowFBD[[x]], 
                     file = paste0(neutralFBDLowDir, 
                                   "neutral_low_fbd_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralMidFBD[[x]], 
                     file = paste0(neutralFBDMidDir, 
                                   "neutral_mid_fbd_", x, ".nex"),
                     format = "standard")))
  invisible(lapply(1:nReps, function(x)
    write.nexus.data(neutralHighFBD[[x]], 
                     file = paste0(neutralFBDHighDir, 
                                   "neutral_high_fbd_", x, ".nex"),
                     format = "standard")))
  
  # save fossil records as tsv
  invisible(lapply(1:nReps, function(x)
    write_tsv(fossilsList[[x]], 
              file = paste0(fossilsDir, "fossils_", x, ".tsv"))))
  
  # save ultrametric trees
  invisible(lapply(1:nReps, function(x)
    write.nexus(treeUltraList[[x]], 
                file = paste0(treeUltraDir, "tree_ultra_", x, ".nex"))))
  # and fbd trees
  invisible(lapply(1:nReps, function(x)
    write.nexus(treeFBDList[[x]], 
                file = paste0(treeFBDDir, "tree_fbd_", x, ".nex"))))

  # save complete trees as well
  invisible(lapply(1:nReps, function(x)
    write.nexus(treeList[[x]], 
                file = paste0(treeDir, "tree_", x, ".nex"))))
}

# create function to run simulations for a list of parameters
simulate <- function(seeds, nReps, comb, key, simDir, N, psi) {
  ### recover parameters from key
  # ref
  ref <- key$ref[comb]
  
  # speciation rates
  lambda0 <- key$lambda0[comb]
  lambda1 <- key$lambda1[comb]
  lambda2 <- key$lambda2[comb]
  lambda3 <- key$lambda3[comb]
  
  # extinction rates
  mu0 <- key$mu0[comb]
  mu1 <- key$mu1[comb]
  mu2 <- key$mu2[comb]
  mu3 <- key$mu3[comb]
  
  # transition rates
  q01 <- key$q01[comb]
  q10 <- key$q10[comb]
  q02 <- key$q02[comb]
  q20 <- key$q20[comb]
  q03 <- key$q03[comb]
  q30 <- key$q30[comb]
  q12 <- key$q12[comb]
  q21 <- key$q21[comb]
  q13 <- key$q13[comb]
  q31 <- key$q31[comb]
  q23 <- key$q23[comb]
  q32 <- key$q32[comb]
  
  # fossil sampling rates
  psi <- psi
  
  ## create directories
  # base directory for simulations, if it wasn't created
  smart.dir.create(simDir)
  
  # create directory for this combination of simulations
  combDir <- paste0(simDir, comb, "_", ref, "/")
  smart.dir.create(combDir)
  
  # run simulations
  simReps <- lapply(1:nReps, function(x) {
    print(paste0("psi: ", psi, " comb: ", ref, " rep: ", x, 
                 " seed: ", seeds[x]))
    set.seed(seeds[x])
    simulate_rep(lambda0, lambda1, lambda2, lambda3, 
                 mu0, mu1, mu2, mu3,
                 q01, q10, q02, q20, q03, q30,
                 q12, q21, q13, q31, q23, q32, N, psi)
  })
  
  # save them
  invisible(save_sims(nReps, simReps, combDir))
}

###
# number of reps to run each combination of parameters
nReps <- 100

# psi ref
psiRefs <- c("1_low_psi", "2_med_psi", "3_high_psi")

# psi vector
psiVar <- c(0.01, 0.05, 0.1)

# for each psi
for (i in 1:length(psiVar)) {
  # simulations directory
  simDir <- paste0("/Users/petrucci/Documents/research/ssetests_chap1/", 
                   "simulation/hisse/replicates/", psiRefs[i], "/")
  smart.dir.create(simDir)
  
  # run simulations for each combination of parameters
  for (comb in 1:nrow(key)) {
    # get a seed for each rep
    seeds <- runif(nReps, (comb - 1)*nReps, comb*nReps)
    
    # create directory for this combination of simulations
    combDir <- paste0(simDir, comb, "_", key$ref[comb], "/")
    smart.dir.create(combDir)
  
    # save seeds
    save(seeds, file = paste0(combDir, "seeds.RData"))
    
    # simulate the reps with this parameter combination
    simulate(seeds, nReps, comb, key, simDir, N, psiVar[i])
  }
}
