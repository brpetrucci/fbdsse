#####################################
##            Chapter 1            ##
##   BiSSE + FBD Simulation Study  ##
##            Plotting             ##
##    Bruno do Rosario Petrucci    ##
#####################################

### 
## packages

# ggplot2
library(ggplot2)
library(glue)
library(ggridges)
library(forcats)
library(reshape2)

# coda
library(coda)

###
## reading logs

# create list of logs
logs <- list()

# create a list of ESS values
ess <- list()

# and create a dataframe to hold comb values for each ref
refs_df <- data.frame(matrix(nrow = 0, ncol = 4))

# number of parameter combinations
n_refs <- 76

# number of reps
n_reps <- 100

# vectors of combination names
psiRefs <- c("1_low_psi", "2_med_psi", "3_high_psi")
parRefs <- c("1_base", "2_high_lambda1", "3_high_mu0", "4_low_q10")
modelRefs <- c("fbd", "bisse", "both")
traitRefs <- c("real", "low", "mid", "high")

# base directory
base_dir <- "/Users/petrucci/Documents/research/ssetests_chap1/low_sample_size/"

# for each combination
for (i in 1:n_refs) {
  # source the ref in question
  source(paste0(base_dir, "analysis/bisse/refs/refs_", i, ".Rev"))
  
  # if model is 1 (fbd), trait is not necessary, same for 2 (bisse) and psi
  if (modelComb == 1) {
    traitComb <- 1
  } else if (modelComb == 2) {
    psiComb <- 1
  }
  
  # psi reference
  psiRef <- psiRefs[psiComb]
  
  # parameter reference
  parRef <- parRefs[parComb]
  
  # model reference
  modelRef <- modelRefs[modelComb]
  
  # trait reference
  traitRef <- traitRefs[traitComb]
  
  # append to refs_df
  refs_df <- rbind(refs_df, c(psiComb, parComb, modelComb, traitComb))
  
  # if model is 1, no traitRef field is necessary
  if (modelComb == 1) {
    # get the directory where the reps are
    target_dir <- paste0(base_dir, "output/bisse/", psiRef, "/", parRef, "/",
                         modelRef, "/")
  } else {
    # get the directory where the reps are
    target_dir <- paste0(base_dir, "output/bisse/", psiRef, "/", parRef, "/",
                         modelRef, "/", traitRef, "/")
  }
  
  # get each rep
  reps <- lapply(1:n_reps, function(x) {
    table <- read.delim(paste0(target_dir, "rep_", x, "/bisse.log"), 
                        header = TRUE, sep = "\t");
    table[, 5:ncol(table)]
    })
  
  # get the ess vales
  ess_reps <- lapply(1:n_reps, function(x) effectiveSize(reps[[x]]))
  col_names <- names(ess_reps[[1]])
  ess_reps <- do.call(rbind.data.frame, ess_reps)
  colnames(ess_reps) <- col_names
  
  # put it in the list of logs
  logs[[i]] <- reps
  ess[[i]] <- ess_reps
}

# name refs
colnames(refs_df) <- c("psiComb", "parComb", "modelComb", "traitComb")

###
## testing accuracy - FBD

# get log indices for the FBD reps
fbd_refs <- refs_df[which(refs_df[, 3] == 1), ]
fbd_refs$comb <- as.numeric(rownames(fbd_refs))

# add true values to the data frame
fbd_refs$lambda1 <- 0.1
fbd_refs$lambda2 <- c(0.1, 0.2)[(fbd_refs$parComb == 2) + 1]
fbd_refs$mu1 <- c(0.03, 0.06)[(fbd_refs$parComb == 3) + 1]
fbd_refs$mu2 <- 0.03
fbd_refs$psi <- c(0.01, 0.05, 0.1)[fbd_refs$psiComb]

# make data frames for 95% CI and median, and one for mean
fbd_low <- fbd_med <- fbd_high <- fbd_mean <- 
  data.frame(matrix(nrow = 0, ncol = 4))
fbd_cov <- data.frame(matrix(nrow = 0, ncol = 5))

# iterate through logs
for (i in 1:nrow(fbd_refs)) {
  # get ref
  ref <- fbd_refs$comb[i]
  
  # and reps
  for (j in 1:n_reps) {
    # get log
    log <- logs[[ref]][[j]]
    
    # apply burnout
    log <- log[(nrow(log)/5):nrow(log), ]
    
    # get quantiles
    quants <- do.call(rbind.data.frame, lapply(c(0.025, 0.5, 0.975), 
                                                    function(x)
      unlist(lapply(1:ncol(log), function(y) quantile(log[, y], x)))))
    quants <- cbind(rep(ref, 3), quants)
    colnames(quants) <- c("comb", "lambda", "mu", "psi")
    
    # fill data frames
    fbd_low <- rbind(fbd_low, quants[1, ])
    fbd_med <- rbind(fbd_med, quants[2, ])
    fbd_high <- rbind(fbd_high, quants[3, ])
    fbd_mean <- rbind(fbd_mean, c(ref, colMeans(log)))
  }
  
  # low and high dfs for this combination
  fbd_low_i <- fbd_low[fbd_low[, 1] == ref, ]
  fbd_high_i <- fbd_high[fbd_high[, 1] == ref, ]
  
  # coverage
  fbd_cov <- rbind(fbd_cov, 
                   c(sum(fbd_low_i[, 2] < fbd_refs$lambda1[i] & 
                     fbd_high_i[, 2] > fbd_refs$lambda1[i]),
                   sum(fbd_low_i[, 2] < fbd_refs$lambda2[i] &
                     fbd_high_i[, 2] > fbd_refs$lambda2[i]),
                   sum(fbd_low_i[, 3] < fbd_refs$mu1[i] &
                     fbd_high_i[, 3] > fbd_refs$mu1[i]),
                   sum(fbd_low_i[, 3] < fbd_refs$mu2[i] &
                     fbd_high_i[, 3] > fbd_refs$mu2[i]),
                   sum(fbd_low_i[, 4] < fbd_refs$psi[i] &
                     fbd_high_i[, 4] > fbd_refs$psi[i])))
  colnames(fbd_cov) <- c("lambda1", "lambda2", "mu1", "mu2", "psi")
}

# normalize and name fbd_cov
fbd_cov <- fbd_cov / n_reps
rownames(fbd_cov) <- rownames(fbd_refs)

# name fbd_mean
colnames(fbd_mean) <- colnames(fbd_low)

# mean boxplot
fbd_mean_bplot_lambda <- ggplot(fbd_mean, aes(factor(comb), lambda)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = fbd_refs, aes(x = factor(comb), y = lambda1), 
             col = "red", size = 2) + 
  geom_point(data = fbd_refs, aes(x = factor(comb), y = lambda2), 
             col = "red", size = 2) + 
  labs(x = "Parameter combination") +
  theme_bw()
fbd_mean_bplot_mu <- ggplot(fbd_mean, aes(factor(comb), mu)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = fbd_refs, aes(x = factor(comb), y = mu1), 
             col = "red", size = 2) + 
  geom_point(data = fbd_refs, aes(x = factor(comb), y = mu2), 
             col = "red", size = 2) + 
  ylim(0, 0.1) +
  labs(x = "Parameter combination") +
  theme_bw()
fbd_mean_bplot_psi <- 
  ggplot(fbd_mean, aes(reorder(factor(comb), psi, FUN = mean), psi)) + 
  geom_boxplot(outlier.shape = NA) + geom_point(data = fbd_refs,
                              aes(x = reorder(factor(comb), psi), y = psi),
                              col = "red", size = 2) + 
  labs(x = "Parameter combination") + 
  theme_bw()

###
## testing accuracy - BiSSE (real trait only)

# get log indices for the BiSSE reps with real traits
bisse_refs <- refs_df[which(refs_df[, 3] == 2 & refs_df[, 4] == 1), ]
bisse_refs$comb <- as.numeric(rownames(bisse_refs))

# add true values to the data frame
bisse_refs$lambda1 <- 0.1
bisse_refs$lambda2 <- c(0.1, 0.2)[(bisse_refs$parComb == 2) + 1]
bisse_refs$mu1 <- c(0.03, 0.06)[(bisse_refs$parComb == 3) + 1]
bisse_refs$mu2 <- 0.03
bisse_refs$q01 <- 0.01
bisse_refs$q10 <- c(0.01, 0.005)[(bisse_refs$parComb == 4) + 1]

# make data frames for 95% CI and median, one for mean, and one for coverage
bisse_low <- bisse_med <- bisse_high <- bisse_mean <- bisse_cov <- bisse_pce <-
  data.frame(matrix(nrow = 0, ncol = 6))

# iterate through logs
for (i in 1:nrow(bisse_refs)) {
  # get ref
  ref <- bisse_refs$comb[i]
  
  # and reps
  for (j in 1:n_reps) {
    # get log
    log <- logs[[ref]][[j]]
    
    # apply burnout and take out pi
    log <- log[(nrow(log)/5):nrow(log), -c(5, 6)]
    
    # get quantiles
    quants <- do.call(rbind.data.frame, lapply(c(0.025, 0.5, 0.975), 
                                               function(x)
      unlist(lapply(1:ncol(log), function(y) quantile(log[, y], x)))))
    quants <- cbind(rep(ref, 3), quants)
    colnames(quants) <- c("comb", "lambda1", "lambda2", 
                          "mu1", "mu2", "q01", "q10")
    
    # fill data frames
    bisse_low <- rbind(bisse_low, quants[1, ])
    bisse_med <- rbind(bisse_med, quants[2, ])
    bisse_high <- rbind(bisse_high, quants[3, ])
    bisse_mean <- rbind(bisse_mean, c(ref, colMeans(log)))
  }
  
  
  # low and high dfs for this combination
  bisse_low_i <- bisse_low[bisse_low[, 1] == ref, ]
  bisse_high_i <- bisse_high[bisse_high[, 1] == ref, ]
  bisse_mean_i <- bisse_mean[bisse_mean[, 1] == ref, ]
  
  # coverage
  bisse_cov <- rbind(bisse_cov, 
                   c(sum(bisse_low_i[, 2] < bisse_refs$lambda1[i] & 
                           bisse_high_i[, 2] > bisse_refs$lambda1[i]),
                     sum(bisse_low_i[, 3] < bisse_refs$lambda2[i] &
                           bisse_high_i[, 3] > bisse_refs$lambda2[i]),
                     sum(bisse_low_i[, 4] < bisse_refs$mu1[i] &
                           bisse_high_i[, 4] > bisse_refs$mu1[i]),
                     sum(bisse_low_i[, 5] < bisse_refs$mu2[i] &
                           bisse_high_i[, 5] > bisse_refs$mu2[i]),
                     sum(bisse_low_i[, 6] < bisse_refs$q01[i] &
                           bisse_high_i[, 6] > bisse_refs$q01[i]),
                     sum(bisse_low_i[, 7] < bisse_refs$q10[i] &
                           bisse_high_i[, 7] > bisse_refs$q10[i])))
  colnames(bisse_cov) <- c("lambda1", "lambda2", "mu1", "mu2", "q01", "q10")
  
  # percent error
  bisse_pce <- rbind(bisse_pce,
                     abs(c(mean((bisse_mean_i[, 2] - bisse_refs$lambda1[i]) / 
                                  bisse_refs$lambda1[i]),
                       mean((bisse_mean_i[, 3] - bisse_refs$lambda2[i]) / 
                              bisse_refs$lambda2[i]),
                       mean((bisse_mean_i[, 4] - bisse_refs$mu1[i]) / 
                              bisse_refs$mu1[i]),
                       mean((bisse_mean_i[, 5] - bisse_refs$mu2[i]) / 
                              bisse_refs$mu2[i]),
                       mean((bisse_mean_i[, 6] - bisse_refs$q01[i]) /
                              bisse_refs$q01[i]),
                       mean((bisse_mean_i[, 7] - bisse_refs$q10[i]) /
                              bisse_refs$q10[i]))))
  colnames(bisse_pce) <- colnames(bisse_cov)
}

# normalize and name bisse_cov
bisse_cov <- bisse_cov / n_reps
rownames(bisse_cov) <- rownames(bisse_pce) <- rownames(bisse_refs)

# name bisse_mean
colnames(bisse_mean) <- colnames(bisse_low)

# round PCE
bisse_pce <- round(bisse_pce, digits = 3)

# get mean precision per comb
precision_q01 <- unlist(lapply(unique(bisse_mean$comb), 
                               function(x) 
                                 mean(bisse_high$q01[bisse_high$comb == x] - 
                                        bisse_low$q01[bisse_low$comb == x])))
names(precision_q01) <- unique(bisse_mean$comb)
precision_q10 <- unlist(lapply(unique(bisse_mean$comb), 
                               function(x) 
                                 mean(bisse_high$q10[bisse_high$comb == x] - 
                                        bisse_low$q10[bisse_low$comb == x])))
names(precision_q10) <- unique(bisse_mean$comb)

# mean boxplot
bisse_mean_bplot_lambda1 <- ggplot(bisse_mean, aes(factor(comb), lambda1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = lambda1), 
             col = "red", size = 2) +
  labs(x = "Parameter combination") +
  theme_bw()
bisse_mean_bplot_lambda2 <- ggplot(bisse_mean, aes(factor(comb), lambda2)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = lambda2), 
             col = "red", size = 2) +
  labs(x = "Parameter combination") +
  theme_bw()
bisse_mean_bplot_mu1 <- ggplot(bisse_mean, aes(factor(comb), mu1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = mu1), 
             col = "red", size = 2) + 
  ylim(0, 0.07) +
  labs(x = "Parameter combination") +
  theme_bw()
bisse_mean_bplot_mu2 <- ggplot(bisse_mean, aes(factor(comb), mu2)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = mu2), 
             col = "red", size = 2) + 
  ylim(0, 0.061) +
  labs(x = "Parameter combination") +
  theme_bw()
bisse_mean_bplot_q01 <- ggplot(bisse_mean, aes(factor(comb), q01)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = q01), 
             col = "red", size = 2) + 
  labs(x = "Parameter combination") +
  theme_bw()
bisse_mean_bplot_q10 <- ggplot(bisse_mean, aes(factor(comb), q10)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = bisse_refs, aes(x = factor(comb), y = q10), 
             col = "red", size = 2) + 
  labs(x = "Parameter combination") +
  theme_bw()

# ridge plots
bisse_ridge_lambda1 <- ggplot(bisse_mean[bisse_mean$comb %in% c(7, 11),], 
                              aes(x = lambda1, y = comb, 
                                  fill = factor(comb == 11),
                                  color = factor(comb == 11))) +
  geom_density_ridges() +
  geom_vline(xintercept = 0.1, color = "black", lwd = 1) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 2),
                    values = c("#E69F00", "#56B4E9")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  labs(title = expression("Mean "*lambda['0']*" estimates"), 
       x = expression("Estimated "*lambda['0'])) +
  geom_text(x = 0.04, y = 10, 
            label = paste0("Cov = ", bisse_cov$lambda1[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.04, y = 16, 
            label = paste0("Cov = ", bisse_cov$lambda1[2]),
            color = "#56B4E9", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.04, y = 9, 
            label = paste0("MPE = ", bisse_pce$lambda1[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.04, y = 15, 
            label = paste0("MPE = ", bisse_pce$lambda1[2]),
            color = "#56B4E9", size = 5, check_overlap = TRUE) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

bisse_ridge_lambda2 <- ggplot(bisse_mean[bisse_mean$comb %in% c(7, 11),], 
                              aes(x = lambda2, y = comb, 
                                  fill = factor(comb == 11),
                                  color = factor(comb == 11))) +
  geom_density_ridges() +
  geom_vline(xintercept = c(0.1, 0.2), color = c("#E69F00", "#56B4E9"), 
             lwd = c(1, 1)) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 2),
                    values = c("#E69F00", "#56B4E9")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  labs(title = expression("Mean "*lambda['1']*" estimates"), 
       x = expression("Estimated "*lambda['1'])) +
  geom_text(x = 0.28, y = 10, 
            label = paste0("Cov = ", bisse_cov$lambda2[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.28, y = 15, 
            label = paste0("Cov = ", bisse_cov$lambda2[2]),
            color = "#56B4E9", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.28, y = 9, 
            label = paste0("MPE = ", bisse_pce$lambda2[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.28, y = 14, 
            label = paste0("MPE = ", bisse_pce$lambda2[2]),
            color = "#56B4E9", size = 5, check_overlap = TRUE) +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

bisse_ridge_mu1 <- ggplot(bisse_mean[bisse_mean$comb %in% c(7, 15),], 
                              aes(x = mu1, y = comb, 
                                  fill = factor(comb == 15),
                                  color = factor(comb == 15))) +
  geom_density_ridges() +
  geom_vline(xintercept = c(0.03, 0.06), color = c("#E69F00", "#009E73"), 
             lwd = c(1, 1)) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 2),
                    values = c("#E69F00", "#009E73")) + 
  scale_color_manual(values = c("#E69F00", "#009E73")) +
  labs(title = expression("Mean "*mu['0']*" estimates"), 
       x = expression("Estimated "*mu['0'])) +
  geom_text(x = 0.1, y = 14, 
            label = paste0("Cov = ", bisse_cov$mu1[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.1, y = 24, 
            label = paste0("Cov = ", bisse_cov$mu1[3]),
            color = "#009E73", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.1, y = 13, 
            label = paste0("MPE = ", bisse_pce$mu1[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.1, y = 23, 
            label = paste0("MPE = ", bisse_pce$mu1[3]),
            color = "#009E73", size = 5, check_overlap = TRUE) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

bisse_ridge_mu2 <- ggplot(bisse_mean[bisse_mean$comb %in% c(7, 15),], 
                              aes(x = mu2, y = comb, 
                                  fill = factor(comb == 15),
                                  color = factor(comb == 15))) +
  geom_density_ridges() +
  geom_vline(xintercept = 0.03, color = "black", lwd = 1) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 3),
                    values = c("#E69F00", "#009E73")) + 
  scale_color_manual(values = c("#E69F00", "#009E73")) +
  labs(title = expression("Mean "*mu['1']*" estimates"), 
       x = expression("Estimated "*mu['1'])) +
  geom_text(x = 0.08, y = 14, 
            label = paste0("Cov = ", bisse_cov$mu2[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.08, y = 26, 
            label = paste0("Cov = ", bisse_cov$mu2[3]),
            color = "#009E73", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.08, y = 13, 
            label = paste0("MPE = ", bisse_pce$mu2[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.08, y = 25, 
            label = paste0("MPE = ", bisse_pce$mu2[3]),
            color = "#009E73", size = 5, check_overlap = TRUE) +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

bisse_ridge_q01 <- ggplot(bisse_mean[bisse_mean$comb %in% c(7, 19),], 
                          aes(x = q01, y = comb, 
                              fill = factor(comb == 19),
                              color = factor(comb == 19))) +
  geom_density_ridges() +
  geom_vline(xintercept = 0.01, color = "black", lwd = 1) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 4),
                    values = c("#E69F00", "#CC79A7")) + 
  scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  labs(title = expression("Mean "*q['01']*" estimates"), 
       x = expression("Estimated "*q['01'])) +
  geom_text(x = 0.025, y = 16, 
            label = paste0("Cov = ", bisse_cov$q01[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.025, y = 38, 
            label = paste0("Cov = ", bisse_cov$q01[4]),
            color = "#CC79A7", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.025, y = 14.5, 
            label = paste0("MPE = ", bisse_pce$q01[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.025, y = 36.5, 
            label = paste0("MPE = ", bisse_pce$q01[4]),
            color = "#CC79A7", size = 5, check_overlap = TRUE) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

bisse_ridge_q10 <- ggplot(bisse_mean[bisse_mean$comb %in% c(7, 19),], 
                          aes(x = q10, y = comb, 
                              fill = factor(comb == 19),
                              color = factor(comb == 19))) +
  geom_density_ridges() +
  geom_vline(xintercept = c(0.01, 0.005), color = c("#E69F00", "#CC79A7"), 
             lwd = c(1, 1)) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 4),
                    values = c("#E69F00", "#CC79A7")) + 
  scale_color_manual(values = c("#E69F00", "#CC79A7")) +
  labs(title = expression("Mean "*q['10']*" estimates"), 
       x = expression("Estimated "*q['10'])) +
  geom_text(x = 0.03, y = 16, 
            label = paste0("Cov = ", bisse_cov$q10[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.03, y = 38, 
            label = paste0("Cov = ", bisse_cov$q10[4]),
            color = "#CC79A7", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.03, y = 14.5, 
            label = paste0("MPE = ", bisse_pce$q10[1]),
            color = "#E69F00", size = 5, check_overlap = TRUE) +
  geom_text(x = 0.03, y = 36.5, 
            label = paste0("MPE = ", bisse_pce$q10[4]),
            color = "#CC79A7", size = 5, check_overlap = TRUE) +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# comparison plots Maddison 2007
bisse_comp_lambda <- ggplot(bisse_mean[bisse_mean$comb %in% c(7, 11),], 
                            aes(lambda1, lambda2, 
                                color = factor(comb == 11))) +
  geom_point() +
  geom_hline(yintercept = c(0.1, 0.2), color = c("#E69F00", "#56B4E9")) +
  scale_color_manual(name = expression("True "*lambda['1']),
                     labels = c(0.1, 0.2),
                     values = c("#E69F00", "#56B4E9")) + 
  labs(title = "Speciation rate estimates", 
       x = expression("Estimated "*lambda['0']),
       y = expression("Estimated "*lambda['1'])) +
  theme_bw() + 
  theme(title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))


bisse_comp_mu <- ggplot(bisse_mean[bisse_mean$comb %in% c(7, 15),], 
                        aes(mu1, mu2, 
                            color = factor(comb == 15))) +
  geom_point() +
  geom_vline(xintercept = c(0.03, 0.06), color = c("#E69F00", "#56B4E9")) +
  scale_color_manual(name = expression("True "*mu['0']),
                     labels = c(0.03, 0.06),
                     values = c("#E69F00", "#56B4E9")) + 
  labs(title = "Extinction rate estimates", 
       x = expression("Estimated "*mu['0']),
       y = expression("Estimated "*mu['1'])) +
  theme_bw() + 
  theme(title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

bisse_comp_q <- ggplot(bisse_mean[bisse_mean$comb %in% c(7, 19),], 
                       aes(q01, q10, 
                           color = factor(comb == 19))) +
  geom_point() +
  geom_hline(yintercept = c(0.01, 0.005), color = c("#E69F00", "#56B4E9")) +
  scale_color_manual(name = expression("True "*q['10']),
                     labels = c(0.01, 0.005),
                     values = c("#E69F00", "#56B4E9")) + 
  labs(title = "Transition rate estimates", 
       x = expression("Estimated "*q['01']),
       y = expression("Estimated "*q['10'])) +
  theme_bw() + 
  theme(title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

###
## testing accuracy - BiSSE + FBD

# get log indices for the combined reps with real traits
both_refs <- refs_df[which((refs_df[, 3] == 3 | refs_df[, 3] == 2) & refs_df[, 4] == 1), ]
both_refs$comb <- as.numeric(rownames(both_refs))
both_refs$comb_seq <- 1:nrow(both_refs)

# add true values to the data frame
both_refs$lambda1 <- 0.1
both_refs$lambda2 <- c(0.1, 0.2)[(both_refs$parComb == 2) + 1]
both_refs$mu1 <- c(0.03, 0.06)[(both_refs$parComb == 3) + 1]
both_refs$mu2 <- 0.03
both_refs$q01 <- 0.01
both_refs$q10 <- c(0.01, 0.005)[(both_refs$parComb == 4) + 1]
both_refs$psi <- c(0, 0.01, 0.05, 0.1)[c(rep(1, 4), fbd_refs$psiComb + 1)]

# make data frames for 95% CI and median, one for mean, and one for coverage
both_low <- both_med <- both_high <- both_mean <- both_cov <- both_pce <-
  data.frame(matrix(nrow = 0, ncol = 10))

# iterate through logs
for (i in 1:nrow(both_refs)) {
  # get ref
  ref <- both_refs$comb[i]
  
  # and reps
  for (j in 1:n_reps) {
    # get log
    log <- logs[[ref]][[j]]
    
    # check if it's extant
    if (ref %in% c(7, 11, 15, 19)) {
      # add psi columns to log
      log <- cbind(log[, 1:4], rep(0, nrow(log)), rep(0, nrow(log)), log[, 5:8])
    }
    
    # apply burnout and take out pi
    log <- log[(nrow(log)/5):nrow(log), -c(7, 8)]
    
    # get quantiles
    quants <- do.call(rbind.data.frame, lapply(c(0.025, 0.5, 0.975), 
                                               function(x)
                                                 unlist(lapply(1:ncol(log), function(y) quantile(log[, y], x)))))
    quants <- cbind(rep(ref, 3), quants, rep(i, 3))
    colnames(quants) <- c("comb", "lambda1", "lambda2", 
                          "mu1", "mu2", "psi1", "psi2", "q01", "q10", "comb_seq")
    
    # fill data frames
    both_low <- rbind(both_low, quants[1, ])
    both_med <- rbind(both_med, quants[2, ])
    both_high <- rbind(both_high, quants[3, ])
    both_mean <- rbind(both_mean, c(ref, colMeans(log), i))
  }
  
  # low and high dfs for this combination
  both_low_i <- both_low[both_low[, 1] == ref, ]
  both_high_i <- both_high[both_high[, 1] == ref, ]
  both_mean_i <- both_mean[both_mean[, 1] == ref, ]
  
  # coverage
  both_cov <- rbind(both_cov, 
                     c(sum(both_low_i[, 2] < both_refs$lambda1[i] & 
                             both_high_i[, 2] > both_refs$lambda1[i]),
                       sum(both_low_i[, 3] < both_refs$lambda2[i] &
                             both_high_i[, 3] > both_refs$lambda2[i]),
                       sum(both_low_i[, 4] < both_refs$mu1[i] &
                             both_high_i[, 4] > both_refs$mu1[i]),
                       sum(both_low_i[, 5] < both_refs$mu2[i] &
                             both_high_i[, 5] > both_refs$mu2[i]),
                       sum(both_low_i[, 6] < both_refs$psi[i] &
                             both_high_i[, 6] > both_refs$psi[i]),
                       sum(both_low_i[, 7] < both_refs$psi[i] &
                             both_high_i[, 7] > both_refs$psi[i]),
                       sum(both_low_i[, 8] < both_refs$q01[i] &
                             both_high_i[, 8] > both_refs$q01[i]),
                       sum(both_low_i[, 9] < both_refs$q10[i] &
                             both_high_i[, 9] > both_refs$q10[i])))
  colnames(both_cov) <- c("lambda1", "lambda2", "mu1", "mu2", 
                          "psi1", "psi2", "q01", "q10")
  
  # percent error
  both_pce <- rbind(both_pce,
                     abs(c(mean((both_mean_i[, 2] - both_refs$lambda1[i]) / 
                                  both_refs$lambda1[i]),
                           mean((both_mean_i[, 3] - both_refs$lambda2[i]) / 
                                  both_refs$lambda2[i]),
                           mean((both_mean_i[, 4] - both_refs$mu1[i]) / 
                                  both_refs$mu1[i]),
                           mean((both_mean_i[, 5] - both_refs$mu2[i]) / 
                                  both_refs$mu2[i]),
                           mean((both_mean_i[, 6] - both_refs$psi[i]) /
                                  both_refs$psi[i]),
                           mean((both_mean_i[, 7] - both_refs$psi[i]) /
                                  both_refs$psi[i]),
                           mean((both_mean_i[, 8] - both_refs$q01[i]) /
                                  both_refs$q01[i]),
                           mean((both_mean_i[, 9] - both_refs$q10[i]) /
                                  both_refs$q10[i]))))
  colnames(both_pce) <- colnames(both_cov)
}

# normalize and name both_cov
both_cov <- both_cov / n_reps
rownames(both_cov) <- rownames(both_pce) <- rownames(both_refs)

# name both_mean
colnames(both_mean) <- colnames(both_low)

# round PCE
both_pce <- round(both_pce, digits = 3)

# transition rate precision
precision_q01 <- unlist(lapply(unique(both_mean$comb), 
                               function(x) 
                                 mean(both_high$q01[both_high$comb == x] - 
                                        both_low$q01[both_low$comb == x])))
names(precision_q01) <- unique(both_mean$comb)
precision_q10 <- unlist(lapply(unique(both_mean$comb), 
                               function(x) 
                                 mean(both_high$q10[both_high$comb == x] - 
                                        both_low$q10[both_low$comb == x])))
names(precision_q10) <- unique(both_mean$comb)

# add true psi to the mean df
both_mean$true_psi <- unlist(lapply(1:nrow(both_mean), function(x)
                              both_refs$psi[both_refs$comb == both_mean$comb[x]]))

# mean boxplot
both_mean_bplot_lambda1 <- ggplot(both_mean, aes(factor(comb), lambda1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = both_refs, aes(x = factor(comb), y = lambda1), 
             col = "red", size = 2) +
  labs(x = "Parameter combination") +
  theme_bw()
both_mean_bplot_lambda2 <- ggplot(both_mean, aes(factor(comb), lambda2)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = both_refs, aes(x = factor(comb), y = lambda2), 
             col = "red", size = 2) +
  labs(x = "Parameter combination") +
  theme_bw()
both_mean_bplot_mu1 <- ggplot(both_mean, aes(factor(comb), mu1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = both_refs, aes(x = factor(comb), y = mu1), 
             col = "red", size = 2) + 
  labs(x = "Parameter combination") +
  theme_bw()
both_mean_bplot_mu2 <- ggplot(both_mean, aes(factor(comb), mu2)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = both_refs, aes(x = factor(comb), y = mu2), 
             col = "red", size = 2) + 
  labs(x = "Parameter combination") +
  theme_bw()
both_mean_bplot_q01 <- ggplot(both_mean, aes(factor(comb), q01)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = both_refs, aes(x = factor(comb), y = q01), 
             col = "red", size = 2) + 
  labs(x = "Parameter combination") +
  theme_bw()
both_mean_bplot_q10 <- ggplot(both_mean, aes(factor(comb), q10)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = both_refs, aes(x = factor(comb), y = q10), 
             col = "red", size = 2) + 
  labs(x = "Parameter combination") +
  theme_bw()

# psi plots
both_psi <- both_mean[both_mean$comb %in% c(23:25, 35:37, 47:49, 59:61), ]
both_psi$sce[both_psi$comb %in% c(23:25)] <- 1
both_psi$sce[both_psi$comb %in% c(35:37)] <- 2
both_psi$sce[both_psi$comb %in% c(47:49)] <- 3
both_psi$sce[both_psi$comb %in% c(59:61)] <- 4

both_psi_refs <- both_refs[both_refs$comb %in% c(23:25, 35:37, 47:49, 59:61), ]
both_psi_refs$sce[both_psi_refs$comb %in% c(23:25)] <- 1
both_psi_refs$sce[both_psi_refs$comb %in% c(35:37)] <- 2
both_psi_refs$sce[both_psi_refs$comb %in% c(47:49)] <- 3
both_psi_refs$sce[both_psi_refs$comb %in% c(59:61)] <- 4
color_values <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

both_mean_bplot_psi1 <- 
  ggplot(both_psi, aes(reorder(factor(comb), true_psi, FUN = mean), psi1,
                        fill = factor(sce))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = both_psi_refs,
             aes(x = reorder(factor(comb), psi), y = psi),
             col = "red", size = 2) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 2, 3, 4),
                    values = color_values) +
  labs(title = expression("Mean posterior estimates of "*psi['0']),
       x = "Parameter combination",
       y = expression(psi['0'])) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

# ridge plots
both_mean_lambda <- both_mean[both_mean$comb %in% c(7, 11, 23:25, 35:37), ]
both_mean_lambda$comb[both_mean_lambda$comb %in% c(7, 23:25)] <- 1
both_mean_lambda$comb[both_mean_lambda$comb %in% c(11, 35:37)] <- 2

both_ridge_lambda1 <- ggplot(both_mean_lambda, 
                             aes(x = lambda1, y = comb, 
                                 fill = factor(comb == 2),
                                 color = factor(comb == 2))) +
  geom_density_ridges() +
  geom_vline(xintercept = 0.1, color = "black", lwd = 1) +
  xlim(0, 0.21) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 2),
                    values = c("#E69F00", "#56B4E9", 
                               "#E69F00", "#56B4E9",
                               "#E69F00", "#56B4E9",
                               "#E69F00", "#56B4E9")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9", 
                                "#E69F00", "#56B4E9",
                                "#E69F00", "#56B4E9",
                                "#E69F00", "#56B4E9")) + 
  labs(title = expression("Mean "*lambda['0']*" estimates"), 
       x = expression("Estimated "*lambda['0'])) +
  facet_grid(. ~ fct_inorder(glue('psi*" = {true_psi}"')), labeller = label_parsed) +
  theme_bw() +
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14))

both_ridge_lambda2 <- ggplot(both_mean_lambda, 
                           aes(x = lambda2, y = comb, 
                               fill = factor(comb == 2),
                               color = factor(comb == 2))) +
  geom_density_ridges() +
  geom_vline(xintercept = c(0.1, 0.2), color = c("#E69F00", "#56B4E9", 
                                                 "#E69F00", "#56B4E9", 
                                                 "#E69F00", "#56B4E9"),
             lwd = c(1, 1, 1, 1, 1, 1)) +
  scale_fill_manual(name = "Scenario",
                     labels = c(1, 2),
                     values = c("#E69F00", "#56B4E9", 
                                "#E69F00", "#56B4E9",
                                "#E69F00", "#56B4E9")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9", 
                                "#E69F00", "#56B4E9",
                                "#E69F00", "#56B4E9")) + 
  labs(title = expression("Mean "*lambda['1']*" estimates"), 
       x = expression("Estimated "*lambda['1'])) +
  facet_grid(. ~ glue('psi*" = {true_psi}"'), labeller = label_parsed) +
  geom_text(data = data.frame(label = c(paste0("Cov = ", both_cov$lambda2[1]),
                                        paste0("Cov = ", both_cov$lambda2[2]),
                                        paste0("Cov = ", both_cov$lambda2[3])),
                              comb = 1,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.33, y = 1.75, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("MPE = ", both_pce$lambda2[1]),
                                        paste0("MPE = ", both_pce$lambda2[2]),
                                        paste0("MPE = ", both_pce$lambda2[3])),
                              comb = 1,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.33, y = 1.5, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("Cov = ", both_cov$lambda2[4]),
                                        paste0("Cov = ", both_cov$lambda2[5]),
                                        paste0("Cov = ", both_cov$lambda2[6])),
                              comb = 2,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.33, y = 4, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("MPE = ", both_pce$lambda2[4]),
                                        paste0("MPE = ", both_pce$lambda2[5]),
                                        paste0("MPE = ", both_pce$lambda2[6])),
                              comb = 2,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.33, y = 3.75, label = label),
            size = 4.5) +
  theme_bw() +
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14))

both_mean_mu <- both_mean[both_mean$comb %in% c(23:25, 47:49), ]
both_mean_mu$comb[both_mean_mu$comb %in% 23:25] <- 1
both_mean_mu$comb[both_mean_mu$comb %in% 47:49] <- 2

both_ridge_mu1 <- ggplot(both_mean_mu, 
                             aes(x = mu1, y = comb, 
                                 fill = factor(comb == 2),
                                 color = factor(comb == 2))) +
  geom_density_ridges() +
  geom_vline(xintercept = c(0.03, 0.06), color = c("#E69F00", "#009E73", 
                                                 "#E69F00", "#009E73", 
                                                 "#E69F00", "#009E73"),
             lwd = c(1, 1, 1, 1, 1, 1)) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 3),
                    values = c("#E69F00", "#009E73", 
                               "#E69F00", "#009E73",
                               "#E69F00", "#009E73")) + 
  scale_color_manual(values = c("#E69F00", "#009E73", 
                                "#E69F00", "#009E73",
                                "#E69F00", "#009E73")) + 
  labs(title = expression("Mean "*mu['0']*" estimates"), 
       x = expression("Estimated "*mu['0'])) +
  facet_grid(. ~ glue('psi*" = {true_psi}"'), labeller = label_parsed) +
  geom_text(data = data.frame(label = c(paste0("Cov = ", both_cov$mu1[1]),
                                        paste0("Cov = ", both_cov$mu1[2]),
                                        paste0("Cov = ", both_cov$mu1[3])),
                              comb = 1,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.15, y = 1.75, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("MPE = ", both_pce$mu1[1]),
                                        paste0("MPE = ", both_pce$mu1[2]),
                                        paste0("MPE = ", both_pce$mu1[3])),
                              comb = 1,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.15, y = 1.5, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("Cov = ", both_cov$mu1[7]),
                                        paste0("Cov = ", both_cov$mu1[8]),
                                        paste0("Cov = ", both_cov$mu1[9])),
                              comb = 2,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.15, y = 4, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("MPE = ", both_pce$mu1[7]),
                                        paste0("MPE = ", both_pce$mu1[8]),
                                        paste0("MPE = ", both_pce$mu1[9])),
                              comb = 2,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.15, y = 3.75, label = label),
            size = 4.5) +
  theme_bw() +
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14))

both_ridge_mu2 <- ggplot(both_mean_mu, 
                             aes(x = mu2, y = comb, 
                                 fill = factor(comb == 2),
                                 color = factor(comb == 2))) +
  geom_density_ridges() +
  geom_vline(xintercept = 0.03, color = "black", lwd = 1) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 3),
                    values = c("#E69F00", "#009E73", 
                               "#E69F00", "#009E73",
                               "#E69F00", "#009E73")) + 
  scale_color_manual(values = c("#E69F00", "#009E73", 
                                "#E69F00", "#009E73",
                                "#E69F00", "#009E73")) + 
  labs(title = expression("Mean "*mu['1']*" estimates"), 
       x = expression("Estimated "*mu['1'])) +
  facet_grid(. ~ glue('psi*" = {true_psi}"'), labeller = label_parsed) +
  geom_text(data = data.frame(label = c(paste0("Cov = ", both_cov$mu2[1]),
                                        paste0("Cov = ", both_cov$mu2[2]),
                                        paste0("Cov = ", both_cov$mu2[3])),
                              comb = 1,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.14, y = 1.75, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("MPE = ", both_pce$mu2[1]),
                                        paste0("MPE = ", both_pce$mu2[2]),
                                        paste0("MPE = ", both_pce$mu2[3])),
                              comb = 1,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.14, y = 1.5, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("Cov = ", both_cov$mu2[7]),
                                        paste0("Cov = ", both_cov$mu2[8]),
                                        paste0("Cov = ", both_cov$mu2[9])),
                              comb = 2,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.14, y = 4, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("MPE = ", both_pce$mu2[7]),
                                        paste0("MPE = ", both_pce$mu2[8]),
                                        paste0("MPE = ", both_pce$mu2[9])),
                              comb = 2,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.14, y = 3.75, label = label),
            size = 4.5) +
  theme_bw() +
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14))

both_mean_q <- both_mean[both_mean$comb %in% c(23:25, 59:61), ]
both_mean_q$comb[both_mean_q$comb %in% 23:25] <- 1
both_mean_q$comb[both_mean_q$comb %in% 59:61] <- 2

both_ridge_q01 <- ggplot(both_mean_q, 
                             aes(x = q01, y = comb, 
                                 fill = factor(comb == 2),
                                 color = factor(comb == 2))) +
  geom_density_ridges() +
  geom_vline(xintercept = 0.01, color = "black", lwd = 1) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 4),
                    values = c("#E69F00", "#CC79A7", 
                               "#E69F00", "#CC79A7",
                               "#E69F00", "#CC79A7")) + 
  scale_color_manual(values = c("#E69F00", "#CC79A7", 
                                "#E69F00", "#CC79A7",
                                "#E69F00", "#CC79A7")) + 
  labs(title = expression("Mean "*q['01']*" estimates"), 
       x = expression("Estimated "*q['01'])) +
  facet_grid(. ~ glue('psi*" = {true_psi}"'), labeller = label_parsed) +
  geom_text(data = data.frame(label = c(paste0("Cov = ", both_cov$q01[1]),
                                        paste0("Cov = ", both_cov$q01[2]),
                                        paste0("Cov = ", both_cov$q01[3])),
                              comb = 1,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.025, y = 1.75, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("MPE = ", both_pce$q01[1]),
                                        paste0("MPE = ", both_pce$q01[2]),
                                        paste0("MPE = ", both_pce$q01[3])),
                              comb = 1,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.025, y = 1.5, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("Cov = ", both_cov$q01[10]),
                                        paste0("Cov = ", both_cov$q01[11]),
                                        paste0("Cov = ", both_cov$q01[12])),
                              comb = 2,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.025, y = 4, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("MPE = ", both_pce$q01[10]),
                                        paste0("MPE = ", both_pce$q01[11]),
                                        paste0("MPE = ", both_pce$q01[12])),
                              comb = 2,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.025, y = 3.75, label = label),
            size = 4.5) +
  theme_bw() +
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14))

both_ridge_q10 <- ggplot(both_mean_q, 
                             aes(x = q10, y = comb, 
                                 fill = factor(comb == 2),
                                 color = factor(comb == 2))) +
  geom_density_ridges() +
  geom_vline(xintercept = c(0.01, 0.005), color = c("#E69F00", "#CC79A7", 
                                                 "#E69F00", "#CC79A7", 
                                                 "#E69F00", "#CC79A7"),
             lwd = c(1, 1, 1, 1, 1, 1)) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 4),
                    values = c("#E69F00", "#CC79A7", 
                               "#E69F00", "#CC79A7",
                               "#E69F00", "#CC79A7")) + 
  scale_color_manual(values = c("#E69F00", "#CC79A7", 
                                "#E69F00", "#CC79A7",
                                "#E69F00", "#CC79A7")) + 
  labs(title = expression("Mean "*q['10']*" estimates"), 
       x = expression("Estimated "*q['10'])) +
  facet_grid(. ~ glue('psi*" = {true_psi}"'), labeller = label_parsed) +
  geom_text(data = data.frame(label = c(paste0("Cov = ", both_cov$q10[1]),
                                        paste0("Cov = ", both_cov$q10[2]),
                                        paste0("Cov = ", both_cov$q10[3])),
                              comb = 1,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.025, y = 1.75, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("MPE = ", both_pce$q10[1]),
                                        paste0("MPE = ", both_pce$q10[2]),
                                        paste0("MPE = ", both_pce$q10[3])),
                              comb = 1,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.025, y = 1.5, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("Cov = ", both_cov$q10[10]),
                                        paste0("Cov = ", both_cov$q10[11]),
                                        paste0("Cov = ", both_cov$q10[12])),
                              comb = 2,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.025, y = 4, label = label),
            size = 4.5) +
  geom_text(data = data.frame(label = c(paste0("MPE = ", both_pce$q10[10]),
                                        paste0("MPE = ", both_pce$q10[11]),
                                        paste0("MPE = ", both_pce$q10[12])),
                              comb = 2,
                              true_psi = c(0.01, 0.05, 0.1)),
            mapping = aes(x = 0.025, y = 3.75, label = label),
            size = 4.5) +
  theme_bw() +
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14))

# plots equivalent to Maddison 2007 - means
both_comp_lambda <- ggplot(both_mean[both_mean$comb %in% c(23:25, 35:37), ], 
                           aes(lambda1, lambda2, 
                               color = factor(comb %in% 35:37))) +
  geom_point() +
  geom_hline(yintercept = c(0.1, 0.2), color = c("#E69F00", "#56B4E9", 
                                                   "#E69F00", "#56B4E9", 
                                                   "#E69F00", "#56B4E9")) +
  scale_color_manual(name = expression("True "*lambda['1']),
                     labels = c(0.1, 0.2),
                     values = c("#E69F00", "#56B4E9", "#E69F00", "#56B4E9",
                                "#E69F00", "#56B4E9")) + 
  labs(title = expression("Speciation rate estimates"), 
       x = expression("Estimated "*lambda['0']),
       y = expression("Estimated "*lambda['1'])) +
  facet_grid(. ~ glue('psi*" = {true_psi}"'), labeller = label_parsed) +
  theme_bw() +
  theme(title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14))

both_comp_mu <- ggplot(both_mean[both_mean$comb %in% c(23:25, 47:49), ], 
                       aes(mu1, mu2, 
                           color = factor(comb %in% 47:49))) +
  geom_point() +
  geom_vline(xintercept = c(0.03, 0.06), color = c("#E69F00", "#56B4E9", 
                                                   "#E69F00", "#56B4E9", 
                                                   "#E69F00", "#56B4E9")) +
  scale_color_manual(name = expression("True "*mu['0']),
                     labels = c(0.03, 0.06),
                     values = c("#E69F00", "#56B4E9", "#E69F00", "#56B4E9",
                                "#E69F00", "#56B4E9")) + 
  labs(title = expression("Extinction rate estimates"), 
       x = expression("Estimated "*mu['0']),
       y = expression("Estimated "*mu['1']), 
       grid = expression("True "*psi)) +
  facet_grid(. ~ glue('psi*" = {true_psi}"'), labeller = label_parsed) +
  theme_bw() +
  theme(title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14))


both_comp_q <- ggplot(both_mean, aes(q01, q10, 
                                     color = factor(comb %in% 59:61))) +
  geom_point() +
  geom_hline(yintercept = c(0.01, 0.005), color = c("#E69F00", "#56B4E9", 
                                                   "#E69F00", "#56B4E9", 
                                                   "#E69F00", "#56B4E9")) +
  scale_color_manual(name = expression("True "*q['10']),
                     labels = c(0.01, 0.005),
                     values = c("#E69F00", "#56B4E9", "#E69F00", "#56B4E9",
                                "#E69F00", "#56B4E9")) + 
  labs(title = expression("Transition rate estimates"), 
       x = expression("Estimated "*q['01']),
       y = expression("Estimated "*q['10']), 
       grid = expression("True "*psi)) +
  facet_grid(. ~ glue('psi*" = {true_psi}"'), labeller = label_parsed) +
  theme_bw() +
  theme(title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14))

# final complete facet_grid
both_mean_total <- both_mean[, -c(6,7)]
both_mean_total$comb[both_mean_total$comb %in% c(7, 23:25)] <- 1
both_mean_total$comb[both_mean_total$comb %in% c(11, 35:37)] <- 2
both_mean_total$comb[both_mean_total$comb %in% c(15, 47:49)] <- 3
both_mean_total$comb[both_mean_total$comb %in% c(19, 59:61)] <- 4

rates_mean <- melt(both_mean_total, id.vars = c("comb_seq", "comb", "true_psi"))

rates_mean <- rates_mean[!(rates_mean$comb == 2 & !(rates_mean$variable %in% c("lambda1", "lambda2"))), ]
rates_mean <- rates_mean[!(rates_mean$comb == 3 & !(rates_mean$variable %in% c("mu1", "mu2"))), ]
rates_mean <- rates_mean[!(rates_mean$comb == 4 & !(rates_mean$variable %in% c("q01", "q10"))), ]
rates_mean$comb_seq <- rates_mean$comb
rates_mean$comb_seq[rates_mean$comb_seq != 1] <- 2

rates_mean$variable <- factor(rates_mean$variable, 
                              levels = unique(rates_mean$variable),
                              labels = c(expression(lambda['0']), expression(lambda['1']),
                                         expression(mu['0']), expression(mu['1']),
                                         expression(q['01']), expression(q['10'])))
rates_mean$true_psi <- factor(rates_mean$true_psi,
                              levels = unique(rates_mean$true_psi),
                              labels = fct_inorder(glue('psi*" = {unique(rates_mean$true_psi)}"')))

# now for coverage
both_cov_total <- both_cov[, -c(5, 6)]
both_cov_total$comb <- as.integer(rownames(both_cov_total))
both_cov_total$comb_seq <- 1:16
both_cov_total$true_psi <- both_refs$psi

both_cov_total$comb[both_cov_total$comb %in% c(7, 23:25)] <- 1
both_cov_total$comb[both_cov_total$comb %in% c(11, 35:37)] <- 2
both_cov_total$comb[both_cov_total$comb %in% c(15, 47:49)] <- 3
both_cov_total$comb[both_cov_total$comb %in% c(19, 59:61)] <- 4

total_covs <- melt(both_cov_total, id.vars = c("comb_seq", "comb", "true_psi"))

total_covs <- total_covs[!(total_covs$comb == 2 & !(total_covs$variable %in% c("lambda1", "lambda2"))), ]
total_covs <- total_covs[!(total_covs$comb == 3 & !(total_covs$variable %in% c("mu1", "mu2"))), ]
total_covs <- total_covs[!(total_covs$comb == 4 & !(total_covs$variable %in% c("q01", "q10"))), ]
total_covs$comb_seq <- total_covs$comb
total_covs$comb_seq[total_covs$comb_seq != 1] <- 2

total_covs$variable <- factor(total_covs$variable, 
                              levels = unique(total_covs$variable),
                              labels = c(expression(lambda['0']), expression(lambda['1']),
                                         expression(mu['0']), expression(mu['1']),
                                         expression(q['01']), expression(q['10'])))
total_covs$true_psi <- factor(total_covs$true_psi,
                              levels = unique(total_covs$true_psi),
                              labels = fct_inorder(glue('psi*" = {unique(total_covs$true_psi)}"')))
total_covs$label <- paste0("CP = ", total_covs$value, " ")

true_rates <- data.frame(rate = c(rep(0.1, 8), rep(c(0.1, 0.2), 4),
                                  rep(c(0.03, 0.06), 4), rep(0.03, 8),
                                  rep(0.01, 8), rep(c(0.01, 0.005), 4)),
                         variable = c(rep("lambda1", 8), rep("lambda2", 8),
                                      rep("mu1", 8), rep("mu2", 8),
                                      rep("q01", 8), rep("q10", 8)),
                         comb = c(rep(c(1, 2), 8), rep(c(1, 3), 8),
                                  rep(c(1, 4), 8)),
                         true_psi = factor(rep(c(0, 0, 0.01, 0.01, 0.05, 0.05, 0.1, 0.1), 6)))
true_rates$true_psi <- factor(true_rates$true_psi,
                              levels = unique(true_rates$true_psi),
                              labels = fct_inorder(glue('psi*" = {unique(true_rates$true_psi)}"')))
true_rates$variable <- factor(true_rates$variable, 
                              levels = unique(true_rates$variable),
                              labels = c(expression(lambda['0']), expression(lambda['1']),
                                         expression(mu['0']), expression(mu['1']),
                                         expression(q['01']), expression(q['10'])))

color_values <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
  
rates_ridge <- ggplot(rates_mean, 
                      aes(x = value, 
                          y = factor(comb_seq), 
                          fill = factor(comb),
                          color = factor(comb))) +
  geom_density_ridges() +
  geom_vline(data = true_rates, aes(xintercept = rate, color = factor(comb)), lwd = 1) +
  scale_x_continuous(n.breaks = 4, limits = c(0, NA)) +
  scale_y_discrete(expand = expansion(mult = c(0.3, 2))) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 2, 3, 4),
                    values = color_values) +
  scale_color_manual(values = color_values) +
  #labs(title = expression("Mean rate estimates")) +
  facet_grid(true_psi ~ variable, 
             labeller = label_parsed, 
             scales = "free") +
  theme_bw() +
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

# make a pce value for rates_mean
rates_mean$pce <- rep(0, nrow(rates_mean))
total_covs$max <- rep(0, nrow(total_covs))
for (i in 1:nrow(rates_mean)) {
  comb <- rates_mean$comb[i]
  tpsi <- rates_mean$true_psi[i]
  var <- rates_mean$variable[i]
  
  true_rate <- true_rates$rate[true_rates$comb == comb &
                                 true_rates$true_psi == tpsi &
                                 true_rates$variable == var]
  
  rates_mean$pce[i] <- (rates_mean$value[i] - true_rate) / true_rate
  
  if (i %% 100 == 0) {
    total_covs$max[i / 100] <- max(rates_mean$pce[(i - 100 + 1):i])
  }
}

pce_ridge <- ggplot(rates_mean, 
                    aes(x = pce, 
                        y = factor(comb_seq), 
                        #fill = factor(comb),
                        fill = factor(comb))) +
  geom_boxplot() +
  #geom_vline(data = true_rates, aes(xintercept = rate, color = factor(comb)), lwd = 1) +
  #scale_x_continuous(n.breaks = 4, limits = c(0, NA)) +
  scale_y_discrete(expand = expansion(mult = c(0.2, 0.4))) +
  scale_fill_manual(name = "Scenario",
                     labels = c(1, 2, 3, 4),
                     values = color_values) +
  geom_text(data = total_covs,
            mapping = aes(x = Inf, y = factor(comb_seq - 0.5), label = label,
                          hjust = 1, color = factor(comb)),
            size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = color_values) +
  #labs(title = expression("Percent error of mean estimates")) +
  facet_grid(true_psi ~ variable, 
             labeller = label_parsed, 
             scales = "free") +
  theme_bw() +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))


cov_bar <- ggplot(total_covs, 
                  aes(x = value, 
                      y = factor(comb_seq), 
                      fill = factor(comb))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 2, 3, 4),
                    values = color_values) +
  labs(title = expression("Coverage of posterior distributions")) +
  facet_grid(true_psi ~ variable, 
             labeller = label_parsed, 
             scales = "free") +
  theme_bw() +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

# psi plots
both_psi <- both_mean[both_mean$comb %in% c(23:25, 35:37, 47:49, 59:61), ]
both_psi$sce[both_psi$comb %in% c(23:25)] <- 1
both_psi$sce[both_psi$comb %in% c(35:37)] <- 2
both_psi$sce[both_psi$comb %in% c(47:49)] <- 3
both_psi$sce[both_psi$comb %in% c(59:61)] <- 4

both_psi_refs <- both_refs[both_refs$comb %in% c(23:25, 35:37, 47:49, 59:61), ]
both_psi_refs$sce[both_psi_refs$comb %in% c(23:25)] <- 1
both_psi_refs$sce[both_psi_refs$comb %in% c(35:37)] <- 2
both_psi_refs$sce[both_psi_refs$comb %in% c(47:49)] <- 3
both_psi_refs$sce[both_psi_refs$comb %in% c(59:61)] <- 4
color_values <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

both_mean_bplot_psi1 <- 
  ggplot(both_psi, aes(reorder(factor(comb), true_psi, FUN = mean), psi1,
                       fill = factor(sce))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = both_psi_refs,
             aes(x = reorder(factor(comb), psi), y = psi),
             col = "#E34234", size = 2,
             show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 0.15)) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 2, 3, 4),
                    values = color_values) +
  guides(color = "none") +
  labs(x = "Parameter combination",
       y = expression(psi['0'])) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

both_mean_bplot_psi2 <- 
  ggplot(both_psi, aes(reorder(factor(comb), true_psi, FUN = mean), psi2,
                       fill = factor(sce))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = both_psi_refs,
             aes(x = reorder(factor(comb), psi), y = psi),
             col = "#E34234", size = 2,
             show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 0.15)) + 
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 2, 3, 4),
                    values = color_values) +
  guides(color = "none") +
  labs(x = "Parameter combination",
       y = expression(psi['1'])) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

###
## testing false positive - BiSSE 

# get log indices for the BiSSE reps
fp_bisse_refs <- refs_df[which(refs_df[, 3] == 2 & refs_df[, 2] != 4), ]
fp_bisse_refs$comb <- as.numeric(rownames(fp_bisse_refs))

# add true values to the data frame
fp_bisse_refs$lambda1 <- 0.1
fp_bisse_refs$lambda2 <- c(0.1, 0.2)[(fp_bisse_refs$parComb == 2) + 1]
fp_bisse_refs$mu1 <- c(0.03, 0.06)[(fp_bisse_refs$parComb == 3) + 1]
fp_bisse_refs$mu2 <- 0.03

# make data frames for 95% CI and median, one for mean, and one for coverage
fp_bisse_mode <- data.frame(matrix(nrow = 0, ncol = 4))

# iterate through logs
for (i in 1:nrow(fp_bisse_refs)) {
  # get ref
  ref <- fp_bisse_refs$comb[i]
  
  # and reps
  for (j in 1:n_reps) {
    # get log
    log <- logs[[ref]][[j]]
    
    # apply burnout and take out pi
    log <- log[(nrow(log)/5):nrow(log), 1:2]
    
    # modes
    modes <- unlist(lapply(1:ncol(log), function(x) 
      density(log[, x])$x[which.max(density(log[, x])$y)]))
    
    # pvalue
    post_prob <- mean(log[, 2] > log[, 1])
      #ifelse(mean(log[, 2] > log[, 1]) == 1, 151, 
          #      mean(log[, 2] > log[, 1]) / (1 - mean(log[, 2] > log[, 1])))
    
    # add modes to data frame
    fp_bisse_mode <- rbind(fp_bisse_mode, c(ref, modes, post_prob))
  }
  
  colnames(fp_bisse_mode) <- c("comb", "lambda1", "lambda2", "post_prob")
}

# boxplot of modes
fp_bisse_mode_lambda_bplot <- ggplot(fp_bisse_mode, aes(factor(comb), lambda2 - lambda1)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x = "Parameter combination") +
  theme_bw()

# change df a bit
fp_bisse_mode <- fp_bisse_mode[fp_bisse_mode$comb %in% 7:14, ]
fp_bisse_mode$parComb <- unlist(lapply(1:nrow(fp_bisse_mode), 
        function(x) 
          fp_bisse_refs$parComb[fp_bisse_refs$comb == fp_bisse_mode$comb[x]]))
fp_bisse_mode$parComb <- factor(fp_bisse_mode$parComb, levels = c(2, 1))
fp_bisse_mode$traitComb <- unlist(lapply(1:nrow(fp_bisse_mode), 
        function(x) 
          fp_bisse_refs$traitComb[fp_bisse_refs$comb == fp_bisse_mode$comb[x]]))

# make a labeller for facet_grid
trait_labs <- c("Effect trait", "q = 0.01", "q = 0.1", "q = 1")
names(trait_labs) <- 1:4
par_labs <- c("q = 0.01", "No shifts")
names(par_labs) <- 2:1

# plot facet_grid
fp_bisse_hist <- ggplot(fp_bisse_mode, aes(post_prob)) +
  geom_histogram(aes(y = stat(count / sum(count)),
                     fill = parComb, color = parComb), binwidth = 0.1) +
  facet_grid(parComb ~ traitComb, 
             labeller = labeller(parComb = par_labs, traitComb = trait_labs)) +
  scale_fill_manual(name = "Scenario", values = c("#56B4E9", "#E69F00")) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  labs(#title = expression("Posterior probability of "*lambda['1']*" > "*lambda['0']),
       x = "Posterior probability", 
       y = "Proportion of simulations") +
  guides(color = "none") +
  theme_bw() +
  theme(title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

###
## testing false positive - BiSSE + FBD

# get log indices for the BiSSE reps
fp_both_refs <- refs_df[which(refs_df[, 3] == 3 & refs_df[, 2] != 4), ]
fp_both_refs$comb <- as.numeric(rownames(fp_both_refs))

# add true values to the data frame
fp_both_refs$lambda1 <- 0.1
fp_both_refs$lambda2 <- c(0.1, 0.2)[(fp_both_refs$parComb == 2) + 1]
fp_both_refs$mu1 <- c(0.03, 0.06)[(fp_both_refs$parComb == 3) + 1]
fp_both_refs$mu2 <- 0.03
fp_both_refs$psi <- c(0.01, 0.05, 0.1)[fbd_refs$psiComb]

# make data frames for 95% CI and median, one for mean, and one for coverage
fp_both_mode <- data.frame(matrix(nrow = 0, ncol = 9))

# iterate through logs
for (i in 1:nrow(fp_both_refs)) {
  # get ref
  ref <- fp_both_refs$comb[i]
  
  # and reps
  for (j in 1:n_reps) {
    # get log
    log <- logs[[ref]][[j]]
    
    # apply burnout and take out pi
    log <- log[(nrow(log)/5):nrow(log), 1:6]
    
    # modes
    modes <- unlist(lapply(1:ncol(log), function(x) 
      density(log[, x])$x[which.max(density(log[, x])$y)]))
    
    # Bayes factor
    lambda_post_prob <- mean(log[, 2] > log[, 1]) #/ (1 - mean(log[, 2] > log[, 1]))
    mu_post_prob <- mean(log[, 3] > log[, 4]) #/ (1 - mean(log[, 4] > log[, 3]))
    
    # add modes to data frame
    fp_both_mode <- rbind(fp_both_mode, c(ref, modes, 
                                          lambda_post_prob, mu_post_prob))
  }
  
  colnames(fp_both_mode) <- c("comb", "lambda1", "lambda2", "mu1", "mu2",
                              "psi1", "psi2",
                              "lambda_post_prob", "mu_post_prob")
}

# boxplot of modes
fp_both_mode_lambda_bplot <- ggplot(fp_both_mode, aes(factor(comb), lambda2 - lambda1)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x = "Parameter combination") +
  theme_bw()
fp_both_mode_mu_bplot <- ggplot(fp_both_mode, aes(factor(comb), mu2 - mu1)) + 
  geom_boxplot(outlier.shape = NA) + 
  labs(x = "Parameter combination") +
  theme_bw()

# boxplot of sampling to make sure it's not wrong
fp_both_mode_bplot_psi1 <- 
  ggplot(fp_both_mode, aes(reorder(factor(comb), psi1, FUN = mean), psi1)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = fp_both_refs,
             aes(x = reorder(factor(comb), psi), y = psi),
             col = "red", size = 2) + 
  labs(x = "Parameter combination") + 
  theme_bw()
fp_both_mode_bplot_psi2 <- 
  ggplot(fp_both_mode, aes(reorder(factor(comb), psi2, FUN = mean), psi2)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data = fp_both_refs,
             aes(x = reorder(factor(comb), psi), y = psi),
             col = "red", size = 2) + 
  labs(x = "Parameter combination") + 
  theme_bw()

## lambda first
# change df a bit
fp_both_mode_lambda <- fp_both_mode[fp_both_mode$comb %in% 23:46, ]
fp_both_mode_lambda$parComb <- unlist(lapply(1:nrow(fp_both_mode_lambda), 
      function(x) 
        fp_both_refs$parComb[fp_both_refs$comb == fp_both_mode_lambda$comb[x]]))
fp_both_mode_lambda$parComb <- factor(fp_both_mode_lambda$parComb, 
                                      levels = c(2, 1))
fp_both_mode_lambda$traitComb <- unlist(lapply(1:nrow(fp_both_mode_lambda), 
      function(x) 
        fp_both_refs$traitComb[fp_both_refs$comb == fp_both_mode_lambda$comb[x]]))
fp_both_mode_lambda$psiComb <- unlist(lapply(1:nrow(fp_both_mode_lambda),
      function(x)
        fp_both_refs$psiComb[fp_both_refs$comb == fp_both_mode_lambda$comb[x]]))

# make a labeller for facet_grid
trait_labs <- c("Effect trait", "q = 0.01", "q = 0.1", "q = 1")
names(trait_labs) <- 1:4
par_labs <- c("q = 0.01", "No shifts")
names(par_labs) <- 2:1

# plot facet_grid
fp_both_lambda_low <- 
  ggplot(fp_both_mode_lambda[fp_both_mode_lambda$psiComb == 1, ], 
      aes(lambda_post_prob)) +
geom_histogram(aes(y = stat(count / sum(count)), 
                   fill = parComb, color = parComb), binwidth = 0.1) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
facet_grid(parComb ~ traitComb, 
           labeller = labeller(parComb = par_labs, 
                               traitComb = trait_labs)) + 
labs(#title = expression("Posterior probability of "*lambda['1']*" > "*lambda['0']*", "*psi*" = "*0.01),
     x = "Posterior probability", 
     y = "Proportion of simulations") +
theme_bw() +
theme(legend.position = "none",
      title = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      strip.text.x = element_text(size = 14),
      strip.text.y = element_text(size = 14))

fp_both_lambda_mid <- 
  ggplot(fp_both_mode_lambda[fp_both_mode_lambda$psiComb == 2, ], 
         aes(lambda_post_prob)) +
  geom_histogram(aes(y = stat(count / sum(count)), 
                     fill = parComb, color = parComb), binwidth = 0.1) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  facet_grid(parComb ~ traitComb, 
             labeller = labeller(parComb = par_labs, 
                                 traitComb = trait_labs)) + 
  labs(#title = expression("Posterior probability of "*lambda['1']*" > "*lambda['0']*", "*psi*" = "*0.05),
       x = "Posterior probability", 
       y = "Proportion of simulations") +
  theme_bw() +
  theme(legend.position = "none",
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

fp_both_lambda_high <- 
  ggplot(fp_both_mode_lambda[fp_both_mode_lambda$psiComb == 3, ], 
         aes(lambda_post_prob)) +
  geom_histogram(aes(y = stat(count / sum(count)), 
                     fill = parComb, color = parComb), binwidth = 0.1) +
  scale_fill_manual(name = "Scenario",
                    values = c("#56B4E9", "#E69F00")) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  facet_grid(parComb ~ traitComb, 
             labeller = labeller(parComb = par_labs, 
                                 traitComb = trait_labs)) + 
  labs(#title = expression("Posterior probability of "*lambda['1']*" > "*lambda['0']*", "*psi*" = "*0.1),
       x = "Posterior probability", 
       y = "Proportion of simulations") +
  theme_bw() +
  guides(color = "none") +
  theme(title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

## mu
# change df a bit
fp_both_mode_mu <- fp_both_mode[fp_both_mode$comb %in% 23:46, ]
fp_both_mode_mu$parComb <- unlist(lapply(1:nrow(fp_both_mode_mu), 
                                             function(x) 
                                               fp_both_refs$parComb[fp_both_refs$comb == fp_both_mode_mu$comb[x]]))
fp_both_mode_mu$parComb <- factor(fp_both_mode_mu$parComb, 
                                      levels = c(2, 1))
fp_both_mode_mu$traitComb <- unlist(lapply(1:nrow(fp_both_mode_mu), 
                                               function(x) 
                                                 fp_both_refs$traitComb[fp_both_refs$comb == fp_both_mode_mu$comb[x]]))
fp_both_mode_mu$psiComb <- unlist(lapply(1:nrow(fp_both_mode_mu),
                                             function(x)
                                               fp_both_refs$psiComb[fp_both_refs$comb == fp_both_mode_mu$comb[x]]))

# make a labeller for facet_grid
trait_labs <- c("Effect trait", "q = 0.01", "q = 0.1", "q = 1")
names(trait_labs) <- 1:4
par_labs <- c("q = 0.01", "No shifts")
names(par_labs) <- 2:1

# plot facet_grid
fp_both_mu_low <- 
  ggplot(fp_both_mode_mu[fp_both_mode_mu$psiComb == 1, ], 
         aes(mu_post_prob)) +
  geom_histogram(aes(y = stat(count / sum(count)), 
                     fill = parComb, color = parComb), binwidth = 0.1) +
  scale_fill_manual(values = c("#009E73", "#56B4E9")) +
  scale_color_manual(values = c("black", "#0072B2")) +
  facet_grid(parComb ~ traitComb, 
             labeller = labeller(parComb = par_labs, 
                                 traitComb = trait_labs)) + 
  labs(#title = expression("Posterior probability of "*mu['0']*" > "*mu['1']*", "*psi*" = "*0.01),
       x = "Posterior probability", 
       y = "Proportion of simulations") +
  theme_bw() +
  theme(legend.position = "none",
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

fp_both_mu_mid <- 
  ggplot(fp_both_mode_mu[fp_both_mode_mu$psiComb == 2, ], 
         aes(mu_post_prob)) +
  geom_histogram(aes(y = stat(count / sum(count)), 
                     fill = parComb, color = parComb), binwidth = 0.1) +
  scale_fill_manual(values = c("#009E73", "#56B4E9")) +
  scale_color_manual(values = c("black", "#0072B2")) +
  facet_grid(parComb ~ traitComb, 
             labeller = labeller(parComb = par_labs, 
                                 traitComb = trait_labs)) + 
  labs(#title = expression("Posterior probability of "*mu['0']*" > "*mu['1']*", "*psi*" = "*0.05),
       x = "Posterior probability", 
       y = "Proportion of simulations") +
  theme_bw() +
  theme(legend.position = "none",
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

fp_both_mu_high <- 
  ggplot(fp_both_mode_mu[fp_both_mode_mu$psiComb == 3, ], 
         aes(mu_post_prob)) +
  geom_histogram(aes(y = stat(count / sum(count)), 
                     fill = parComb, color = parComb), binwidth = 0.1) +
  scale_fill_manual(values = c("#009E73", "#56B4E9")) +
  scale_color_manual(values = c("black", "#0072B2")) +
  facet_grid(parComb ~ traitComb, 
             labeller = labeller(parComb = par_labs, 
                                 traitComb = trait_labs)) + 
  labs(#title = expression("Posterior probability of "*mu['0']*" > "*mu['1']*", "*psi*" = "*0.1),
       x = "Posterior probability", 
       y = "Proportion of simulations") +
  theme_bw() +
  theme(legend.position = "none",
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))
