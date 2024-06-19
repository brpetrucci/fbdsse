#####################################
##            Chapter 1            ##
##   BiSSE + FBD Simulation Study  ##
##            Plotting             ##
##    Bruno do Rosario Petrucci    ##
#####################################

### 
## packages

# ggplot2 and other graphics packages
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
## testing accuracy - BiSSE + FBD

# get log indices for modelComb 2 or 3 (BiSSE and BiSSE+FBD)
# and traitComb 1 (focus trait, not neutral traits)
both_refs <- refs_df[which((refs_df[, 3] == 3 | refs_df[, 3] == 2) & 
                             refs_df[, 4] == 1), ]
both_refs$comb <- as.numeric(rownames(both_refs))
both_refs$comb_seq <- 1:nrow(both_refs)

# add true values to the data frame
both_refs$lambda1 <- 0.1
both_refs$lambda2 <- c(0.1, 0.2)[(both_refs$parComb == 2) + 1]
both_refs$mu1 <- c(0.03, 0.06)[(both_refs$parComb == 3) + 1]
both_refs$mu2 <- 0.03
both_refs$q01 <- 0.01
both_refs$q10 <- c(0.01, 0.005)[(both_refs$parComb == 4) + 1]
both_refs$psi <- c(0, 0.01, 0.05, 0.1)[c(rep(1, 4), both_refs$psiComb[5:16] + 1)]

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
    quants <- do.call(rbind.data.frame, 
                      lapply(c(0.025, 0.5, 0.975), 
                             function(x)
                               unlist(lapply(1:ncol(log), 
                                             function(y) quantile(log[, y], x)))))
    quants <- cbind(rep(ref, 3), quants, rep(i, 3))
    colnames(quants) <- c("comb", "lambda1", "lambda2", 
                          "mu1", "mu2", "psi1", "psi2", 
                          "q01", "q10", "comb_seq")
    
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

# add true psi to the mean df
both_mean$true_psi <- unlist(lapply(1:nrow(both_mean), function(x)
  both_refs$psi[both_refs$comb == both_mean$comb[x]]))

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

# final complete facet_grid
both_mean_total <- both_mean[, -c(6,7)]
both_mean_total$comb[both_mean_total$comb %in% c(7, 23:25)] <- 1
both_mean_total$comb[both_mean_total$comb %in% c(11, 35:37)] <- 2
both_mean_total$comb[both_mean_total$comb %in% c(15, 47:49)] <- 3
both_mean_total$comb[both_mean_total$comb %in% c(19, 59:61)] <- 4

rates_mods <- both_mean_total
rates_mods$tau1 <- rates_mods$lambda1 + rates_mods$mu1
rates_mods$tau2 <- rates_mods$lambda2 + rates_mods$mu2
rates_mods$eps1 <- rates_mods$mu1 / rates_mods$lambda1
rates_mods$eps2 <- rates_mods$mu2 / rates_mods$lambda2
rates_mods_total <- rates_mods[, -c(2:8)]

rates_mean <- melt(both_mean_total, id.vars = c("comb_seq", "comb", "true_psi"))
mod_rates_mean <- melt(rates_mods_total, id.vars = c("comb", "true_psi"))

rates_mean <- rates_mean[!(rates_mean$comb == 2 & !(rates_mean$variable %in% c("lambda1", "lambda2"))), ]
rates_mean <- rates_mean[!(rates_mean$comb == 3 & !(rates_mean$variable %in% c("mu1", "mu2"))), ]
rates_mean <- rates_mean[!(rates_mean$comb == 4 & !(rates_mean$variable %in% c("q01", "q10"))), ]
rates_mean$comb_seq <- rates_mean$comb
rates_mean$comb_seq[rates_mean$comb_seq != 1] <- 2

mod_rates_mean <- mod_rates_mean[!(mod_rates_mean$comb == 4), ]

rates_mean$variable <- factor(rates_mean$variable, 
                              levels = unique(rates_mean$variable),
                              labels = c(expression(lambda['0']), expression(lambda['1']),
                                         expression(mu['0']), expression(mu['1']),
                                         expression(q['01']), expression(q['10'])))
rates_mean$true_psi <- factor(rates_mean$true_psi,
                              levels = unique(rates_mean$true_psi),
                              labels = fct_inorder(glue('psi*" = {unique(rates_mean$true_psi)}"')))

mod_rates_mean$variable <- factor(mod_rates_mean$variable, 
                                  levels = unique(mod_rates_mean$variable),
                                  labels = c(expression(tau['0']), expression(tau['1']),
                                             expression(epsilon['0']), expression(epsilon['1'])))
mod_rates_mean$true_psi <- factor(mod_rates_mean$true_psi,
                                  levels = unique(mod_rates_mean$true_psi),
                                  labels = fct_inorder(glue('psi*" = {unique(mod_rates_mean$true_psi)}"')))

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

mod_true_rates <- data.frame(comb = c(1, 2, 3),
                             tau1 = c(0.1 + 0.03, 0.1 + 0.03, 0.1 + 0.06),
                             tau2 = c(0.1 + 0.03, 0.2 + 0.03, 0.1 + 0.03),
                             eps1 = c(0.03/0.1, 0.03/0.1, 0.06/0.1),
                             eps2 = c(0.03/0.1, 0.03/0.2, 0.03/0.1))
mod_true_rates <- do.call("rbind", replicate(4, mod_true_rates, simplify = FALSE))
mod_true_rates$true_psi <- factor(rep(c(0, 0.01, 0.05, 0.1), 3))
mod_true_rates <- melt(mod_true_rates, id.vars = c("comb", "true_psi"))

true_rates$true_psi <- factor(true_rates$true_psi,
                              levels = unique(true_rates$true_psi),
                              labels = fct_inorder(glue('psi*" = {unique(true_rates$true_psi)}"')))
true_rates$variable <- factor(true_rates$variable, 
                              levels = unique(true_rates$variable),
                              labels = c(expression(lambda['0']), expression(lambda['1']),
                                         expression(mu['0']), expression(mu['1']),
                                         expression(q['01']), expression(q['10'])))

mod_true_rates$true_psi <- factor(mod_true_rates$true_psi,
                                  levels = unique(mod_true_rates$true_psi),
                                  labels = fct_inorder(glue('psi*" = {unique(mod_true_rates$true_psi)}"')))
mod_true_rates$variable <- factor(mod_true_rates$variable, 
                                  levels = unique(mod_true_rates$variable),
                                  labels = c(expression(tau['0']), expression(tau['1']),
                                             expression(epsilon['0']), expression(epsilon['1'])))


color_values <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

rates_ridge <- ggplot(rates_mean, 
                      aes(x = value, 
                          y = factor(comb_seq), 
                          fill = factor(comb),
                          color = factor(comb))) +
  geom_density_ridges(alpha = 0.75) +
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
dev.copy2pdf(file = paste0("/Users/petrucci/Documents/research/ssetests_chap1/", 
                           "low_sample_size/paper/figures/rates_ridge_opaque.pdf"),
             out.type = "cairo", width = 12.5, height = 7.84)

mod_rates_ridge <- ggplot(mod_rates_mean, 
                          aes(x = value, 
                              y = factor(comb), 
                              fill = factor(comb),
                              color = factor(comb))) +
  geom_density_ridges(alpha = 0.75) +
  geom_vline(data = mod_true_rates, aes(xintercept = value, color = factor(comb)), lwd = 1,
             linetype = rep(c("solid", "dashed", "dotdash"), 16)) +
  scale_x_continuous(n.breaks = 4, limits = c(0, NA)) +
  scale_y_discrete(expand = expansion(mult = c(0.3, 1))) +
  scale_fill_manual(name = "Scenario",
                    labels = c(1, 2, 3),
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
## testing false positive - BiSSE + FBD

# get log indices for the BiSSE reps
fp_both_refs <- refs_df[which(refs_df[, 3] == 3 & refs_df[, 2] != 4), ]
fp_both_refs$comb <- as.numeric(rownames(fp_both_refs))

# add true values to the data frame
fp_both_refs$lambda1 <- 0.1
fp_both_refs$lambda2 <- c(0.1, 0.2)[(fp_both_refs$parComb == 2) + 1]
fp_both_refs$mu1 <- c(0.03, 0.06)[(fp_both_refs$parComb == 3) + 1]
fp_both_refs$mu2 <- 0.03
fp_both_refs$psi <- c(0.01, 0.05, 0.1)[both_refs$psiComb[5:16]]

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