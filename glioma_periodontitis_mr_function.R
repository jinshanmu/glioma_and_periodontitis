library(ggplot2)
library(MendelianRandomization)
library(MRPRESSO)
library(patchwork)
library(RMediation)
library(TwoSampleMR)

# glioma data
glioma_melin_2017_new <- data.table::fread("glioma_melin_2017_new.tsv.gz")
glioma_melin_2017_new$phenotype <- "glioma"
gli_exp <- format_data(dat = glioma_melin_2017_new[glioma_melin_2017_new$P < 5e-05, ], 
                       type = "exposure", 
                       phenotype_col = "phenotype", 
                       snp_col = "SNP", 
                       beta_col = "BETA", 
                       se_col = "SE", 
                       eaf_col = "FRQ", 
                       effect_allele_col = "A2", 
                       other_allele_col = "A1", 
                       pval_col = "P", 
                       chr_col = "CHR", 
                       pos_col = "BP")

# GBM data
gbm_melin_2017_new <- data.table::fread("gbm_melin_2017_new.tsv.gz")
gbm_melin_2017_new$phenotype <- "GBM"
gbm_exp <- format_data(dat = gbm_melin_2017_new[gbm_melin_2017_new$P < 5e-05, ], 
                       type = "exposure", 
                       phenotype_col = "phenotype", 
                       snp_col = "SNP", 
                       beta_col = "BETA", 
                       se_col = "SE", 
                       eaf_col = "FRQ", 
                       effect_allele_col = "A2", 
                       other_allele_col = "A1", 
                       pval_col = "P", 
                       chr_col = "CHR", 
                       pos_col = "BP")

# non-GBM glioma data
non_gbm_melin_2017_new <- data.table::fread("non_gbm_melin_2017_new.tsv.gz")
non_gbm_melin_2017_new$phenotype <- "non-GBM glioma"
non_gbm_exp <- format_data(dat = non_gbm_melin_2017_new[non_gbm_melin_2017_new$P < 5e-05, ], 
                           type = "exposure", 
                           phenotype_col = "phenotype", 
                           snp_col = "SNP", 
                           beta_col = "BETA", 
                           se_col = "SE", 
                           eaf_col = "FRQ", 
                           effect_allele_col = "A2", 
                           other_allele_col = "A1", 
                           pval_col = "P", 
                           chr_col = "CHR", 
                           pos_col = "BP")

# periodontitis data
EUR_perio_excl_HCHSSOL_new <- data.table::fread("EUR_perio_excl_HCHSSOL_new.tsv.gz")
EUR_perio_excl_HCHSSOL_new$phenotype <- "periodontitis"

# primary MR analysis
mr_pri_ana <- function(exp, out, exp_sample_size) {
  # exposure
  repeat {
    tryCatch(
      {
        withCallingHandlers(
          {
            exp <- clump_data(dat = exp, 
                              clump_kb = 1000, 
                              clump_r2 = 0.1)
            break
          }, 
          message = function(m){if(grepl("Server error", m$message)){stop()}}
        )
      }, 
      error = function(e){message("Retry")}
    )
  }
  exp$samplesize.exposure <- exp_sample_size
  exp$R2 <- exp$beta.exposure ^ 2 / (exp$beta.exposure ^ 2 + exp$samplesize.exposure * exp$se.exposure ^ 2)
  exp$F <- (exp$samplesize.exposure - 2) * exp$R2 / (1 - exp$R2)
  exp <- exp[exp$F >= 10, ]
  
  # outcome
  out <- format_data(dat = out, 
                     type = "outcome", 
                     snps = exp$SNP, 
                     phenotype_col = "phenotype", 
                     snp_col = "SNP", 
                     beta_col = "BETA", 
                     se_col = "SE", 
                     eaf_col = "FRQ", 
                     effect_allele_col = "A2", 
                     other_allele_col = "A1", 
                     pval_col = "P", 
                     samplesize_col = "N", 
                     chr_col = "CHR", 
                     pos_col = "BP")
  
  # harmonization
  exp_out <- harmonise_data(exposure_dat = exp, 
                            outcome_dat = out, 
                            action = 2)
  exp_out <- exp_out[exp_out$pval.outcome >= 5e-05, ]
  
  # outlier removal
  exp_out_cd <- cooks.distance(lm(formula = beta.outcome ~ beta.exposure + 0, data = exp_out[exp_out$mr_keep == TRUE, ], weights = se.outcome ^ (-2)))
  exp_out_cd_large <- as.numeric(names(which(exp_out_cd > 4 / nrow(exp_out[exp_out$mr_keep == TRUE, ]))))
  exp_out_sr <- rstudent(lm(formula = beta.outcome ~ beta.exposure + 0, data = exp_out[exp_out$mr_keep == TRUE, ], weights = se.outcome ^ (-2)))
  exp_out_sr_large <- as.numeric(names(which(abs(exp_out_sr) > 2)))
  exp_out <- exp_out[-c(exp_out_cd_large, exp_out_sr_large), ]
  rownames(exp_out) <- 1 : nrow(exp_out)
  
  # MR
  method_list <- c("mr_egger_regression", 
                   "mr_simple_median", 
                   "mr_weighted_median", 
                   "mr_ivw", 
                   "mr_simple_mode", 
                   "mr_weighted_mode")
  mr_res_exp_out <- mr(exp_out, method_list = method_list)
  mr_res_or_exp_out <- generate_odds_ratios(mr_res_exp_out)
  
  return(list(exp_out = exp_out, 
              mr_res_or_exp_out = mr_res_or_exp_out))
}

# supplementary MR analysis
mr_sup_ana <- function(exp_out, exp_name, out_name, mr_res_exp_out) {
  # MR-PRESSO
  mr_presso_res_exp_out <- mr_presso(BetaOutcome = "beta.outcome", 
                                     BetaExposure = "beta.exposure", 
                                     SdOutcome = "se.outcome", 
                                     SdExposure = "se.exposure", 
                                     OUTLIERtest = TRUE, 
                                     DISTORTIONtest = TRUE, 
                                     data = exp_out, 
                                     NbDistribution = 1000, 
                                     SignifThreshold = 0.05, 
                                     seed = 123)
  
  # MR-LASSO
  mr_input_exp_out <- mr_input(bx = exp_out$beta.exposure, 
                               bxse = exp_out$se.exposure, 
                               by = exp_out$beta.outcome, 
                               byse = exp_out$se.outcome, 
                               exposure = exp_name, 
                               outcome = out_name, 
                               snps = exp_out$SNP)
  mr_lasso_res_exp_out <- mr_lasso(mr_input_exp_out)
  
  # pleiotropy test
  pleiotropy_exp_out <- mr_pleiotropy_test(exp_out)
  
  # heterogeneity test
  heterogeneity_exp_out <- mr_heterogeneity(exp_out)
  
  # forest plot
  method_list <- c("mr_egger_regression", 
                   "mr_simple_median", 
                   "mr_weighted_median", 
                   "mr_ivw", 
                   "mr_simple_mode", 
                   "mr_weighted_mode")
  single_snp_exp_out <- mr_singlesnp(exp_out, all_method = method_list)
  forest_plot_mr_exp_out <- mr_forest_plot(single_snp_exp_out)
  forest_plot_mr_exp_out <- forest_plot_mr_exp_out[[1]] + 
    theme_classic() + 
    theme(legend.position = "none") + 
    xlab(paste("MR effect size for", exp_name, "on", out_name))
  
  # scatter plot
  scatter_plot_mr_exp_out <- mr_scatter_plot(mr_res_exp_out, exp_out)
  scatter_plot_mr_exp_out <- scatter_plot_mr_exp_out[[1]] + 
    theme_classic() + 
    theme(legend.direction = "vertical", 
          legend.position = "right", 
          legend.background = element_blank()) + 
    xlab(paste("SNP effect on", exp_name)) + 
    ylab(paste("SNP effect on", out_name))
  scatter_plot_mr_exp_out <- scatter_plot_mr_exp_out + guides(col = guide_legend(ncol = 1, title = "MR method"))
  
  # funnel plot
  funnel_plot_mr_exp_out <- mr_funnel_plot(single_snp_exp_out)
  funnel_plot_mr_exp_out <- funnel_plot_mr_exp_out[[1]] + 
    theme_classic() + 
    theme(legend.direction = "vertical", 
          legend.position = "right", 
          legend.background = element_blank())
  funnel_plot_mr_exp_out <- funnel_plot_mr_exp_out + guides(col = guide_legend(ncol = 1, title = "MR method"))
  
  # density plot
  density_plot_mr_exp_out <- mr_density_plot(single_snp_exp_out, mr_res_exp_out)
  density_plot_mr_exp_out <- density_plot_mr_exp_out[[1]] + 
    theme_classic() + 
    xlab("MR estimate per SNP") + 
    ylab("Density")
  density_plot_mr_exp_out$labels$colour <- "MR method"
  density_plot_mr_exp_out$labels$size <- expression(1/SE[IV])
  
  # leave-one-out plot
  leaveoneout_exp_out <- mr_leaveoneout(exp_out)
  leaveoneout_plot_mr_exp_out <- mr_leaveoneout_plot(leaveoneout_exp_out)
  leaveoneout_plot_mr_exp_out <- leaveoneout_plot_mr_exp_out[[1]] + 
    theme_classic() + 
    theme(legend.position = "none") + 
    xlab(paste("MR leave-one-out analysis for", exp_name, "on", out_name))
  leaveoneout_plot_mr_exp_out_snp_rank <- order(leaveoneout_plot_mr_exp_out$data$b[1 : (nrow(leaveoneout_plot_mr_exp_out$data) - 2)])
  leaveoneout_plot_mr_exp_out_snp <- leaveoneout_plot_mr_exp_out$data$SNP[1 : (nrow(leaveoneout_plot_mr_exp_out$data) - 2)][leaveoneout_plot_mr_exp_out_snp_rank]
  leaveoneout_plot_mr_exp_out_snp <- c("All - Inverse variance weighted", "", as.character(leaveoneout_plot_mr_exp_out_snp))
  leaveoneout_plot_mr_exp_out$data$SNP <- as.character(leaveoneout_plot_mr_exp_out$data$SNP)
  leaveoneout_plot_mr_exp_out$data$SNP[which(leaveoneout_plot_mr_exp_out$data$SNP == "All")] <- "All - Inverse variance weighted"
  leaveoneout_plot_mr_exp_out$data$SNP <- ordered(leaveoneout_plot_mr_exp_out$data$SNP, levels = leaveoneout_plot_mr_exp_out_snp)
  
  # patch plot
  patch_plot_mr_exp_out <- wrap_plots(forest_plot_mr_exp_out, scatter_plot_mr_exp_out, leaveoneout_plot_mr_exp_out, funnel_plot_mr_exp_out, ncol = 2, nrow = 2) +
    plot_annotation(tag_levels = "A")
  
  return(list(mr_presso_res_exp_out = mr_presso_res_exp_out, 
              mr_lasso_res_exp_out = mr_lasso_res_exp_out, 
              pleiotropy_exp_out = pleiotropy_exp_out, 
              heterogeneity_exp_out = heterogeneity_exp_out, 
              forest_plot_mr_exp_out = forest_plot_mr_exp_out, 
              scatter_plot_mr_exp_out = scatter_plot_mr_exp_out, 
              funnel_plot_mr_exp_out = funnel_plot_mr_exp_out, 
              density_plot_mr_exp_out = density_plot_mr_exp_out, 
              leaveoneout_plot_mr_exp_out = leaveoneout_plot_mr_exp_out, 
              patch_plot_mr_exp_out = patch_plot_mr_exp_out))
}

mr_pri_ana_gli_per <- mr_pri_ana(exp = gli_exp, 
                                 out = EUR_perio_excl_HCHSSOL_new, 
                                 exp_sample_size = 12496 + 18190)
mr_sup_ana_gli_per <- mr_sup_ana(exp_out = mr_pri_ana_gli_per$exp_out, 
                                 exp_name = "glioma", 
                                 out_name = "periodontitis", 
                                 mr_res_exp_out = mr_pri_ana_gli_per$mr_res_or_exp_out)

mr_pri_ana_gbm_per <- mr_pri_ana(exp = gbm_exp, 
                                 out = EUR_perio_excl_HCHSSOL_new, 
                                 exp_sample_size = 6191 + 18190)
mr_sup_ana_gbm_per <- mr_sup_ana(exp_out = mr_pri_ana_gbm_per$exp_out, 
                                 exp_name = "GBM", 
                                 out_name = "periodontitis", 
                                 mr_res_exp_out = mr_pri_ana_gbm_per$mr_res_or_exp_out)

mr_pri_ana_non_gbm_per <- mr_pri_ana(exp = non_gbm_exp, 
                                     out = EUR_perio_excl_HCHSSOL_new, 
                                     exp_sample_size = 5819 + 18190)
mr_sup_ana_non_gbm_per <- mr_sup_ana(exp_out = mr_pri_ana_non_gbm_per$exp_out, 
                                     exp_name = "non-GBM", 
                                     out_name = "periodontitis", 
                                     mr_res_exp_out = mr_pri_ana_non_gbm_per$mr_res_or_exp_out)

# mediation analysis on immune cell traits
med_ana <- function(exp, exp_sample_size, out, mr_res_or_exp_out) {
  repeat {
    tryCatch(
      {
        withCallingHandlers(
          {
            exp <- clump_data(dat = exp, 
                              clump_kb = 1000, 
                              clump_r2 = 0.1)
            break
          }, 
          message = function(m){if(grepl("Server error", m$message)){stop()}}
        )
      }, 
      error = function(e){message("Retry")}
    )
  }
  exp$samplesize.exposure <- exp_sample_size
  exp$R2 <- exp$beta.exposure ^ 2 / (exp$beta.exposure ^ 2 + exp$samplesize.exposure * exp$se.exposure ^ 2)
  exp$F <- (exp$samplesize.exposure - 2) * exp$R2 / (1 - exp$R2)
  exp <- exp[exp$F >= 10, ]
  
  med_sig_exp <- c()
  exp_med_b <- c()
  exp_med_list <- list()
  for (i in 1391 : 2121) {
    repeat {
      tryCatch({
        withCallingHandlers({
          med_out <- extract_outcome_data(snps = exp$SNP, 
                                          outcomes = paste("ebi-a-GCST9000", i, sep = ""), 
                                          access_token = NULL)
          break
        }, 
        message = function(m){if(grepl("Server error", m$message)){stop()}}
        )
      }, 
      error = function(e){message("Retry")}
      )
    }
    try({
      exp_med <- harmonise_data(exposure_dat = exp, 
                                outcome_dat = med_out, 
                                action = 2)
      exp_med <- exp_med[exp_med$pval.outcome >= 5e-05, ]
      
      exp_med_cd <- cooks.distance(lm(formula = beta.outcome ~ beta.exposure + 0, data = exp_med[exp_med$mr_keep == TRUE, ], weights = se.outcome ^ (-2)))
      exp_med_cd_large <- as.numeric(names(which(exp_med_cd > 4 / nrow(exp_med[exp_med$mr_keep == TRUE, ]))))
      exp_med_sr <- rstudent(lm(formula = beta.outcome ~ beta.exposure + 0, data = exp_med[exp_med$mr_keep == TRUE, ], weights = se.outcome ^ (-2)))
      exp_med_sr_large <- as.numeric(names(which(abs(exp_med_sr) > 2)))
      exp_med <- exp_med[-c(exp_med_cd_large, exp_med_sr_large), ]
      rownames(exp_med) <- 1 : nrow(exp_med)
      
      mr_res_exp_med <- mr(exp_med, method_list = "mr_ivw")
      mr_res_or_exp_med <- generate_odds_ratios(mr_res_exp_med)
      
      if(mr_res_or_exp_med$pval < 0.05) {
        med_sig_exp <- c(med_sig_exp, mr_res_or_exp_med$id.outcome)
        exp_med_b <- c(exp_med_b, mr_res_or_exp_med$b)
        exp_med_list <- c(exp_med_list, list(exp_med))
      }
    })
  }
  
  med_sig_exp_out <- c()
  mr_res_or_exp_med_list <- list()
  mr_res_or_med_out_list <- list()
  mr_res_sup_exp_med_list <- list()
  mr_res_sup_med_out_list <- list()
  med_eff_list <- list()
  med_pro_list <- list()
  for (j in med_sig_exp) {
    repeat {
      tryCatch({
        withCallingHandlers({
          med_exp <- extract_instruments(outcomes = j, 
                                         p1 = 5e-05, 
                                         clump = TRUE, 
                                         p2 = 5e-05, 
                                         r2 = 0.1, 
                                         kb = 1000, 
                                         access_token = NULL)
          break
        }, 
        message = function(m){if(grepl("Server error", m$message)){stop()}}
        )
      }, 
      error = function(e){message("Retry")}
      )
    }
    try({
      med_exp$R2 <- med_exp$beta.exposure ^ 2 / (med_exp$beta.exposure ^ 2 + med_exp$samplesize.exposure * med_exp$se.exposure ^ 2)
      med_exp$F <- (med_exp$samplesize.exposure - 2) * med_exp$R2 / (1 - med_exp$R2)
      med_exp <- med_exp[med_exp$F >= 10, ]
      out_out <- format_data(dat = out, 
                             type = "outcome", 
                             snps = med_exp$SNP, 
                             phenotype_col = "phenotype", 
                             snp_col = "SNP", 
                             beta_col = "BETA", 
                             se_col = "SE", 
                             eaf_col = "FRQ", 
                             effect_allele_col = "A2", 
                             other_allele_col = "A1", 
                             pval_col = "P", 
                             samplesize_col = "N", 
                             chr_col = "CHR", 
                             pos_col = "BP")
      med_out <- harmonise_data(exposure_dat = med_exp, 
                                outcome_dat = out_out, 
                                action = 2)
      med_out <- med_out[med_out$pval.outcome >= 5e-05, ]
      
      med_out_cd <- cooks.distance(lm(formula = beta.outcome ~ beta.exposure + 0, data = med_out[med_out$mr_keep == TRUE, ], weights = se.outcome ^ (-2)))
      med_out_cd_large <- as.numeric(names(which(med_out_cd > 4 / nrow(med_out[med_out$mr_keep == TRUE, ]))))
      med_out_sr <- rstudent(lm(formula = beta.outcome ~ beta.exposure + 0, data = med_out[med_out$mr_keep == TRUE, ], weights = se.outcome ^ (-2)))
      med_out_sr_large <- as.numeric(names(which(abs(med_out_sr) > 2)))
      med_out <- med_out[-c(med_out_cd_large, med_out_sr_large), ]
      rownames(med_out) <- 1 : nrow(med_out)
      
      method_list <- c("mr_egger_regression", 
                       "mr_simple_median", 
                       "mr_weighted_median", 
                       "mr_ivw", 
                       "mr_simple_mode", 
                       "mr_weighted_mode")
      mr_res_med_out <- mr(med_out, method_list = method_list)
      mr_res_or_med_out <- generate_odds_ratios(mr_res_med_out)
      
      if(mr_res_or_med_out$pval[mr_res_or_med_out$method == "Inverse variance weighted"] < 0.05 & exp_med_b[med_sig_exp == j] * mr_res_or_med_out$b[mr_res_or_med_out$method == "Inverse variance weighted"] / mr_res_or_exp_out$b[mr_res_or_exp_out$method == "Inverse variance weighted"] > 0) {
        med_sig_exp_out <- c(med_sig_exp_out, j)
        mr_res_exp_med <- mr(exp_med_list[[which(med_sig_exp == j)]], method_list = method_list)
        mr_res_or_exp_med <- generate_odds_ratios(mr_res_exp_med)
        mr_res_sup_exp_med <- mr_sup_ana(exp_out = exp_med_list[[which(med_sig_exp == j)]], 
                                         exp_name = unique(mr_res_or_exp_med$exposure), 
                                         out_name = strsplit(unique(mr_res_or_exp_med$outcome), " || ", fixed = TRUE)[[1]][1], 
                                         mr_res_exp_out = mr_res_or_exp_med)
        mr_res_sup_med_out <- mr_sup_ana(exp_out = med_out, 
                                         exp_name = strsplit(unique(mr_res_or_exp_med$outcome), " || ", fixed = TRUE)[[1]][1], 
                                         out_name = unique(mr_res_or_med_out$outcome), 
                                         mr_res_exp_out = mr_res_or_med_out)
        med_eff <- medci(mu.x = mr_res_or_exp_med$b[mr_res_or_exp_med$method == "Inverse variance weighted"], 
                         mu.y = mr_res_or_med_out$b[mr_res_or_med_out$method == "Inverse variance weighted"], 
                         se.x = mr_res_or_exp_med$se[mr_res_or_exp_med$method == "Inverse variance weighted"], 
                         se.y = mr_res_or_med_out$se[mr_res_or_med_out$method == "Inverse variance weighted"])
        med_pro <- sapply(med_eff, function(x){x / mr_res_or_exp_out$b[mr_res_or_exp_out$method == "Inverse variance weighted"]})
        mr_res_or_exp_med_list <- c(mr_res_or_exp_med_list, list(mr_res_or_exp_med))
        mr_res_or_med_out_list <- c(mr_res_or_med_out_list, list(mr_res_or_med_out))
        mr_res_sup_exp_med_list <- c(mr_res_sup_exp_med_list, list(mr_res_sup_exp_med))
        mr_res_sup_med_out_list <- c(mr_res_sup_med_out_list, list(mr_res_sup_med_out))
        med_eff_list <- c(med_eff_list, list(med_eff))
        med_pro_list <- c(med_pro_list, list(med_pro))
      }
    })
  }
  
  return(list(med_sig_exp = med_sig_exp, 
              med_sig_exp_out = med_sig_exp_out, 
              mr_res_or_exp_med_list = mr_res_or_exp_med_list, 
              mr_res_or_med_out_list = mr_res_or_med_out_list, 
              mr_res_sup_exp_med_list = mr_res_sup_exp_med_list, 
              mr_res_sup_med_out_list = mr_res_sup_med_out_list, 
              med_eff_list = med_eff_list, 
              med_pro_list = med_pro_list))
}

med_ana_gli_per <- med_ana(exp = gli_exp, 
                           exp_sample_size = 12496 + 18190, 
                           out = EUR_perio_excl_HCHSSOL_new, 
                           mr_res_or_exp_out = mr_pri_ana_gli_per$mr_res_or_exp_out)

med_ana_gbm_per <- med_ana(exp = gbm_exp, 
                           exp_sample_size = 6191 + 18190, 
                           out = EUR_perio_excl_HCHSSOL_new, 
                           mr_res_or_exp_out = mr_pri_ana_gbm_per$mr_res_or_exp_out)

med_ana_non_gbm_per <- med_ana(exp = non_gbm_exp, 
                               exp_sample_size = 5819 + 18190, 
                               out = EUR_perio_excl_HCHSSOL_new, 
                               mr_res_or_exp_out = mr_pri_ana_non_gbm_per$mr_res_or_exp_out)


