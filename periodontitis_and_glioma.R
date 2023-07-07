# devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
# devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
library(MendelianRandomization)
library(ggplot2)

# periodontitis → glioma
# extract SNP associated with exposure
repeat{
  tryCatch({
    withCallingHandlers({
      per_exp <- extract_instruments(outcomes = "ebi-a-GCST003484", 
                                     p1 = 5e-06, 
                                     clump = TRUE, 
                                     p2 = 5e-06, 
                                     r2 = 0.01, 
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
repeat{
  tryCatch({
    withCallingHandlers({
      per_n <- ieugwasr::gwasinfo("ebi-a-GCST003484")$sample_size
      break
    }, 
    message = function(m){if(grepl("Server error", m$message)){stop()}}
    )
  }, 
  error = function(e){message("Retry")}
  )
}
per_exp$R2 <- per_exp$beta.exposure ^ 2 / (per_exp$beta.exposure ^ 2 + per_n * per_exp$se.exposure ^ 2)
per_exp$F <- (per_n - 2) * per_exp$R2 / (1 - per_exp$R2)
per_exp <- subset(per_exp, F >= 10)

# extract outcome
repeat{
  tryCatch({
    withCallingHandlers({
      gli_out <- extract_outcome_data(snps = per_exp$SNP, 
                                      outcomes = "ieu-a-1013", 
                                      access_token = NULL)
      break
    }, 
    message = function(m){if(grepl("Server error", m$message)){stop()}}
    )
  }, 
  error = function(e){message("Retry")}
  )
}

# combine data
per_gli <- harmonise_data(exposure_dat = per_exp, 
                          outcome_dat = gli_out, 
                          action = 2)

# exclude SNP associated with outcome
per_gli <- subset(per_gli, pval.outcome >= 5e-06)

# MR-LASSO
mr_input_per_gli <- mr_input(bx = per_gli$beta.exposure, 
                             bxse = per_gli$se.exposure, 
                             by = per_gli$beta.outcome, 
                             byse = per_gli$se.outcome, 
                             exposure = "periodontitis", 
                             outcome = "glioma", 
                             snps = per_gli$SNP)
mrlasso_res_per_gli <- mr_lasso(mr_input_per_gli, lambda = 0.5)
per_gli <- per_gli[per_gli$SNP %in% mrlasso_res_per_gli@ValidSNPs, ]

# MR-PRESSO
mrpresso_res_per_gli <- mr_presso(BetaOutcome = "beta.outcome", 
                                  BetaExposure = "beta.exposure", 
                                  SdOutcome = "se.outcome", 
                                  SdExposure = "se.exposure", 
                                  OUTLIERtest = TRUE, 
                                  DISTORTIONtest = TRUE, 
                                  data = per_gli, 
                                  NbDistribution = 1000, 
                                  SignifThreshold = 0.05, 
                                  seed = 123)

# MR
method_list <- c("mr_egger_regression", 
                 "mr_weighted_median", 
                 "mr_ivw", 
                 "mr_simple_mode", 
                 "mr_weighted_mode")
mr_res_per_gli <- mr(per_gli, method_list = method_list)
mr_res_or_per_gli <- generate_odds_ratios(mr_res_per_gli)

# pleiotropy test
pleiotropy_per_gli <- mr_pleiotropy_test(per_gli)

# heterogeneity test
heterogeneity_per_gli <- mr_heterogeneity(per_gli)

# forest plot
single_snp_per_gli <- mr_singlesnp(per_gli, all_method = method_list)
forest_plot_mr_per_gli <- mr_forest_plot(single_snp_per_gli)
forest_plot_mr_per_gli <- forest_plot_mr_per_gli$`ebi-a-GCST003484.ieu-a-1013` + 
  theme_classic() + 
  theme(legend.position = "none") + 
  xlab("MR effect size for periodontitis on glioma")

# scatter plot
scatter_plot_mr_per_gli <- mr_scatter_plot(mr_res_per_gli, per_gli)
scatter_plot_mr_per_gli <- scatter_plot_mr_per_gli$`ebi-a-GCST003484.ieu-a-1013` + 
  theme_classic() + 
  theme(legend.direction = "horizontal", 
        legend.position = "top", 
        legend.background = element_blank()) + 
  xlab("SNP effect on periodontitis") + 
  ylab("SNP effect on glioma")
scatter_plot_mr_per_gli <- scatter_plot_mr_per_gli + guides(col = guide_legend(nrow = 1))

# funnel plot
funnel_plot_mr_per_gli <- mr_funnel_plot(single_snp_per_gli)
funnel_plot_mr_per_gli <- funnel_plot_mr_per_gli$`ebi-a-GCST003484.ieu-a-1013` + 
  theme_classic() + 
  theme(legend.direction = "horizontal", 
        legend.position = "top", 
        legend.background = element_blank())
funnel_plot_mr_per_gli <- funnel_plot_mr_per_gli + guides(col = guide_legend(nrow = 1))

# dnesity plot
density_plot_mr_per_gli <- mr_density_plot(single_snp_per_gli, mr_res_per_gli)
density_plot_mr_per_gli <- density_plot_mr_per_gli$`ebi-a-GCST003484.ieu-a-1013` + 
  theme_classic()

# leave-one-out plot
leaveoneout_per_gli <- mr_leaveoneout(per_gli)
leaveoneout_plot_mr_per_gli <- mr_leaveoneout_plot(leaveoneout_per_gli)
leaveoneout_plot_mr_per_gli <- leaveoneout_plot_mr_per_gli$`ebi-a-GCST003484.ieu-a-1013` + 
  theme_classic() + 
  theme(legend.position = "none") + 
  xlab("MR leave-one-out sensitivity analysis for periodontitis on glioma")

# glioma → periodontitis
# extract SNP associated with exposure
repeat{
  tryCatch({
    withCallingHandlers({
      gli_exp <- extract_instruments(outcomes = "ieu-a-1013", 
                                     p1 = 5e-06, 
                                     clump = TRUE, 
                                     p2 = 5e-06, 
                                     r2 = 0.01, 
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
repeat{
  tryCatch({
    withCallingHandlers({
      gli_n <- ieugwasr::gwasinfo("ieu-a-1013")$sample_size
      break
    }, 
    message = function(m){if(grepl("Server error", m$message)){stop()}}
    )
  }, 
  error = function(e){message("Retry")}
  )
}
gli_exp$R2 <- gli_exp$beta.exposure ^ 2 / (gli_exp$beta.exposure ^ 2 + gli_n * gli_exp$se.exposure ^ 2)
gli_exp$F <- (gli_n - 2) * gli_exp$R2 / (1 - gli_exp$R2)
gli_exp <- subset(gli_exp, F >= 10)

# extract outcome
repeat{
  tryCatch({
    withCallingHandlers({
      per_out <- extract_outcome_data(snps = gli_exp$SNP, 
                                      outcomes = "ebi-a-GCST003484", 
                                      access_token = NULL)
      break
    }, 
    message = function(m){if(grepl("Server error", m$message)){stop()}}
    )
  }, 
  error = function(e){message("Retry")}
  )
}

# combine data
gli_per <- harmonise_data(exposure_dat = gli_exp, 
                          outcome_dat = per_out, 
                          action = 2)

# exclude SNP associated with outcome
gli_per <- subset(gli_per, pval.outcome >= 5e-06)

# MR-LASSO
mr_input_gli_per <- mr_input(bx = gli_per$beta.exposure, 
                             bxse = gli_per$se.exposure, 
                             by = gli_per$beta.outcome, 
                             byse = gli_per$se.outcome, 
                             exposure = "glioma", 
                             outcome = "periodontitis", 
                             snps = gli_per$SNP)
mrlasso_res_gli_per <- mr_lasso(mr_input_gli_per, lambda = 0.5)
gli_per <- gli_per[gli_per$SNP %in% mrlasso_res_gli_per@ValidSNPs, ]

# MR-PRESSO
mrpresso_res_gli_per <- mr_presso(BetaOutcome = "beta.outcome", 
                                  BetaExposure = "beta.exposure", 
                                  SdOutcome = "se.outcome", 
                                  SdExposure = "se.exposure", 
                                  OUTLIERtest = TRUE, 
                                  DISTORTIONtest = TRUE, 
                                  data = gli_per, 
                                  NbDistribution = 1000, 
                                  SignifThreshold = 0.05, 
                                  seed = 123)

# MR
method_list <- c("mr_egger_regression", 
                 "mr_weighted_median", 
                 "mr_ivw", 
                 "mr_simple_mode", 
                 "mr_weighted_mode")
mr_res_gli_per <- mr(gli_per, method_list = method_list)
mr_res_or_gli_per <- generate_odds_ratios(mr_res_gli_per)

# pleiotropy test
pleiotropy_gli_per <- mr_pleiotropy_test(gli_per)

# heterogeneity test
heterogeneity_gli_per <- mr_heterogeneity(gli_per)

# forest plot
single_snp_gli_per <- mr_singlesnp(gli_per, all_method = method_list)
forest_plot_mr_gli_per <- mr_forest_plot(single_snp_gli_per)
forest_plot_mr_gli_per <- forest_plot_mr_gli_per$`ieu-a-1013.ebi-a-GCST003484` + 
  theme_classic() + 
  theme(legend.position = "none") + 
  xlab("MR effect size for glioma on periodontitis")

# scatter plot
scatter_plot_mr_gli_per <- mr_scatter_plot(mr_res_gli_per, gli_per)
scatter_plot_mr_gli_per <- scatter_plot_mr_gli_per$`ieu-a-1013.ebi-a-GCST003484` + 
  theme_classic() + 
  theme(legend.direction = "horizontal", 
        legend.position = "top", 
        legend.background = element_blank()) + 
  xlab("SNP effect on glioma") + 
  ylab("SNP effect on periodontitis")
scatter_plot_mr_gli_per <- scatter_plot_mr_gli_per + guides(col = guide_legend(nrow = 1))

# funnel plot
funnel_plot_mr_gli_per <- mr_funnel_plot(single_snp_gli_per)
funnel_plot_mr_gli_per <- funnel_plot_mr_gli_per$`ieu-a-1013.ebi-a-GCST003484` + 
  theme_classic() + 
  theme(legend.direction = "horizontal", 
        legend.position = "top", 
        legend.background = element_blank())
funnel_plot_mr_gli_per <- funnel_plot_mr_gli_per + guides(col = guide_legend(nrow = 1))

# density plot
density_plot_mr_gli_per <- mr_density_plot(single_snp_gli_per, mr_res_gli_per)
density_plot_mr_gli_per <- density_plot_mr_gli_per$`ieu-a-1013.ebi-a-GCST003484` + 
  theme_classic()

# leave-one-out plot
leaveoneout_gli_per <- mr_leaveoneout(gli_per)
leaveoneout_plot_mr_gli_per <- mr_leaveoneout_plot(leaveoneout_gli_per)
leaveoneout_plot_mr_gli_per <- leaveoneout_plot_mr_gli_per$`ieu-a-1013.ebi-a-GCST003484` + 
  theme_classic() + 
  theme(legend.position = "none") + 
  xlab("MR leave-one-out sensitivity analysis for glioma on periodontitis")


