library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(httr)
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

glioma_melin_2017 <- readxl::read_xlsx("glioma_melin_2017.xlsx")
glioma_melin_2017_maf <- c()
for (i in glioma_melin_2017$SNP) {
  q <- GET(paste("http://grch37.rest.ensembl.org/", "variation/", "human/", i, "?pops=1", sep = ""),
           content_type("application/json"))
  message_for_status(q)
  for (j in content(q)$populations) {
    if (j$population == "1000GENOMES:phase_3:EUR" & j$frequency < 0.5) {
      glioma_melin_2017_maf <- c(glioma_melin_2017_maf, j$frequency)
    }
  }
}
glioma_melin_2017$MAF <- glioma_melin_2017_maf
glioma_melin_2017_new <- format_sumstats(glioma_melin_2017,
                                         ref_genome = "GRCh37",
                                         impute_beta = TRUE,
                                         impute_se = TRUE,
                                         allele_flip_frq = FALSE,
                                         dbSNP = "144",
                                         nThread = 4,
                                         save_path = "glioma_melin_2017_new.tsv.gz",
                                         return_data = TRUE)

gbm_melin_2017 <- readxl::read_xlsx("gbm_melin_2017.xlsx")
gbm_melin_2017_maf <- c()
for (i in gbm_melin_2017$SNP) {
  q <- GET(paste("http://grch37.rest.ensembl.org/", "variation/", "human/", i, "?pops=1", sep = ""),
           content_type("application/json"))
  message_for_status(q)
  for (j in content(q)$populations) {
    if (j$population == "1000GENOMES:phase_3:EUR" & j$frequency < 0.5) {
      gbm_melin_2017_maf <- c(gbm_melin_2017_maf, j$frequency)
    }
  }
}
gbm_melin_2017$MAF <- gbm_melin_2017_maf
gbm_melin_2017_new <- format_sumstats(gbm_melin_2017,
                                      ref_genome = "GRCh37",
                                      impute_beta = TRUE,
                                      impute_se = TRUE,
                                      allele_flip_frq = FALSE,
                                      dbSNP = "144",
                                      nThread = 4,
                                      save_path = "gbm_melin_2017_new.tsv.gz",
                                      return_data = TRUE)

non_gbm_melin_2017 <- readxl::read_xlsx("non_gbm_melin_2017.xlsx")
non_gbm_melin_2017_maf <- c()
for (i in non_gbm_melin_2017$SNP) {
  q <- GET(paste("http://grch37.rest.ensembl.org/", "variation/", "human/", i, "?pops=1", sep = ""),
           content_type("application/json"))
  message_for_status(q)
  for (j in content(q)$populations) {
    if (j$population == "1000GENOMES:phase_3:EUR" & j$frequency < 0.5) {
      non_gbm_melin_2017_maf <- c(non_gbm_melin_2017_maf, j$frequency)
    }
  }
}
non_gbm_melin_2017$MAF <- non_gbm_melin_2017_maf
non_gbm_melin_2017_new <- format_sumstats(non_gbm_melin_2017,
                                          ref_genome = "GRCh37",
                                          impute_beta = TRUE,
                                          impute_se = TRUE,
                                          allele_flip_frq = FALSE,
                                          dbSNP = "144",
                                          nThread = 4,
                                          save_path = "non_gbm_melin_2017_new.tsv.gz",
                                          return_data = TRUE)

EUR_perio_excl_HCHSSOL <- data.table::fread("EUR_perio_excl_HCHSSOL.txt")
EUR_perio_excl_HCHSSOL_new <- tidyr::separate(data = EUR_perio_excl_HCHSSOL,
                                              col = "MarkerName",
                                              into = c("CHR", "BP"),
                                              sep = ":")
EUR_perio_excl_HCHSSOL_new$MAF <- EUR_perio_excl_HCHSSOL_new$MAC / (2 * EUR_perio_excl_HCHSSOL_new$N)
EUR_perio_excl_HCHSSOL_new <- dplyr::rename(EUR_perio_excl_HCHSSOL_new,
                                            CHR = CHR,
                                            BP = BP,
                                            A1 = Allele2,
                                            A2 = Allele1,
                                            BETA = Effect,
                                            SE = StdErr,
                                            P = `P-value`,
                                            N = N,
                                            MAF = MAF)
EUR_perio_excl_HCHSSOL_new <- dplyr::select(EUR_perio_excl_HCHSSOL_new, CHR, BP, A1, A2, BETA, SE, P, N, MAF)
EUR_perio_excl_HCHSSOL_new <- format_sumstats(EUR_perio_excl_HCHSSOL_new,
                                              ref_genome = "GRCh37",
                                              allele_flip_frq = FALSE,
                                              dbSNP = "144",
                                              nThread = 4,
                                              save_path = "EUR_perio_excl_HCHSSOL_new.tsv.gz",
                                              return_data = TRUE)


