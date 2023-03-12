rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
set.seed(123232)
# uncomment the following source code to use different generating mechanism corresponding to,
# from small to large deviation from proxy relevance assumption
source("data_generating_small.R")
#source("data_generating_medium.R")
#source("data_generating_large.R")
require(lava)
library(nleqslv)
library(MASS)
library(KRLS)
library(numDeriv)
library(gmm)




true_beta <- c(-1.9488706, 1)
source("proxicausal.R")

N = 4000
rep_num = 1000

POR_est_result <- c()
PIPW_est_result <- c()
PDR_est_result <- c()


POR_sd_result <- c()
PIPW_sd_result <- c()
PDR_sd_result <- c()

POR_covering <- c()
PIPW_covering <- c()
PDR_covering <- c()
for(rep in 1:rep_num) {
  df <- data_gen(N, para_set)

  source("simu_POR.R")
  POR_est_result <- rbind(POR_est_result, POR_est)
  POR_sd_result <- rbind(POR_sd_result, POR_sd_est)
  POR_covering <- rbind(POR_covering, true_beta <= POR_est + 1.96 * POR_sd_est & true_beta >= POR_est - 1.96 * POR_sd_est)


  source("simu_PIPW.R")
  PIPW_est_result <- rbind(PIPW_est_result, PIPW_est)
  PIPW_sd_result <- rbind(PIPW_sd_result, PIPW_sd_est)
  PIPW_covering <- rbind(PIPW_covering, true_beta <= PIPW_est + 1.96 * PIPW_sd_est & true_beta >= PIPW_est - 1.96 * PIPW_sd_est)


  source("simu_PDR.R")
  PDR_est_result <- rbind(PDR_est_result, PDR_est)
  PDR_sd_result <- rbind(PDR_sd_result, PDR_sd_est)
  PDR_covering <- rbind(PDR_covering, true_beta <= PDR_est + 1.96 * PDR_sd_est &
                          true_beta >= PDR_est - 1.96 * PDR_sd_est)


  print(rep)
}







stats_result <- list(POR_bias = apply(POR_est_result, 2, mean) - true_beta, 
                     PIPW_bias = apply(PIPW_est_result, 2, mean) - true_beta, 
                     PDR_bias = apply(PDR_est_result, 2, mean) - true_beta, 
                     POR_see = apply(POR_est_result, 2, sd), 
                     PIPW_see = apply(PIPW_est_result, 2, sd), 
                     PDR_see = apply(PDR_est_result, 2, sd), 
                     POR_sd = apply(POR_sd_result, 2, mean),
                     PIPW_sd = apply(PIPW_sd_result, 2, mean),
                     PDR_sd = apply(PDR_sd_result, 2, mean),
                     POR_covering = apply(POR_covering, 2, mean),
                     PIPW_covering = apply(PIPW_covering, 2, mean),
                     PDR_covering = apply(PDR_covering, 2, mean)
)

# uncomment the following source code to use different generating mechanism corresponding to,
# from small to large deviation from proxy relevance assumption
save(stats_result, file = paste0("result_small", N, ".rda"))
#save(stats_result, file = paste0("result_medium", N, ".rda"))
#save(stats_result, file = paste0("result_large", N, ".rda"))




