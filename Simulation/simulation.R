rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
set.seed(30)
source("data_generating.R")
require(lava)
library(nleqslv)
library(MASS)
library(KRLS)
library(numDeriv)
library(gmm)


N = 1000000
true_beta_result <- c()



true_beta <- c(-1.9488706, 1)
source("proxicausal.R")

N = 500
rep_num = 1000
round = 100


POR_est_result <- c()
PIPW_est_result <- c()
POR_WOR_est_result <- c()
PIPW_WIPW_est_result <- c()
PDR_est_result <- c()
PDR_WOR_est_result <- c()
PDR_WIPW_est_result <- c()
PDR_BW_est_result <- c()
DR_est_result <- c()

POR_sd_result <- c()
PIPW_sd_result <- c()
POR_WOR_sd_result <- c()
PIPW_WIPW_sd_result <- c()
PDR_sd_result <- c()
PDR_WOR_sd_result <- c()
PDR_WIPW_sd_result <- c()
PDR_BW_sd_result <- c()
DR_sd_result <- c()

POR_covering <- c()
PIPW_covering <- c()
POR_WOR_covering <- c()
PIPW_WIPW_covering <- c()
PDR_covering <- c()
PDR_WOR_covering <- c()
PDR_WIPW_covering <- c()
PDR_BW_covering <- c()
DR_covering <- c()

POR_converging_result <- c()
PIPW_converging_result <- c()
POR_WOR_converging_result <- c()
PIPW_WIPW_converging_result <- c()
PDR_converging_result <- c()
PDR_WOR_converging_result <- c()
PDR_WIPW_converging_result <- c()
PDR_BW_converging_result <- c()
DR_converging_result <- c()
for(rep in 1:rep_num) {
  df <- data_gen(N, para_set)

  source("simu_POR.R")
  POR_est_result <- rbind(POR_est_result, POR_est)
  POR_sd_result <- rbind(POR_sd_result, POR_sd_est)
  POR_covering <- rbind(POR_covering, true_beta <= POR_est + 1.96 * POR_sd_est & true_beta >= POR_est - 1.96 * POR_sd_est)
  POR_converging_result <- c(POR_converging_result, POR_converging)
  

  source("simu_PIPW.R")
  PIPW_est_result <- rbind(PIPW_est_result, PIPW_est)
  PIPW_sd_result <- rbind(PIPW_sd_result, PIPW_sd_est)
  PIPW_covering <- rbind(PIPW_covering, true_beta <= PIPW_est + 1.96 * PIPW_sd_est & true_beta >= PIPW_est - 1.96 * PIPW_sd_est)
  PIPW_converging_result <- c(PIPW_converging_result, PIPW_converging)

  source("simu_PDR.R")
  PDR_est_result <- rbind(PDR_est_result, PDR_est)
  PDR_sd_result <- rbind(PDR_sd_result, PDR_sd_est)
  PDR_covering <- rbind(PDR_covering, true_beta <= PDR_est + 1.96 * PDR_sd_est &
                          true_beta >= PDR_est - 1.96 * PDR_sd_est)
  PDR_converging_result <- c(PDR_converging_result, POR_converging & PIPW_converging)
  

  source("simu_POR_WOR.R")
  POR_WOR_est_result <- rbind(POR_WOR_est_result, POR_WOR_est)
  POR_WOR_sd_result <- rbind(POR_WOR_sd_result, POR_WOR_sd_est)
  POR_WOR_covering <- rbind(POR_WOR_covering, true_beta <= POR_WOR_est + 1.96 * POR_WOR_sd_est &
                              true_beta >= POR_WOR_est - 1.96 * POR_WOR_sd_est)
  POR_WOR_converging_result <- c(POR_WOR_converging_result, POR_WOR_converging)
  

  source("simu_PDR_WOR.R")
  PDR_WOR_est_result <- rbind(PDR_WOR_est_result, PDR_WOR_est)
  PDR_WOR_sd_result <- rbind(PDR_WOR_sd_result, PDR_WOR_sd_est)
  PDR_WOR_covering <- rbind(PDR_WOR_covering, true_beta <= PDR_WOR_est + 1.96 * PDR_WOR_sd_est &
                              true_beta >= PDR_WOR_est - 1.96 * PDR_WOR_sd_est)
  PDR_WOR_converging_result <- c(PDR_WOR_converging_result, POR_WOR_converging & PIPW_converging)
  

  source("simu_PIPW_WIPW.R")
  PIPW_WIPW_est_result <- rbind(PIPW_WIPW_est_result, PIPW_WIPW_est)
  PIPW_WIPW_sd_result <- rbind(PIPW_WIPW_sd_result, PIPW_WIPW_sd_est)
  PIPW_WIPW_covering <- rbind(PIPW_WIPW_covering, true_beta <= PIPW_WIPW_est + 1.96 * PIPW_WIPW_sd_est &
                                true_beta >= PIPW_WIPW_est - 1.96 * PIPW_WIPW_sd_est)
  PIPW_WIPW_converging_result <- c(PIPW_WIPW_converging_result, PIPW_WIPW_converging)
  

  source("simu_PDR_WIPW.R")
  PDR_WIPW_est_result <- rbind(PDR_WIPW_est_result, PDR_WIPW_est)
  PDR_WIPW_sd_result <- rbind(PDR_WIPW_sd_result, PDR_WIPW_sd_est)
  PDR_WIPW_covering <- rbind(PDR_WIPW_covering, true_beta <= PDR_WIPW_est + 1.96 * PDR_WIPW_sd_est &
                               true_beta >= PDR_WIPW_est - 1.96 * PDR_WIPW_sd_est)
  PDR_WIPW_converging_result <- c(PDR_WIPW_converging_result, POR_converging & PIPW_WIPW_converging)
  

  source("simu_PDR_BW.R")
  PDR_BW_est_result <- rbind(PDR_BW_est_result, PDR_BW_est)
  PDR_BW_sd_result <- rbind(PDR_BW_sd_result, PDR_BW_sd_est)
  PDR_BW_covering <- rbind(PDR_BW_covering, true_beta <= PDR_BW_est + 1.96 * PDR_BW_sd_est &
                               true_beta >= PDR_BW_est - 1.96 * PDR_BW_sd_est)
  PDR_BW_converging_result <- c(PDR_BW_converging_result, POR_WOR_converging & PIPW_WIPW_converging)
  wrong_t1_recording <- rbind(wrong_t1_recording, t1)
  
  source("simu_DR.R")
  DR_est_result <- rbind(DR_est_result, DR_est)
  DR_sd_result <- rbind(DR_sd_result, DR_sd_est)
  DR_covering <- rbind(DR_covering, true_beta <= DR_est + 1.96 * DR_sd_est &
                              true_beta >= DR_est - 1.96 * DR_sd_est)
  DR_converging_result <- c(DR_converging_result, TRUE)
  print(rep)
}






running_result <- list(POR_est_result = POR_est_result,
                       POR_WOR_est_result = POR_WOR_est_result,
                       PIPW_est_result = PIPW_est_result,
                       PIPW_WIPW_est_result = PIPW_WIPW_est_result,
                       PDR_est_result = PDR_est_result,
                       PDR_WOR_est_result = PDR_WOR_est_result,
                       PDR_WIPW_est_result = PDR_WIPW_est_result,
                       PDR_BW_est_result = PDR_BW_est_result,
                       DR_est_result = DR_est_result,
                       POR_sd_result = POR_sd_result,
                       POR_WOR_sd_result = POR_WOR_sd_result,
                       PIPW_sd_result = PIPW_sd_result,
                       PIPW_WIPW_sd_result = PIPW_WIPW_sd_result,
                       PDR_sd_result = PDR_sd_result,
                       PDR_WOR_sd_result = PDR_WOR_sd_result,
                       PDR_WIPW_sd_result = PDR_WIPW_sd_result,
                       PDR_BW_sd_result = PDR_BW_sd_result,
                       DR_sd_result = DR_sd_result,
                       POR_covering = POR_covering,
                       POR_WOR_covering = POR_WOR_covering,
                       PIPW_covering = PIPW_covering,
                       PIPW_WIPW_covering = PIPW_WIPW_covering,
                       PDR_covering = PDR_covering,
                       PDR_WOR_covering = PDR_WOR_covering,
                       PDR_WIPW_covering = PDR_WIPW_covering,
                       PDR_BW_covering = PDR_BW_covering,
                       DR_covering = DR_covering,
                       POR_converging_result = POR_converging_result,
                       POR_WOR_converging_result = POR_WOR_converging_result,
                       PIPW_converging_result = PIPW_converging_result,
                       PIPW_WIPW_converging_result = PIPW_WIPW_converging_result,
                       PDR_converging_result = PDR_converging_result,
                       PDR_WOR_converging_result = PDR_WOR_converging_result,
                       PDR_WIPW_converging_result = PDR_WIPW_converging_result,
                       PDR_BW_converging_result = PDR_BW_converging_result,
                       DR_converging_result = DR_converging_result)


stats_result <- with(running_result, list(POR_bias = apply(POR_est_result, 2, mean) - true_beta, 
                     POR_WOR_bias = apply(POR_WOR_est_result, 2, mean) - true_beta, 
                     PIPW_bias = apply(PIPW_est_result, 2, mean) - true_beta, 
                     PIPW_WIPW_bias = apply(PIPW_WIPW_est_result, 2, mean) - true_beta,
                     PDR_bias = apply(PDR_est_result, 2, mean) - true_beta, 
                     PDR_WOR_bias = apply(PDR_WOR_est_result, 2, mean) - true_beta, 
                     PDR_WIPW_bias = apply(PDR_WIPW_est_result, 2, mean) - true_beta, 
                     PDR_BW_bias = apply(PDR_BW_est_result, 2, mean) - true_beta,
                     DR_bias = apply(DR_est_result, 2, mean) - true_beta,
                     POR_see = apply(POR_est_result, 2, sd), 
                     POR_WOR_see = apply(POR_WOR_est_result, 2, sd), 
                     PIPW_see = apply(PIPW_est_result, 2, sd), 
                     PIPW_WIPW_see = apply(PIPW_WIPW_est_result, 2, sd),
                     PDR_see = apply(PDR_est_result, 2, sd), 
                     PDR_WOR_see = apply(PDR_WOR_est_result, 2, sd), 
                     PDR_WIPW_see = apply(PDR_WIPW_est_result, 2, sd),
                     PDR_BW_see = apply(PDR_BW_est_result, 2, sd),
                     DR_see = apply(DR_est_result, 2, sd),
                     POR_sd = apply(POR_sd_result, 2, mean),
                     POR_WOR_sd = apply(POR_WOR_sd_result, 2, mean),
                     PIPW_sd = apply(PIPW_sd_result, 2, mean),
                     PIPW_WIPW_sd = apply(PIPW_WIPW_sd_result, 2, mean),
                     PDR_sd = apply(PDR_sd_result, 2, mean),
                     PDR_WOR_sd = apply(PDR_WOR_sd_result, 2, mean),
                     PDR_WIPW_sd = apply(PDR_WIPW_sd_result, 2, mean),
                     PDR_BW_sd = apply(PDR_BW_sd_result, 2, mean),
                     DR_sd = apply(DR_sd_result, 2, mean),
                     POR_covering = apply(POR_covering, 2, mean),
                     POR_WOR_covering = apply(POR_WOR_covering, 2, mean),
                     PIPW_covering = apply(PIPW_covering, 2, mean),
                     PIPW_WIPW_covering = apply(PIPW_WIPW_covering, 2, mean),
                     PDR_covering = apply(PDR_covering, 2, mean),
                     PDR_WOR_covering = apply(PDR_WOR_covering, 2, mean),
                     PDR_WIPW_covering = apply(PDR_WIPW_covering, 2, mean),
                     PDR_BW_covering = apply(PDR_BW_covering, 2, mean),
                     DR_covering = apply(DR_covering, 2, mean),
                     POR_converging = mean(POR_converging_result),
                     POR_WOR_converging = mean(POR_WOR_converging_result),
                     PIPW_converging = mean(PIPW_converging_result),
                     PIPW_WIPW_converging = mean(PIPW_WIPW_converging_result),
                     PDR_converging = mean(PDR_converging_result),
                     PDR_WOR_converging = mean(PDR_WOR_converging_result),
                     PDR_WIPW_converging = mean(PDR_WIPW_converging_result),
                     PDR_BW_converging = mean(PDR_BW_converging_result),
                     DR_converging = mean(DR_converging_result)))

table_result <- with(running_result, 
                     data.frame(
                       Bias = formatC(signif(c(apply(POR_est_result, 2, mean) - true_beta, 
                                apply(POR_WOR_est_result, 2, mean) - true_beta, 
                                apply(PIPW_est_result, 2, mean) - true_beta, 
                                apply(PIPW_WIPW_est_result, 2, mean) - true_beta,
                                apply(PDR_est_result, 2, mean) - true_beta, 
                                apply(PDR_WOR_est_result, 2, mean) - true_beta, 
                                apply(PDR_WIPW_est_result, 2, mean) - true_beta, 
                                apply(PDR_BW_est_result, 2, mean) - true_beta,
                                apply(DR_est_result, 2, mean) - true_beta) * 1000, 3), digits=3,format="fg", flag="#"),
                       SEE = formatC(signif(c(apply(POR_est_result, 2, sd), 
                                      apply(POR_WOR_est_result, 2, sd), 
                                      apply(PIPW_est_result, 2, sd), 
                                      apply(PIPW_WIPW_est_result, 2, sd),
                                      apply(PDR_est_result, 2, sd), 
                                      apply(PDR_WOR_est_result, 2, sd), 
                                      apply(PDR_WIPW_est_result, 2, sd), 
                                      apply(PDR_BW_est_result, 2, sd),
                                      apply(DR_est_result, 2, sd)) * 1000, 3), digits=3,format="fg", flag="#"),
                       SD = formatC(signif(c(apply(POR_sd_result, 2, mean), 
                                     apply(POR_WOR_sd_result, 2, mean), 
                                     apply(PIPW_sd_result, 2, mean), 
                                     apply(PIPW_WIPW_sd_result, 2, mean),
                                     apply(PDR_sd_result, 2, mean), 
                                     apply(PDR_WOR_sd_result, 2, mean), 
                                     apply(PDR_WIPW_sd_result, 2, mean), 
                                     apply(PDR_BW_sd_result, 2, mean),
                                     apply(DR_sd_result, 2, mean)) * 1000, 3), digits=3,format="fg", flag="#"),
                       CP = formatC(round(c(apply(POR_covering, 2, mean),
                                     apply(POR_WOR_covering, 2, mean),
                                     apply(PIPW_covering, 2, mean),
                                     apply(PIPW_WIPW_covering, 2, mean),
                                     apply(PDR_covering, 2, mean),
                                     apply(PDR_WOR_covering, 2, mean),
                                     apply(PDR_WIPW_covering, 2, mean),
                                     apply(PDR_BW_covering, 2, mean),
                                     apply(DR_covering, 2, mean)) * 100, 1), digits=3,format="fg", flag="#")
                       )
              )

save(running_result, stats_result, file = paste0("result", N, ".rda"))



caption = paste0("Simulation results of POR, PIPW and PDR estimators. We report bias ($\\times 10^{-3}$), empirical standard error (SEE) ($\\times 10^{-3}$), average estimated standard error (SD) ($\\times 10^{-3}$), and coverage probability of 95\\% confidence intervals (95\\% CP) of POR with correctly specified outcome confounding bridge functions ($\\widehat \\beta_{\\text{POR}}$), POR with incorrectly specified outcome confounding bridge functions ($\\widehat \\beta_{\\text{POR, WOR}}$), PIPW with correctly specified treatment confounding bridge functions ($\\widehat \\beta_{\\text{PIPW}}$), PIPW with incorrectly specified treatment confounding bridge functions ($\\widehat \\beta_{\\text{PIPW, WIPW}}$), PDR with both outcome and treatment confounding bridge functions correctly specified ($\\widehat \\beta_{\\text{PDR}}$), PDR with incorrectly specified outcome confounding bridge functions ($\\widehat \\beta_{\\text{PDR, WOR}}$), PDR with incorrectly specified treatment confounding bridge functions ($\\widehat \\beta_{\\text{PDR, WIPW}}$), PDR with both outcome and treatment confounding bridge functions wrongly put ($\\widehat \\beta_{\\text{PDR, BW}}$) and a standard doubly robust estimator ($\\widehat \\beta_{\\text{DR}}$) for $\\beta$ in model \\eqref{eq:simumsmm}, for sample size $N = ", N, "$ and $B = 1000$ Monte Carlo samples.")

label = paste0("tab:simubeta", N)
colnames(table_result) = c("Bias", "SEE", "SD", "95\\% CP")
row.names(table_result) = c("$\\widehat \\beta_{0, \\text{POR}}$",
                            "$\\widehat \\beta_{1, \\text{POR}}$",
                            "$\\widehat \\beta_{0, \\text{POR, WOR}}$",
                            "$\\widehat \\beta_{1, \\text{POR, WOR}}$",
                            "$\\widehat \\beta_{0, \\text{PIPW}}$",
                            "$\\widehat \\beta_{1, \\text{PIPW}}$",
                            "$\\widehat \\beta_{0, \\text{PIPW, WIPW}}$",
                            "$\\widehat \\beta_{1, \\text{PIPW, WIPW}}$",
                            "$\\widehat \\beta_{0, \\text{PDR}}$",
                            "$\\widehat \\beta_{1, \\text{PDR}}$",
                            "$\\widehat \\beta_{0, \\text{PDR, WOR}}$",
                            "$\\widehat \\beta_{1, \\text{PDR, WOR}}$",
                            "$\\widehat \\beta_{0, \\text{PDR, WIPW}}$",
                            "$\\widehat \\beta_{1, \\text{PDR, WIPW}}$",
                            "$\\widehat \\beta_{0, \\text{PDR, BW}}$",
                            "$\\widehat \\beta_{1, \\text{PDR, BW}}$",
                            "$\\widehat \\beta_{0, \\text{DR}}$",
                            "$\\widehat \\beta_{1, \\text{DR}}$")


print(xtable(table_result, caption = caption, label = label, align = "|lcccc|"), caption.placement = "top", sanitize.text.function = identity)

