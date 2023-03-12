# rm(list = ls())
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# #set.seed(1234567)
# source("data_generating.R")
# require(lava)
# N = 1000000
# true_beta_result <- c()
# 
# #for(b in 1:100) {
# df <- data_gen(N, para_set)
# 
# true_beta_0 <- mean(intervened_data_gen(N, para_set, a = c(0, 0))$Y)
# 
# true_beta_1 <- mean(intervened_data_gen(N, para_set, a = c(1, 0))$Y) - true_beta_0
# 
# true_beta_2 <- mean(intervened_data_gen(N, para_set, a = c(0, 1))$Y) - true_beta_0
# 
# true_beta_3 <- mean(intervened_data_gen(N, para_set, a = c(1, 1))$Y) - true_beta_0
# 
# true_beta <- c(true_beta_0, true_beta_1, true_beta_2, true_beta_3)
# 
# #true_beta_result <- rbind(true_beta_result, true_beta)
# #}
# #colMeans(true_beta_result)
# 
# 
# true_beta <- c(-1.9488706, 1)
# source("proxicausal.R")
# 
# N = 4000
# rep_num = 1000
# t0_true <- c(0.35, 0.7, 0.7, 0.4)
# t1_true <- c(0.7, 0.7, 0.7, -0.7, -0.7, 0.4, 0.2)
# t0_result <- c()
# t1_result <- c()
# est_result <- c()
# sd_result <- c()
# PIPW_covering <- c()
# for(rep in 1:rep_num) {
#   df <- data_gen(N, para_set)
#   
PIPW_WIPW_df <- df
PIPW_WIPW_df$Z0 <- abs(PIPW_WIPW_df$Z0)#^(1/2) + 1
PIPW_WIPW_df$Z1 <- abs(PIPW_WIPW_df$Z1)#^(1/2) + 1


temp_t1 <- t1
temp_t0 <- t0
temp_nuisance_t1_IF <- nuisance_t1_IF
temp_nuisance_t0_IF <- nuisance_t0_IF
temp_q0_Y <- q0_Y
temp_q0_Z <- q0_Z
temp_q0_W <- q0_W
temp_q1_Y <- q1_Y
temp_q1_Z0 <- q1_Z0
temp_q1_Z1 <- q1_Z1
temp_q1_W <- q1_W

inioptim_t0 <- c(0, 0, 0, 0)
q0_Y <- with(PIPW_WIPW_df, matrix(rep(c(0, 1, 0, 0), each = length(A0)), nrow = length(A0), ncol = 4))
q0_Z <- with(PIPW_WIPW_df, (-1)^(1 - A0) * cbind(1, A0, Z0, X0))
q0_W <- with(PIPW_WIPW_df, (-1)^(1 - A0) * cbind(1, A0, W0, X0))

q0link <- optim(par = inioptim_t0, fn = MSE_func,
                bridge_func = q0bridge, Y = q0_Y, W = q0_W, Z = q0_Z,
                method = "BFGS", hessian = FALSE)
t0 <- q0link$par



q1_Y <- with(PIPW_WIPW_df, matrix(rep(c(0, 1, 0, 0, 0, 0, 0), each = length(A0)), nrow = length(A0), ncol = 7))
q1_Z0 <- with(PIPW_WIPW_df, (-1)^(1 - A1) * cbind(1, A1, A0, Z1, Z0, X1, X0))
q1_Z1 <- with(PIPW_WIPW_df, (-1)^(1 - A1) * cbind(1, A1, A0, Z1, Z0, X1, X0))

q1_W <- with(PIPW_WIPW_df, (-1)^(1 - A1) * cbind(1, A1, A0, W1, W0, X1, X0))


inioptim_t1 <- c(0, 0, 0, 0, 0, 0, 0)
q1link <- optim(par = inioptim_t1, fn = MSE_func_q1,
                bridge_func = q1bridge, Y = q1_Y, W = q1_W, Z0 = q1_Z0, Z1 = q1_Z1, q0 = c(exp(q0_Z %*% t0)), 
                method = "BFGS", hessian = FALSE)
t1 <- q1link$par


PIPW_WIPW_work_df <- data.frame(Y = with(PIPW_WIPW_df, rep(Y, 4)),
                           A = rep(c(0, 1, 1, 2), each = N))

PIPW_WIPW_weights <- with(PIPW_WIPW_df, c((1 - A1) * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)), 
                                (1 - A1) * A0 * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)),
                                A1 * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)), 
                                A1 * A0 * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1))))

PIPW_WIPW_fit <- lm(Y ~ A, weights = PIPW_WIPW_weights, data = PIPW_WIPW_work_df)
weighted_denom_inv <- solve(t(model.matrix(PIPW_WIPW_fit) * PIPW_WIPW_weights) %*% model.matrix(PIPW_WIPW_fit))
PIPW_WIPW_IF <- PIPW_WIPW_weights * model.matrix(PIPW_WIPW_fit) %*% weighted_denom_inv * residuals(PIPW_WIPW_fit)


PIPW_WIPW_est <- PIPW_WIPW_fit$coefficients


t0_res <- with(PIPW_WIPW_df, c(1 + exp(q1_Z1 %*% t1)) * q1_W -
                 matrix(rep(c(0, 1, 0, 0, 0, 0, 0), each = length(A0)), nrow = length(A0), ncol = 7))
partial_t0_mat <- partial_t0(t1 = t1, q1_Y = q1_Y, q1_W = q1_W, q1_Z0 = q1_Z0, q1_Z1 = q1_Z1, q0_Z = q0_Z, q0 = c(exp(q0_Z %*% t0)))

t0_IF <- q0_IF(t0, Y = q0_Y, W = q0_W, Z = q0_Z)

t1_IF <- q1_IF(t1, Y = q1_Y, W = q1_W, Z0 = q1_Z0, Z1 = q1_Z1, q0 = c(exp(q0_Z %*% t0)), t0_IF, partial_t0_mat)


nuisance_t1_IF <- mat_merge(mat_merge(mat_merge(t1_IF, t1_IF), t1_IF), t1_IF)

nuisance_numer_mat_t1 <- with(PIPW_WIPW_df, mat_merge(mat_merge(mat_merge((1 - A1) * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1),
                                                                     (1 - A1) * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)),
                                                           A1 * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)),
                                                 A1 * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)))
nuisance_numer_mat_t1 <- weighted_denom_inv %*% t(model.matrix(PIPW_WIPW_fit)) %*% (nuisance_numer_mat_t1 * residuals(PIPW_WIPW_fit))

# # nuisance_denom_mat_t1_0 <- t(model.matrix(PIPW_WIPW_fit) * with(PIPW_WIPW_df, c((1 - A1) * (1 - A0) * c(exp(q1_Z0 %*% t1)),
# #                                                                       (1 - A1) * A0 * c(exp(q1_Z0 %*% t1)),
# #                                                                       A1 * (1 - A0) * c(exp(q1_Z0 %*% t1)),
# #                                                                       A1 * A0 * c(exp(q1_Z0 %*% t1))))) %*% model.matrix(PIPW_WIPW_fit)
# # 
# # nuisance_denom_mat_t1_1 <- t(model.matrix(PIPW_WIPW_fit) * with(PIPW_WIPW_df, c((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       A1 * A0 * c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1))))) %*% model.matrix(PIPW_WIPW_fit)
# # 
# nuisance_denom_mat_t1 <- array(0, c(ncol(model.matrix(PIPW_WIPW_fit)), ncol(model.matrix(PIPW_WIPW_fit)), ncol(nuisance_t1_IF)))
# 
# temp_design_t1 <- with(PIPW_WIPW_df, cbind((1 - A1) * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1),
#                                       (1 - A1) * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1),
#                                       A1 * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1),
#                                       A1 * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)))
# 
# for (temp_rep in 1:ncol(nuisance_t1_IF)) {
#   nuisance_denom_mat_t1[, , temp_rep] <- t(model.matrix(PIPW_WIPW_fit) * c(temp_design_t1[, temp_rep])) %*% model.matrix(PIPW_WIPW_fit)
# }
# # nuisance_denom_mat_t0 <- t(model.matrix(PIPW_WIPW_fit) * with(PIPW_WIPW_df, c((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1))))) %*% model.matrix(PIPW_WIPW_fit)
# # 
# 
# temp_mat <- weighted_denom_inv %*% t(model.matrix(PIPW_WIPW_fit)) %*% (PIPW_WIPW_weights * residuals(PIPW_WIPW_fit))
# nuisance_denom_mat_t1 <- weighted_denom_inv %*% apply(nuisance_denom_mat_t1, 3, function(x) x %*% temp_mat)
# 
# 
# # nuisance_denom_mat_t1 <- t(model.matrix(PIPW_WIPW_fit) * with(PIPW_WIPW_df, c((1 - A1) * (1 - A0) * c(exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       (1 - A1) * A0 * c(exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       A1 * (1 - A0) * c(exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       A1 * A0 * c(exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1))))) %*% model.matrix(PIPW_WIPW_fit)
# 
# 
# 
nuisance_t0_IF <- mat_merge(mat_merge(mat_merge(t0_IF, t0_IF), t0_IF), t0_IF)

nuisance_numer_mat_t0 <- with(PIPW_WIPW_df, mat_merge(mat_merge(mat_merge((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z,
                                                                     (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z),
                                                           A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z),
                                                 A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z))
nuisance_numer_mat_t0 <- weighted_denom_inv %*% t(model.matrix(PIPW_WIPW_fit)) %*% (nuisance_numer_mat_t0 * residuals(PIPW_WIPW_fit))


# nuisance_denom_mat_t0 <- array(0, c(ncol(model.matrix(PIPW_WIPW_fit)), ncol(model.matrix(PIPW_WIPW_fit)), ncol(nuisance_t0_IF)))
# temp_design_t0 <- with(PIPW_WIPW_df, cbind((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z,
#                                       (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z,
#                                       A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z,
#                                       A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z))
# for (temp_rep in 1:ncol(nuisance_t0_IF)) {
#   nuisance_denom_mat_t0[, , temp_rep] <- t(model.matrix(PIPW_WIPW_fit) * c(temp_design_t0[, temp_rep])) %*% model.matrix(PIPW_WIPW_fit)
# }
# # nuisance_denom_mat_t0 <- t(model.matrix(PIPW_WIPW_fit) * with(PIPW_WIPW_df, c((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1))))) %*% model.matrix(PIPW_WIPW_fit)
# # 
# 
# temp_mat <- weighted_denom_inv %*% t(model.matrix(PIPW_WIPW_fit)) %*% (PIPW_WIPW_weights * residuals(PIPW_WIPW_fit))
# nuisance_denom_mat_t0 <- weighted_denom_inv %*% apply(nuisance_denom_mat_t0, 3, function(x) x %*% temp_mat)
# 


PIPW_WIPW_IF <- PIPW_WIPW_IF + t(nuisance_numer_mat_t1 %*% t(nuisance_t1_IF) + nuisance_numer_mat_t0 %*% t(nuisance_t0_IF))


#PIPW_WIPW_IF <- PIPW_WIPW_IF - t(nuisance_denom_mat_t1 %*% t(nuisance_t1_IF) + nuisance_denom_mat_t0 %*% t(nuisance_t0_IF))
# PIPW_WIPW_IF <- PIPW_WIPW_IF - t(weighted_denom_inv %*% nuisance_denom_mat_t1 %*% weighted_denom_inv %*% t(model.matrix(PIPW_WIPW_fit) * 
#                                                                                                    PIPW_WIPW_weights * residuals(PIPW_WIPW_fit)) %*% mat_merge(q1_Z0, mat_merge(q1_Z0, mat_merge(q1_Z0, q1_Z0))) %*% t(nuisance_t1_IF) +
#                          nuisance_denom_mat_t0 %*% weighted_denom_inv %*% t(model.matrix(PIPW_WIPW_fit) * 
#                                                                               PIPW_WIPW_weights * residuals(PIPW_WIPW_fit)) %*% mat_merge(q0_Z, mat_merge(q0_Z, mat_merge(q0_Z, q0_Z))) %*% t(nuisance_t0_IF))

PIPW_WIPW_IF <- PIPW_WIPW_IF[1:N, ] + PIPW_WIPW_IF[N + 1:N, ] + PIPW_WIPW_IF[2 * N + 1:N, ] + PIPW_WIPW_IF[3 * N + 1:N, ]

PIPW_WIPW_sd_est <- sqrt(colSums(PIPW_WIPW_IF^2))



wrong_t1 <- t1
wrong_t0 <- t0
wrong_nuisance_t1_IF <- nuisance_t1_IF
wrong_nuisance_t0_IF <- nuisance_t0_IF
wrong_q0_Y <- q0_Y
wrong_q0_Z <- q0_Z
wrong_q0_W <- q0_W
wrong_q1_Y <- q1_Y
wrong_q1_Z0 <- q1_Z0
wrong_q1_Z1 <- q1_Z1
wrong_q1_W <- q1_W

t1 <- temp_t1
t0 <- temp_t0
nuisance_t1_IF <- temp_nuisance_t1_IF
nuisance_t0_IF <- temp_nuisance_t0_IF
q0_Y <- temp_q0_Y
q0_Z <- temp_q0_Z
q0_W <- temp_q0_W
q1_Y <- temp_q1_Y
q1_Z0 <- temp_q1_Z0
q1_Z1 <- temp_q1_Z1
q1_W <- temp_q1_W

# est_result <- rbind(est_result, PIPW_WIPW_est)
# sd_result <- rbind(sd_result, PIPW_WIPW_sd_est)
# t0_result <- rbind(t0_result, t0)
# t1_result <- rbind(t1_result, t1)
# PIPW_WIPW_covering <- rbind(PIPW_covering, true_beta <= PIPW_WIPW_est + 1.96 * PIPW_WIPW_sd_est & true_beta >= PIPW_WIPW_est - 1.96 * PIPW_WIPW_sd_est)
# 
# print(rep)
# }
# 
# apply(t0_result, 2, mean)
# apply(t0_result, 2, sd)
# 
# apply(t1_result, 2, mean)
# apply(t1_result, 2, sd)
# 
# apply(est_result, 2, mean) - true_beta
# apply(est_result, 2, sd)
# apply(sd_result, 2, mean)
# apply(PIPW_WIPW_covering, 2, mean)
# 


