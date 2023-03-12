df <- data_gen(N, para_set)

PIPW_df <- df
inioptim_t0 <- c(0, 0, 0, 0)
q0_Y <- with(PIPW_df, matrix(rep(c(0, 1, 0, 0), each = length(A0)), nrow = length(A0), ncol = 4))
q0_Z <- with(PIPW_df, (-1)^(1 - A0) * cbind(1, A0, Z0, X0))
q0_W <- with(PIPW_df, (-1)^(1 - A0) * cbind(1, A0, W0, X0))

q0link <- optim(par = inioptim_t0, fn = MSE_func,
                bridge_func = q0bridge, Y = q0_Y, W = q0_W, Z = q0_Z,
                method = "BFGS", hessian = FALSE)
t0 <- q0link$par

#solve(t(q0_W) %*% (c(exp(q0_Z %*% t0)) * q0_Z)) / q0link$hessian

q1_Y <- with(PIPW_df, matrix(rep(c(0, 1, 0, 0, 0, 0, 0), each = length(A0)), nrow = length(A0), ncol = 7))
q1_Z0 <- with(PIPW_df, (-1)^(1 - A1) * cbind(1, A1, A0, Z1, Z0, X1, X0))
q1_Z1 <- with(PIPW_df, (-1)^(1 - A1) * cbind(1, A1, A0, Z1, Z0, X1, X0))

q1_W <- with(PIPW_df, (-1)^(1 - A1) * cbind(1, A1, A0, W1, W0, X1, X0))


inioptim_t1 <- c(0, 0, 0, 0, 0, 0, 0)
q1link <- optim(par = inioptim_t1, fn = MSE_func_q1,
                bridge_func = q1bridge, Y = q1_Y, W = q1_W, Z0 = q1_Z0, Z1 = q1_Z1, q0 = c(exp(q0_Z %*% t0)), 
                method = "BFGS", hessian = FALSE)
t1 <- q1link$par


# 
# qbridge <- function(tet, x) {
#   para <- tet
#   t0 <- para[1:4]
#   t1 <- para[-(1:4)]
#   A0 <- x[, 1]
#   W0 <- x[, 2]
#   Z0 <- x[, 3]
#   X0 <- x[, 4]
#   A1 <- x[, 5]
#   W1 <- x[, 6]
#   Z1 <- x[, 7]
#   X1 <- x[, 8]
#   #x = A0, W0, Z0, X0, A1, W1, Z1, X1
# 
#   q0_Y <- matrix(rep(c(0, 1, 0, 0), each = length(A0)), nrow = length(A0), ncol = 4)
#   q0_Z <- (-1)^(1 - A0) * cbind(1, A0, Z0, X0)
#   q0_W <- (-1)^(1 - A0) * cbind(1, A0, W0, X0)
#   t0link <- c(1 + exp(q0_Z %*% t0))
# 
#   q1_Z0 <-  (-1)^(1 - A1) * cbind(1, A1, A0, Z1, Z0, X1, X0, 0)
#   q1_Z1 <- (-1)^(1 - A1) * cbind(1, A1, A0, Z1, Z0, X1, X0, (-1)^(1 - A0))
# 
#   q1_W <- (-1)^(1 - A1) * cbind(1, A1, A0, W1, W0, X1, X0, A0 * W0^2, A0 * W1^2, A1 * W0^2, A1 * W1^2)
#   q1_Y <- cbind(matrix(rep(c(0, 1, 0, 0, 0, 0, 0, 0, 0), each = length(A0)), nrow = length(A0), ncol = ncol(q1_W) - 2), W0^2, W1^2)
# 
# 
#   t1link <- c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1))
# 
#   f <- cbind(t0link * q0_W - q0_Y, t1link * q1_W - t0link * q1_Y)
# }
# 
# fit <- gmm(qbridge, x = with(PIPW_df, cbind(A0, W0, Z0, X0, A1, W1, Z1, X1)), rep(0, 4 + 8))
# 

PIPW_work_df <- data.frame(Y = with(PIPW_df, rep(Y, 4)),
                           A = rep(c(0, 1, 1, 2), each = N))

PIPW_weights <- with(PIPW_df, c((1 - A1) * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)), 
                                (1 - A1) * A0 * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)),
                                A1 * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)), 
                                A1 * A0 * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1))))

PIPW_fit <- lm(Y ~ A, weights = PIPW_weights, data = PIPW_work_df)
weighted_denom_inv <- solve(t(model.matrix(PIPW_fit) * PIPW_weights) %*% model.matrix(PIPW_fit))
PIPW_IF <- (PIPW_weights * model.matrix(PIPW_fit)) %*% weighted_denom_inv * residuals(PIPW_fit)

PIPW_est <- PIPW_fit$coefficients


partial_t0_mat <- partial_t0(t1 = t1, q1_Y = q1_Y, q1_W = q1_W, q1_Z0 = q1_Z0, q1_Z1 = q1_Z1, q0_Z = q0_Z, q0 = c(exp(q0_Z %*% t0)))
# 
# data3 <- list()
# data3$q0_Y <- q0_Y
# data3$q0_Z <- q0_Z
# data3$q0_W <- q0_W
# 
# data3$q1_Y <- q1_Y
# data3$q1_Z0 <- q1_Z0
# data3$q1_Z1 <- q1_Z1
# 
# data3$q1_W <- q1_W
# 
# data3$t1 <- t1
# 
# partial_t0_mat <- gradient_appro(partial_t0_numer, t0, data3)
# 


t0_IF <- q0_IF(t0, Y = q0_Y, W = q0_W, Z = q0_Z)

t1_IF <- q1_IF(t1, Y = q1_Y, W = q1_W, Z0 = q1_Z0, Z1 = q1_Z1, q0 = c(exp(q0_Z %*% t0)), t0_IF, partial_t0_mat)



# sqrt(colSums(t0_IF^2))
# partial_t0_mat
# sqrt(colSums(t1_IF^2))
# solve(t(q1_W) %*% (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1))



nuisance_t1_IF <- mat_merge(mat_merge(mat_merge(t1_IF, t1_IF), t1_IF), t1_IF)

nuisance_numer_mat_t1 <- with(PIPW_df, mat_merge(mat_merge(mat_merge((1 - A1) * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1),
                                                                     (1 - A1) * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)),
                                                           A1 * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)),
                                                 A1 * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)))
nuisance_numer_mat_t1 <- weighted_denom_inv %*% t(model.matrix(PIPW_fit)) %*% (nuisance_numer_mat_t1 * residuals(PIPW_fit))

# # nuisance_denom_mat_t1_0 <- t(model.matrix(PIPW_fit) * with(PIPW_df, c((1 - A1) * (1 - A0) * c(exp(q1_Z0 %*% t1)),
# #                                                                       (1 - A1) * A0 * c(exp(q1_Z0 %*% t1)),
# #                                                                       A1 * (1 - A0) * c(exp(q1_Z0 %*% t1)),
# #                                                                       A1 * A0 * c(exp(q1_Z0 %*% t1))))) %*% model.matrix(PIPW_fit)
# # 
# # nuisance_denom_mat_t1_1 <- t(model.matrix(PIPW_fit) * with(PIPW_df, c((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       A1 * A0 * c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1))))) %*% model.matrix(PIPW_fit)
# # 
# nuisance_denom_mat_t1 <- array(0, c(ncol(model.matrix(PIPW_fit)), ncol(model.matrix(PIPW_fit)), ncol(nuisance_t1_IF)))
# 
# temp_design_t1 <- with(PIPW_df, mat_merge(mat_merge(mat_merge((1 - A1) * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1),
#                                       (1 - A1) * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)),
#                                       A1 * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)),
#                                       A1 * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)))
# 
# for (temp_rep in 1:ncol(nuisance_t1_IF)) {
#   nuisance_denom_mat_t1[, , temp_rep] <- t(model.matrix(PIPW_fit) * c(temp_design_t1[, temp_rep])) %*% model.matrix(PIPW_fit)
# }
# # nuisance_denom_mat_t0 <- t(model.matrix(PIPW_fit) * with(PIPW_df, c((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1))))) %*% model.matrix(PIPW_fit)
# # 
# 
# temp_mat <- weighted_denom_inv %*% t(model.matrix(PIPW_fit)) %*% (PIPW_weights * residuals(PIPW_fit))
# nuisance_denom_mat_t1 <- weighted_denom_inv %*% apply(nuisance_denom_mat_t1, 3, function(x) x %*% temp_mat)
# 
# 
# # nuisance_denom_mat_t1 <- t(model.matrix(PIPW_fit) * with(PIPW_df, c((1 - A1) * (1 - A0) * c(exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       (1 - A1) * A0 * c(exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       A1 * (1 - A0) * c(exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)),
# #                                                                       A1 * A0 * c(exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1))))) %*% model.matrix(PIPW_fit)
# 
# 
# 
nuisance_t0_IF <- mat_merge(mat_merge(mat_merge(t0_IF, t0_IF), t0_IF), t0_IF)

nuisance_numer_mat_t0 <- with(PIPW_df, mat_merge(mat_merge(mat_merge((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z,
                                                                     (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z),
                                                           A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z),
                                                 A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z))
nuisance_numer_mat_t0 <- weighted_denom_inv %*% t(model.matrix(PIPW_fit)) %*% (nuisance_numer_mat_t0 * residuals(PIPW_fit))

# 
# nuisance_denom_mat_t0 <- array(0, c(ncol(model.matrix(PIPW_fit)), ncol(model.matrix(PIPW_fit)), ncol(nuisance_t0_IF)))
# temp_design_t0 <- with(PIPW_df, mat_merge(mat_merge(mat_merge((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z,
#                                   (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z),
#                                   A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z),
#                                   A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z))
# for (temp_rep in 1:ncol(nuisance_t0_IF)) {
#   nuisance_denom_mat_t0[, , temp_rep] <- t(model.matrix(PIPW_fit) * c(temp_design_t0[, temp_rep])) %*% model.matrix(PIPW_fit)
# }
# # nuisance_denom_mat_t0 <- t(model.matrix(PIPW_fit) * with(PIPW_df, c((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)),
# #                                                                     A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1))))) %*% model.matrix(PIPW_fit)
# # 
# 
# temp_mat <- weighted_denom_inv %*% t(model.matrix(PIPW_fit)) %*% (PIPW_weights * residuals(PIPW_fit))
# nuisance_denom_mat_t0 <- weighted_denom_inv %*% apply(nuisance_denom_mat_t0, 3, function(x) x %*% temp_mat)
# 


PIPW_IF <- PIPW_IF + t(nuisance_numer_mat_t1 %*% t(nuisance_t1_IF) + nuisance_numer_mat_t0 %*% t(nuisance_t0_IF))


#PIPW_IF <- PIPW_IF - t(nuisance_denom_mat_t1 %*% t(nuisance_t1_IF) + nuisance_denom_mat_t0 %*% t(nuisance_t0_IF))
# PIPW_IF <- PIPW_IF - t(weighted_denom_inv %*% nuisance_denom_mat_t1 %*% weighted_denom_inv %*% t(model.matrix(PIPW_fit) * 
#                                                                                                    PIPW_weights * residuals(PIPW_fit)) %*% mat_merge(q1_Z0, mat_merge(q1_Z0, mat_merge(q1_Z0, q1_Z0))) %*% t(nuisance_t1_IF) +
#                          nuisance_denom_mat_t0 %*% weighted_denom_inv %*% t(model.matrix(PIPW_fit) * 
#                                                                               PIPW_weights * residuals(PIPW_fit)) %*% mat_merge(q0_Z, mat_merge(q0_Z, mat_merge(q0_Z, q0_Z))) %*% t(nuisance_t0_IF))

PIPW_IF <- PIPW_IF[1:N, ] + PIPW_IF[N + 1:N, ] + PIPW_IF[2 * N + 1:N, ] + PIPW_IF[3 * N + 1:N, ]

PIPW_sd_est <- sqrt(colSums(PIPW_IF^2))
#PIPW_sd_est 




