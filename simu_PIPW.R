PIPW_df <- df
inioptim_t0 <- c(0, 0, 0, 0)
q0_Y <- with(PIPW_df, matrix(rep(c(0, 1, 0, 0), each = length(A0)), nrow = length(A0), ncol = 4))
q0_Z <- with(PIPW_df, (-1)^(1 - A0) * cbind(1, A0, Z0, X0))
q0_W <- with(PIPW_df, (-1)^(1 - A0) * cbind(1, A0, W0, X0))

q0link <- optim(par = inioptim_t0, fn = MSE_func,
                bridge_func = q0bridge, Y = q0_Y, W = q0_W, Z = q0_Z,
                method = "BFGS", hessian = FALSE, control = list(maxit = round))
t0 <- q0link$par


q1_Y <- with(PIPW_df, matrix(rep(c(0, 1, 0, 0, 0, 0, 0), each = length(A0)), nrow = length(A0), ncol = 7))
q1_Z0 <- with(PIPW_df, (-1)^(1 - A1) * cbind(1, A1, A0, Z1, Z0, X1, X0))
q1_Z1 <- with(PIPW_df, (-1)^(1 - A1) * cbind(1, A1, A0, Z1, Z0, X1, X0))

q1_W <- with(PIPW_df, (-1)^(1 - A1) * cbind(1, A1, A0, W1, W0, X1, X0))


inioptim_t1 <- c(0, 0, 0, 0, 0, 0, 0)
q1link <- optim(par = inioptim_t1, fn = MSE_func_q1,
                bridge_func = q1bridge, Y = q1_Y, W = q1_W, Z0 = q1_Z0, Z1 = q1_Z1, q0 = c(exp(q0_Z %*% t0)), 
                method = "BFGS", hessian = FALSE, control = list(maxit = round))
t1 <- q1link$par
PIPW_converging <- q1link$convergence == 0


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


t0_IF <- q0_IF(t0, Y = q0_Y, W = q0_W, Z = q0_Z)

t1_IF <- q1_IF(t1, Y = q1_Y, W = q1_W, Z0 = q1_Z0, Z1 = q1_Z1, q0 = c(exp(q0_Z %*% t0)), t0_IF, partial_t0_mat)




nuisance_t1_IF <- mat_merge(mat_merge(mat_merge(t1_IF, t1_IF), t1_IF), t1_IF)

nuisance_numer_mat_t1 <- with(PIPW_df, mat_merge(mat_merge(mat_merge((1 - A1) * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1),
                                                                     (1 - A1) * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)),
                                                           A1 * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)),
                                                 A1 * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0)) * c(exp(q1_Z1 %*% t1)) * q1_Z1)))
nuisance_numer_mat_t1 <- weighted_denom_inv %*% t(model.matrix(PIPW_fit)) %*% (nuisance_numer_mat_t1 * residuals(PIPW_fit))


nuisance_t0_IF <- mat_merge(mat_merge(mat_merge(t0_IF, t0_IF), t0_IF), t0_IF)

nuisance_numer_mat_t0 <- with(PIPW_df, mat_merge(mat_merge(mat_merge((1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z,
                                                                     (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z),
                                                           A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z),
                                                 A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z))
nuisance_numer_mat_t0 <- weighted_denom_inv %*% t(model.matrix(PIPW_fit)) %*% (nuisance_numer_mat_t0 * residuals(PIPW_fit))



PIPW_IF <- PIPW_IF + t(nuisance_numer_mat_t1 %*% t(nuisance_t1_IF) + nuisance_numer_mat_t0 %*% t(nuisance_t0_IF))



PIPW_IF <- PIPW_IF[1:N, ] + PIPW_IF[N + 1:N, ] + PIPW_IF[2 * N + 1:N, ] + PIPW_IF[3 * N + 1:N, ]

PIPW_sd_est <- sqrt(colSums(PIPW_IF^2))




