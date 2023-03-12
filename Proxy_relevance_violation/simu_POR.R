
POR_df <- df
h1_Y <- with(POR_df, Y)
h1_W <- with(POR_df, cbind(1, A1, A0, A1 * A0, W1, W0, X1, X0))
h1_Z <- with(POR_df, cbind(1, A1, A0, A1 * A0, Z1, Z0, X1, X0))
inioptim_b1 <- c(0, 0, 0, 0, 0, 0, 0, 0)
h1link <- optim(par = inioptim_b1, fn = MSE_func,
                bridge_func = h1bridge, Y = h1_Y, W = h1_W, Z = h1_Z, 
                method = "BFGS", hessian = FALSE)
b1 <- h1link$par
h1_a1_1 <- with(POR_df, cbind(1, 1, A0, 1 * A0, W1, W0, X1, X0) %*% b1)
h1_a1_0 <- with(POR_df, cbind(1, 0, A0, 0 * A0, W1, W0, X1, X0) %*% b1)


h0_W <- with(POR_df, cbind(1, A0, W0, X0))
h0_Z <- with(POR_df, cbind(1, A0, Z0, X0))

inioptim_b0 <- c(0, 0, 0, 0)
h0_a1_1 <- optim(par = inioptim_b0, fn = MSE_func,
                    bridge_func = h0bridge, Y = h1_a1_1, W = h0_W, Z = h0_Z, 
                    method = "BFGS", hessian = FALSE)
b0_1 <- h0_a1_1$par
h0_a1_0 <- optim(par = inioptim_b0, fn = MSE_func,
                    bridge_func = h0bridge, Y = h1_a1_0, W = h0_W, Z = h0_Z, 
                    method = "BFGS", hessian = FALSE)
b0_0 <- h0_a1_0$par


POR_work_df <- data.frame(Y = with(POR_df, c(cbind(1, 0, W0, X0) %*% b0_0, cbind(1, 1, W0, X0) %*% b0_0, 
                                                 cbind(1, 0, W0, X0) %*% b0_1, cbind(1, 1, W0, X0) %*% b0_1)),
                          A = rep(c(0, 1, 1, 2), each = N))



POR_fit <- lm(Y ~ A, data = POR_work_df)
POR_IF <- model.matrix(POR_fit) %*% solve(crossprod(model.matrix(POR_fit))) * residuals(POR_fit)


POR_est <- POR_fit$coefficients

b1_IF <- h1_IF(b1, Y = h1_Y, W = h1_W, Z = h1_Z)

partial_b1_1_mat <- partial_b1(h0_Z,  with(POR_df, cbind(1, 1, A0, 1 * A0, W1, W0, X1, X0)))


b0_1_IF <- h0_IF(b0_1, Y = h1_a1_1, W = h0_W, Z = h0_Z, b1_IF, partial_b1_1_mat)


partial_b1_0_mat <- partial_b1(h0_Z, with(POR_df, cbind(1, 0, A0, 0 * A0, W1, W0, X1, X0)))

b0_0_IF <- h0_IF(b0_0, Y = h1_a1_0, W = h0_W, Z = h0_Z, b1_IF, partial_b1_0_mat)


nuisance_mat <- with(POR_df, mat_merge(mat_merge(mat_merge(cbind(1, 0, W0, X0), cbind(1, 1, W0, X0)), cbind(1, 0, W0, X0)), cbind(1, 1, W0, X0)))
nuisance_IF <- mat_merge(mat_merge(mat_merge(b0_0_IF, b0_0_IF), b0_1_IF), b0_1_IF)

POR_IF <- POR_IF + t(solve(crossprod(model.matrix(POR_fit))) %*% t(model.matrix(POR_fit)) %*% nuisance_mat %*% t(nuisance_IF))


POR_IF <- POR_IF[1:N, ] + POR_IF[N + 1:N, ] + POR_IF[2 * N + 1:N, ] + POR_IF[3 * N + 1:N, ]
POR_sd_est <- sqrt(colSums(POR_IF^2)) 
