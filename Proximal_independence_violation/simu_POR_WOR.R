POR_WOR_df <- df
POR_WOR_df$W0 <- abs(POR_WOR_df$W0)^(1/2) + 1
POR_WOR_df$W1 <- abs(POR_WOR_df$W1)^(1/2) + 1

temp_b1 <- b1
temp_b0_0 <- b0_0
temp_b0_1 <- b0_1
temp_nuisance_b1_IF <- nuisance_b1_IF
temp_nuisance_b0_IF <- nuisance_b0_IF


h1_Y <- with(POR_WOR_df, Y)
h1_W <- with(POR_WOR_df, cbind(1, A1, A0, A1 * A0, W1, W0, X1, X0))
h1_Z <- with(POR_WOR_df, cbind(1, A1, A0, A1 * A0, Z1, Z0, X1, X0))
inioptim_b1 <- c(0, 0, 0, 0, 0, 0, 0, 0)
h1link <- optim(par = inioptim_b1, fn = MSE_func,
                bridge_func = h1bridge, Y = h1_Y, W = h1_W, Z = h1_Z, 
                method = "BFGS", hessian = FALSE)
b1 <- h1link$par
h1_a1_1 <- with(POR_WOR_df, cbind(1, 1, A0, 1 * A0, W1, W0, X1, X0) %*% b1)
h1_a1_0 <- with(POR_WOR_df, cbind(1, 0, A0, 0 * A0, W1, W0, X1, X0) %*% b1)


h0_W <- with(POR_WOR_df, cbind(1, A0, W0, X0))
h0_Z <- with(POR_WOR_df, cbind(1, A0, Z0, X0))

inioptim_b0 <- c(0, 0, 0, 0)
h0_a1_1 <- optim(par = inioptim_b0, fn = MSE_func,
                 bridge_func = h0bridge, Y = h1_a1_1, W = h0_W, Z = h0_Z, 
                 method = "BFGS", hessian = FALSE)
b0_1 <- h0_a1_1$par
h0_a1_0 <- optim(par = inioptim_b0, fn = MSE_func,
                 bridge_func = h0bridge, Y = h1_a1_0, W = h0_W, Z = h0_Z, 
                 method = "BFGS", hessian = FALSE)
b0_0 <- h0_a1_0$par




POR_WOR_work_df <- data.frame(Y = with(POR_WOR_df, c(cbind(1, 0, W0, X0) %*% b0_0, cbind(1, 1, W0, X0) %*% b0_0, 
                                             cbind(1, 0, W0, X0) %*% b0_1, cbind(1, 1, W0, X0) %*% b0_1)),
                          A = rep(c(0, 1, 1, 2), each = N))



POR_WOR_fit <- lm(Y ~ A, data = POR_WOR_work_df)
POR_WOR_IF <- model.matrix(POR_WOR_fit) %*% solve(crossprod(model.matrix(POR_WOR_fit))) * residuals(POR_WOR_fit)


POR_WOR_est <- POR_WOR_fit$coefficients

b1_IF <- h1_IF(b1, Y = h1_Y, W = h1_W, Z = h1_Z)

partial_b1_1_mat <- partial_b1(h0_Z,  with(POR_WOR_df, cbind(1, 1, A0, 1 * A0, W1, W0, X1, X0)))


b0_1_IF <- h0_IF(b0_1, Y = h1_a1_1, W = h0_W, Z = h0_Z, b1_IF, partial_b1_1_mat)


partial_b1_0_mat <- partial_b1(h0_Z, with(POR_WOR_df, cbind(1, 0, A0, 0 * A0, W1, W0, X1, X0)))

b0_0_IF <- h0_IF(b0_0, Y = h1_a1_0, W = h0_W, Z = h0_Z, b1_IF, partial_b1_0_mat)

mat_merge <- function(a, b){
  a_r <- nrow(a)
  a_c <- ncol(a)
  b_r <- nrow(b)
  b_c <- ncol(b)
  result <- matrix(0, nrow = a_r + b_r, ncol = a_c + b_c)
  
  result[1:a_r, 1:a_c] = a
  result[a_r + 1:b_r, a_c + 1:b_c] = b
  return(result)
}
nuisance_mat <- with(POR_WOR_df,mat_merge(mat_merge(mat_merge(cbind(1, 0, W0, X0), cbind(1, 1, W0, X0)), cbind(1, 0, W0, X0)), cbind(1, 1, W0, X0)))
nuisance_IF <- mat_merge(mat_merge(mat_merge(b0_0_IF, b0_0_IF), b0_1_IF), b0_1_IF)

POR_WOR_IF <- POR_WOR_IF + t(solve(crossprod(model.matrix(POR_WOR_fit))) %*% t(model.matrix(POR_WOR_fit)) %*% nuisance_mat %*% t(nuisance_IF))


POR_WOR_IF <- POR_WOR_IF[1:N, ] + POR_WOR_IF[N + 1:N, ] + POR_WOR_IF[2 * N + 1:N, ] + POR_WOR_IF[3 * N + 1:N, ]
POR_WOR_sd_est <- sqrt(colSums(POR_WOR_IF^2)) 

wrong_b1 <- b1
wrong_b0_0 <- b0_0
wrong_b0_1 <- b0_1
wrong_nuisance_b1_IF <- nuisance_b1_IF
wrong_nuisance_b0_IF <- nuisance_b0_IF

b1 <- temp_b1
b0_0 <- temp_b0_0
b0_1 <- temp_b0_1
nuisance_b1_IF <- temp_nuisance_b1_IF
nuisance_b0_IF <- temp_nuisance_b0_IF

