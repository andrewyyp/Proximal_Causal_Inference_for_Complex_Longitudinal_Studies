MSE_func <- function(bridge_func, para, Y, W, Z){
  g0 <- bridge_func(para = para, Y = Y, W = W, Z = Z)
  g <- apply(g0, 2, mean)
  gmmf <- sum(g^2)
  return(gmmf)
}

h1bridge <- function(para, Y, W, Z) {
  h1link <- W %*% para
  g1 <- Z
  g <- c(Y - h1link) * g1
  return(g)
}

h0bridge <- function(para, Y, W, Z) {
 
  h0link <- W %*% para
  g0 <- Z
  g <- c(Y - h0link) * g0
  return(g)
}





q0bridge <- function(para, Y, W, Z) { 
  t0link <- 1 + exp(Z %*% para)
  g0 <- W
  g <- c(t0link) * g0 - Y
  return(g)
}


q1bridge <- function(para, Y, W, Z0, Z1, q0) {
  t1link <- 1 + q0 + exp(Z0 %*% para) + q0 * exp(Z1 %*% para)
  g1 <- W
  g <- c(t1link) * g1 - (1 + q0) * Y
  return(g)
}

MSE_func_q1 <- function(bridge_func, para, Y, W, Z0, Z1, q0){
  g0 <- bridge_func(para = para, Y = Y, W = W, Z0 = Z0, Z1 = Z1, q0 = q0)
  g <- apply(g0, 2, mean)
  gmmf <- sum(g^2)
  return(gmmf)
}

partial_b1 <- function(Z0, W1) {
  return(t(Z0) %*% W1)
}

h0_IF <- function(para, Y, W, Z, b1_IF, partial_b1_mat) {
  return(t(solve(t(Z) %*% W) %*% t(Z * c(Y - c(W %*% para)) + t(partial_b1_mat %*% t(b1_IF)))))
}

h1_IF <- function(para, Y, W, Z) {
  return(t(solve(t(Z) %*% W) %*% t(Z * (Y - c(W %*% para)))))
}


partial_t0 <- function(t1, q1_Y, q1_W, q1_Z0, q1_Z1, q0_Z, q0) {
  denom_inv <- solve(t(q1_W) %*% (c(exp(q1_Z0 %*% t1)) * q1_Z0 + q0 * c(exp(q1_Z1 %*% t1)) * q1_Z1))
  numer_mat <- denom_inv %*% t((q1_W * (1 + c(exp(q1_Z1 %*% t1))) - q1_Y) * q0) %*% q0_Z
 
  denom_mat <- 0
  return(numer_mat - denom_mat)
}


partial_t0_numer <- function(para, data){
  q0_Y <- data$q0_Y
  q0_Z <- data$q0_Z
  q0_W <- data$q0_W
  
  q0 <- c(exp(q0_Z %*% para))
  q1_Y <- data$q1_Y
  q1_Z0 <- data$q1_Z0
  q1_Z1 <- data$q1_Z1
  
  q1_W <- data$q1_W

  t1 <- data$t1
  denom_inv <- solve(t(q1_W) %*% (c(exp(q1_Z0 %*% t1)) * q1_Z0 + q0 * c(exp(q1_Z1 %*% t1)) * q1_Z1))
  score <- t(denom_inv %*% t(q1_W * (1 + q0 + c(exp(q1_Z0 %*% t1)) + q0 * c(exp(q1_Z1 %*% t1))) - q1_Y * (1 + q0)))
  
  return(score)
}



q0_IF <- function(para, Y, W, Z) {
  return(-t(solve(t(W) %*% (c(exp(Z %*% para)) * Z)) %*% t(c(1 + exp(Z %*% para)) * W - Y)))
}


q1_IF <- function(para, Y, W, Z0, Z1, q0, t0_IF, partial_t0_mat) {
  return(-t(solve(t(W) %*% (c(exp(Z0 %*% para)) * Z0 + q0 * c(exp(Z1 %*% para)) * Z1)) %*% t(c(1 + q0 + exp(Z0 %*% para) + q0 * exp(Z1 %*% para)) * W - 
                          Y * (1 + q0))) + t(partial_t0_mat %*% t(t0_IF)))
}

trimmed_mean <- function(x, alpha = 0.05) {
  return(mean(x[x <= quantile(x, 1 - alpha/2) & x >= quantile(x, alpha/2)]))
}

trimmed_sd <- function(x) {
  x <- x[x <=  quantile(x, 0.975) & x >= quantile(x, 0.025)]
  return(sd(x))
}


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



score_func <- function(bfun, para, data) {
  score <- apply(bfun(para, data), 2, sum)
  return(score)
}

gradient_appro <- function(bfun, para, data){
  gradient <- jacobian(func = score_func, bfun = bfun, x = para, data = data)
  return(gradient)
}

ipw.var <- function(para, data){
  PIPW_fit <- data1$PIPW_fit
  PIPW_df <- data1$PIPW_df
  
  q0_Y <- data1$q0_Y
  q0_Z <- data1$q0_Z
  q0_W <- data1$q0_W
  
  
  
  
  q1_Y <- data1$q1_Y
  q1_Z0 <- data1$q1_Z0
  q1_Z1 <- data1$q1_Z1
  
  q1_W <- data1$q1_W
  t0 <- para[1:4]
  t1 <- para[-(1:4)]
  PIPW_weights <- with(PIPW_df, c((1 - A1) * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)), 
                                  (1 - A1) * A0 * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)),
                                  A1 * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)), 
                                  A1 * A0 * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1))))
  
  weighted_denom_inv <- data1$weighted_denom_inv
  
  g <- t(weighted_denom_inv %*% (t(model.matrix(PIPW_fit)) * PIPW_weights * rep(PIPW_df$Y, 4))) * length(PIPW_weights)
  return(g)
}



PIPW_var <- function(para, data){
  PIPW_fit <- data$PIPW_fit
  q0_Y <- data$q0_Y
  q0_Z <- data$q0_Z
  q0_W <- data$q0_W
  
  
  
  
  q1_Y <- data$q1_Y
  q1_Z0 <- data$q1_Z0
  q1_Z1 <- data$q1_Z1
  
  q1_W <- data$q1_W
  t0_1 <- para[1:4]
  t1_1 <- para[5:11]
  t0_2 <- para[11 + 1:4]
  t1_2 <- para[11 + 5:11]
  t0_3 <- para[22 + 1:4]
  t1_3 <- para[22 + 5:11]
  t0_4 <- para[33 + 1:4]
  t1_4 <- para[33 + 5:11]
  PIPW_weights <- with(PIPW_df, c((1 - A1) * (1 - A0) * c(1 + exp(q0_Z %*% t0_1) + exp(q1_Z0 %*% t1_1) + exp(q0_Z %*% t0_1) * exp(q1_Z1 %*% t1_1)), 
                                  (1 - A1) * A0 * c(1 + exp(q0_Z %*% t0_2) + exp(q1_Z0 %*% t1_2) + exp(q0_Z %*% t0_2) * exp(q1_Z1 %*% t1_2)),
                                  A1 * (1 - A0) * c(1 + exp(q0_Z %*% t0_3) + exp(q1_Z0 %*% t1_3) + exp(q0_Z %*% t0_3) * exp(q1_Z1 %*% t1_3)), 
                                  A1 * A0 * c(1 + exp(q0_Z %*% t0_4) + exp(q1_Z0 %*% t1_4) + exp(q0_Z %*% t0_4) * exp(q1_Z1 %*% t1_4))))
  
  weighted_denom_inv <- solve(t(model.matrix(PIPW_fit) * PIPW_weights) %*% model.matrix(PIPW_fit))
  
  score <- t(weighted_denom_inv %*% (t(model.matrix(PIPW_fit)) * PIPW_weights * residuals(PIPW_fit))) * length(PIPW_weights)
  return(score)
}


t1_var <- function(para, data){
  PIPW_fit <- data$PIPW_fit
  q0_Y <- data$q0_Y
  q0_Z <- data$q0_Z
  q0_W <- data$q0_W
  q1_Y <- data$q1_Y
  q1_Z0 <- data$q1_Z0
  q1_Z1 <- data$q1_Z1
  q1_W <- data$q1_W
  t0 <- para[1:4]
  t1 <- para[5:11]
  N <- nrow(q1_Y)
  PIPW_weights <- with(PIPW_df, c((1 - A1) * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)), 
                                  (1 - A1) * A0 * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)),
                                  A1 * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)), 
                                  A1 * A0 * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1))))
  
  weighted_denom_inv <- solve(t(model.matrix(PIPW_fit)) %*% model.matrix(PIPW_fit))
  
  score <- t(weighted_denom_inv %*% t(model.matrix(PIPW_fit))) * PIPW_weights * rep(PIPW_df$Y, 4) 
  score <- score[1:N, ] + score[N + 1:N, ] + score[2 * N + 1:N, ] + score[3 * N + 1:N, ]
  return(score)
}



PIPW_var3 <- function(para, data){
  PIPW_fit <- data$PIPW_fit
  q0_Y <- data$q0_Y
  q0_Z <- data$q0_Z
  q0_W <- data$q0_W
  
  
  
  
  q1_Y <- data$q1_Y
  q1_Z0 <- data$q1_Z0
  q1_Z1 <- data$q1_Z1
  
  q1_W <- data$q1_W
  t0_1 <- para[1:4]
  t1_1 <- para[5:11]
  t0_2 <- para[11 + 1:4]
  t1_2 <- para[11 + 5:11]
  t0_3 <- para[22 + 1:4]
  t1_3 <- para[22 + 5:11]
  t0_4 <- para[33 + 1:4]
  t1_4 <- para[33 + 5:11]
  PIPW_weights <- with(PIPW_df, c((1 - A1) * (1 - A0) * c(1 + exp(q0_Z %*% t0_1) + exp(q1_Z0 %*% t1_1) + exp(q0_Z %*% t0_1) * exp(q1_Z1 %*% t1_1)), 
                                  (1 - A1) * A0 * c(1 + exp(q0_Z %*% t0_2) + exp(q1_Z0 %*% t1_2) + exp(q0_Z %*% t0_2) * exp(q1_Z1 %*% t1_2)),
                                  A1 * (1 - A0) * c(1 + exp(q0_Z %*% t0_3) + exp(q1_Z0 %*% t1_3) + exp(q0_Z %*% t0_3) * exp(q1_Z1 %*% t1_3)), 
                                  A1 * A0 * c(1 + exp(q0_Z %*% t0_4) + exp(q1_Z0 %*% t1_4) + exp(q0_Z %*% t0_4) * exp(q1_Z1 %*% t1_4))))
  
  weighted_denom_inv <- solve(t(model.matrix(PIPW_fit) * PIPW_weights) %*% model.matrix(PIPW_fit))
  
  score <- t(weighted_denom_inv %*% (t(model.matrix(PIPW_fit)) * PIPW_weights * residuals(PIPW_fit)))
  return(score)
}
