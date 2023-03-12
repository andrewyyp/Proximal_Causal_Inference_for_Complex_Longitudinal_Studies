para_set <- list(mu_X0 = -0.35,
                 sigma_X0 = 0.5,
                 mu_U10 = 0.35,
                 sigma_U10 = 0.5,
                 mu_U20 = 0.2,
                 sigma_U10 = 0.5,
                 alpha_A0 = c(-0.5, 0.2, 0.25, -0.5),
                 mu_Z0 = c(0.2, 0.5, 0.5, 0.75, -0.8),
                 sigma_Z0 = 0.5,
                 mu_W0 = c(0.2, 0.5, -0.75, 0.75),
                 sigma_W0 = 0.5,
                 mu_X1 = c(0.2, 0.7, 0.7, 0, 0),
                 sigma_X1 = 0.5,
                 mu_U11 = c(0.2, 0.7, 0, 0.7, 0.7),
                 sigma_U11 = 0.5,
                 mu_U21 = c(0.2, 0.7, 0, 0.7, 0.7),
                 sigma_U21 = 0.5,
                 alpha_A1 = c(-0.5, 0.5, 0.15, 0.35, 0.35, 0.15, 0.35, 0.35),
                 mu_Z1 = c(0.2, 0, 0.5, 0.5, -0.75, -0.75, 0.5, 0.5, -0.75, -0.75),
                 sigma_Z1 = 0.5,
                 mu_W1 = c(0.35, 0, -0.5, -0.75, -0.75, 0.5, -0.75, -0.75),
                 sigma_W1 = 0.5,
                 mu_Y = c(-1.3, 1, 1.14, 0, 0, 0, 0.5, -0.35, -0.35, 0.2, -1, -1),
                 sigma_Y = 0.5
)

data_gen <- function(N, para_set) {
  # generate X0, U10
  X0 <- para_set$mu_X0 + rnorm(N, 0, para_set$sigma_X0)
  U10 <- para_set$mu_U10 + rnorm(N, 0, para_set$sigma_U10)
  U20 <- para_set$mu_U11 + rnorm(N, 0, para_set$sigma_U11)
  
  # generate A0
  prop_score_0 <- 1/(1 + exp(-cbind(1, X0, U10, U20) %*% para_set$alpha_A0))
  A0 <- rbinom(N, 1, prop_score_0)
  
  # generate Z0
  Z0 <- cbind(1, A0, X0, U10, U20) %*% para_set$mu_Z0 + rnorm(N, 0, para_set$sigma_Z0)
  
  # generate W0
  W0 <- cbind(1, X0, U10, U20) %*% para_set$mu_W0 + rnorm(N, 0, para_set$sigma_W0)
  
  # generate U11, X1
  X1 <- cbind(1, A0, X0, U10, U20) %*% para_set$mu_X1 + rnorm(N, 0, para_set$sigma_X1)
  U11 <- cbind(1, A0, X0, U10, U20) %*% para_set$mu_U11 + rnorm(N, 0, para_set$sigma_U11)
  U21 <- cbind(1, A0, X0, U10, U20) %*% para_set$mu_U21 + rnorm(N, 0, para_set$sigma_U21)
  
  # generate A1
  prop_score_1 <- 1/(1 + exp(-cbind(1, A0, X0, U10, U20, X1, U11, U21) %*% para_set$alpha_A1))
  A1 <- rbinom(N, 1, prop_score_1)
  
  
  # generate Z0
  Z1 <- cbind(1, Z0, A1, X1, U11, U21, A0, X0, U10, U20) %*% para_set$mu_Z1 + rnorm(N, 0, para_set$sigma_Z1)
  
  
  # generate W0
  W1 <- cbind(1, W0, X1, U11, U21, X0, U10, U20) %*% para_set$mu_W1 + rnorm(N, 0, para_set$sigma_W1)
  
  #generate Y
  Y <- cbind(1, A1, A0, A1 * A0, W1, W0, X1, U11, U21, X0, U10, U20) %*% para_set$mu_Y + rnorm(N, 0, para_set$sigma_Y)
  
  df <- data.frame(X0, A0, Z0, W0, X1, A1, Z1, W1, Y)
  return(df)
}





intervened_data_gen <- function(N, para_set, a = c(1, 1)) {
  # generate X0, U10
  X0 <- para_set$mu_X0 + rnorm(N, 0, para_set$sigma_X0)
  U10 <- para_set$mu_U10 + rnorm(N, 0, para_set$sigma_U10)
  U20 <- para_set$mu_U11 + rnorm(N, 0, para_set$sigma_U11)
  
  # generate A0
  prop_score_0 <- 1/(1 + exp(-cbind(1, X0, U10, U20) %*% para_set$alpha_A0))
  A0 <- rbinom(N, 1, prop_score_0)
  
  # generate Z0
  Z0 <- cbind(1, A0, X0, U10, U20) %*% para_set$mu_Z0 + rnorm(N, 0, para_set$sigma_Z0)
  
  # generate W0
  W0 <- cbind(1, X0, U10, U20) %*% para_set$mu_W0 + rnorm(N, 0, para_set$sigma_W0)
  
  # generate U11, X1
  X1 <- cbind(1, A0, X0, U10, U20) %*% para_set$mu_X1 + rnorm(N, 0, para_set$sigma_X1)
  U11 <- cbind(1, A0, X0, U10, U20) %*% para_set$mu_U11 + rnorm(N, 0, para_set$sigma_U11)
  U21 <- cbind(1, A0, X0, U10, U20) %*% para_set$mu_U21 + rnorm(N, 0, para_set$sigma_U21)
  
  # generate A1
  prop_score_1 <- 1/(1 + exp(-cbind(1, A0, X0, U10, U20, X1, U11, U21) %*% para_set$alpha_A1))
  A1 <- rbinom(N, 1, prop_score_1)
  
  
  # generate Z0
  Z1 <- cbind(1, Z0, A1, X1, U11, U21, A0, X0, U10, U20) %*% para_set$mu_Z1 + rnorm(N, 0, para_set$sigma_Z1)
  
  
  # generate W0
  W1 <- cbind(1, W0, X1, U11, U21, X0, U10, U20) %*% para_set$mu_W1 + rnorm(N, 0, para_set$sigma_W1)
  
  #generate Y
  Y <- cbind(1, A1, A0, A1 * A0, W1, W0, X1, U11, U21, X0, U10, U20) %*% para_set$mu_Y + rnorm(N, 0, para_set$sigma_Y)
  
  df <- data.frame(X0, A0, Z0, W0, X1, A1, Z1, W1, Y)
  return(df)
}





#U_known_est <- mean(U_known_gformula(N, df, a = c(1, 1)))# - mean(U_known_gformula(N, df, a = c(0, 0)))
#U_known_est

#U_unknown_est <- mean(U_unknown_gformula(N, df, a = c(1, 1))) #- mean(U_unknown_gformula(N, df, a = c(0, 0)))
#U_unknown_est

#h0link_para <- proximal_gformula(N, df)
#proximal_est <- mean(cbind(1, rep(1, N), rep(1, N), df$W0, df$X0) %*% h0link_para) - mean(cbind(1, rep(0, N), rep(0, N), df$W0, df$X0) %*% h0link_para)
#proximal_est


#proximal_est <- mean(cbind(1, rep(1, N), rep(1, N), df$W0, df$X0) %*% proximal_gformula(N, df, c(1, 1))) - 
#  mean(cbind(1, rep(0, N), rep(0, N), df$W0, df$X0) %*% proximal_gformula(N, df, c(0, 0)))
#proximal_est