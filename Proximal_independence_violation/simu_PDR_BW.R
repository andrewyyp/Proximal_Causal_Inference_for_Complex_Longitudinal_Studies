PDR_BW_df <- df
PDR_BW_df$W0 <- abs(PDR_BW_df$W0)^(1/2) + 1
PDR_BW_df$W1 <- abs(PDR_BW_df$W1)^(1/2) + 1
PDR_BW_df$Z0 <- abs(PDR_BW_df$Z0)#^(1/2) + 1
PDR_BW_df$Z1 <- abs(PDR_BW_df$Z1)#^(1/2) + 1


temp_b1 <- b1
temp_b0_0 <- b0_0
temp_b0_1 <- b0_1
temp_nuisance_b1_IF <- nuisance_b1_IF
temp_nuisance_b0_IF <- nuisance_b0_IF


b1 <- wrong_b1
b0_0 <- wrong_b0_0
b0_1 <- wrong_b0_1
nuisance_b1_IF <- wrong_nuisance_b1_IF
nuisance_b0_IF <- wrong_nuisance_b0_IF


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

t1 <- wrong_t1
t0 <- wrong_t0
nuisance_t1_IF <- wrong_nuisance_t1_IF
nuisance_t0_IF <- wrong_nuisance_t0_IF
q0_Y <- wrong_q0_Y
q0_Z <- wrong_q0_Z
q0_W <- wrong_q0_W
q1_Y <- wrong_q1_Y
q1_Z0 <- wrong_q1_Z0
q1_Z1 <- wrong_q1_Z1
q1_W <- wrong_q1_W



PDR_BW_work_df <- data.frame(Y = with(PDR_BW_df, c((Y - cbind(1, 0, 0, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) +
                                               (1 - A0) * c(1 + exp(q0_Z %*% t0)) * (cbind(1, 0, 0, 0, W1, W0, X1, X0) %*% b1 - cbind(1, 0, W0, X0) %*% b0_0) +
                                               cbind(1, 0, W0, X0) %*% b0_0,
                                             (Y - cbind(1, 0, 1, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * A0 * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) +
                                               A0 * c(1 + exp(q0_Z %*% t0)) * (cbind(1, 0, 1, 0, W1, W0, X1, X0) %*% b1 - cbind(1, 1, W0, X0) %*% b0_0) +
                                               cbind(1, 1, W0, X0) %*% b0_0,
                                             (Y - cbind(1, 1, 0, 0, W1, W0, X1, X0) %*% b1) * A1 * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) +
                                               (1 - A0) * c(1 + exp(q0_Z %*% t0)) * (cbind(1, 1, 0, 0, W1, W0, X1, X0) %*% b1 - cbind(1, 0, W0, X0) %*% b0_1) +
                                               cbind(1, 0, W0, X0) %*% b0_1,
                                             (Y - cbind(1, 1, 1, 1, W1, W0, X1, X0) %*% b1) * A1 * A0 * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) +
                                               A0 * c(1 + exp(q0_Z %*% t0)) * (cbind(1, 1, 1, 1, W1, W0, X1, X0) %*% b1 - cbind(1, 1, W0, X0) %*% b0_1) +
                                               cbind(1, 1, W0, X0) %*% b0_1)),
                          A = rep(c(0, 1, 1, 2), each = N))


PDR_BW_fit <- lm(Y ~ A, data = PDR_BW_work_df)
PDR_BW_est <- PDR_BW_fit$coefficients
PDR_BW_IF <- model.matrix(PDR_BW_fit) %*% solve(crossprod(model.matrix(PDR_BW_fit))) * residuals(PDR_BW_fit)




nuisance_mat_b1 <- with(PDR_BW_df, mat_merge(mat_merge(mat_merge(cbind(1, 0, 0, 0, W1, W0, X1, X0) * (1 - A0) * (c(1 + exp(q0_Z %*% t0)) - (1 - A1) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1))),
                                                             cbind(1, 0, 1, 0, W1, W0, X1, X0) * A0 * (c(1 + exp(q0_Z %*% t0)) - (1 - A1) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)))), 
                                                   cbind(1, 1, 0, 0, W1, W0, X1, X0) * (1 - A0) * (c(1 + exp(q0_Z %*% t0)) - (1 - A1) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)))), 
                                         cbind(1, 1, 1, 1, W1, W0, X1, X0) * A0 * (c(1 + exp(q0_Z %*% t0)) - (1 - A1) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)))))
#nuisance_b1_IF <- mat_merge(mat_merge(mat_merge(b1_IF, b1_IF), b1_IF), b1_IF)

nuisance_mat_b0 <- with(PDR_BW_df,mat_merge(mat_merge(mat_merge(cbind(1, 0, W0, X0) * (1 - (1 - A0) * c(1 + exp(q0_Z %*% t0))),
                                                            cbind(1, 1, W0, X0) * (1 - A0 * c(1 + exp(q0_Z %*% t0)))), 
                                                  cbind(1, 0, W0, X0) * (1 - (1 - A0) * c(1 + exp(q0_Z %*% t0)))), 
                                        cbind(1, 1, W0, X0) * (1 - A0 * c(1 + exp(q0_Z %*% t0)))))
#nuisance_b0_IF <- mat_merge(mat_merge(mat_merge(b0_0_IF, b0_0_IF), b0_1_IF), b0_1_IF)


nuisance_mat_t0 <- with(PDR_BW_df, mat_merge(mat_merge(mat_merge(c(Y - cbind(1, 0, 0, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z + (1 - A0) * c(cbind(1, 0, 0, 0, W1, W0, X1, X0) %*% b1 - cbind(1, 0, W0, X0) %*% b0_0) * c(exp(q0_Z %*% t0)) * q0_Z,
                                                             c(Y - cbind(1, 0, 1, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z + A0 * c(cbind(1, 0, 1, 0, W1, W0, X1, X0) %*% b1 - cbind(1, 1, W0, X0) %*% b0_0) * c(exp(q0_Z %*% t0)) * q0_Z), 
                                                   c(Y - cbind(1, 1, 0, 0, W1, W0, X1, X0) %*% b1) * A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z + (1 - A0) * c(cbind(1, 1, 0, 0, W1, W0, X1, X0) %*% b1 - cbind(1, 0, W0, X0) %*% b0_0) * c(exp(q0_Z %*% t0)) * q0_Z), 
                                         c(Y - cbind(1, 1, 1, 1, W1, W0, X1, X0) %*% b1) * A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z + A0 * c(cbind(1, 1, 1, 1, W1, W0, X1, X0) %*% b1 - cbind(1, 1, W0, X0) %*% b0_0) * c(exp(q0_Z %*% t0)) * q0_Z))

#nuisance_t0_IF <- mat_merge(mat_merge(mat_merge(t0_IF, t0_IF), t0_IF), t0_IF)

nuisance_mat_t1 <- with(PDR_BW_df, mat_merge(mat_merge(mat_merge(c(Y - cbind(1, 0, 0, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) * q1_Z1),
                                                             c(Y - cbind(1, 0, 1, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) * q1_Z1)), 
                                                   c(Y - cbind(1, 1, 0, 0, W1, W0, X1, X0) %*% b1) * A1 * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) * q1_Z1)), 
                                         c(Y - cbind(1, 1, 1, 1, W1, W0, X1, X0) %*% b1) * A1 * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) * q1_Z1)))
#nuisance_t1_IF <- mat_merge(mat_merge(mat_merge(t1_IF, t1_IF), t1_IF), t1_IF)




PDR_BW_IF <- PDR_BW_IF + t(solve(crossprod(model.matrix(PDR_BW_fit))) %*% t(model.matrix(PDR_BW_fit)) %*% nuisance_mat_b1 %*% t(nuisance_b1_IF) + 
                       solve(crossprod(model.matrix(PDR_BW_fit))) %*% t(model.matrix(PDR_BW_fit)) %*% nuisance_mat_b0 %*% t(nuisance_b0_IF))

PDR_BW_IF <- PDR_BW_IF + t(solve(crossprod(model.matrix(PDR_BW_fit))) %*% t(model.matrix(PDR_BW_fit)) %*% nuisance_mat_t1 %*% t(nuisance_t1_IF) + 
                       solve(crossprod(model.matrix(PDR_BW_fit))) %*% t(model.matrix(PDR_BW_fit)) %*% nuisance_mat_t0 %*% t(nuisance_t0_IF))



PDR_BW_IF <- PDR_BW_IF[1:N, ] + PDR_BW_IF[N + 1:N, ] + PDR_BW_IF[2 * N + 1:N, ] + PDR_BW_IF[3 * N + 1:N, ]

PDR_BW_sd_est <- sqrt(colSums(PDR_BW_IF^2))

