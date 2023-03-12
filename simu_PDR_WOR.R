PDR_WOR_df <- df
PDR_WOR_df$W0 <- transform_W(PDR_WOR_df$W0)
PDR_WOR_df$W1 <- transform_W(PDR_WOR_df$W1)

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


PDR_WOR_work_df <- data.frame(Y = with(PDR_WOR_df, c((Y - cbind(1, 0, 0, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * (1 - A0) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) +
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


PDR_WOR_fit <- lm(Y ~ A, data = PDR_WOR_work_df)
PDR_WOR_est <- PDR_WOR_fit$coefficients
PDR_WOR_IF <- model.matrix(PDR_WOR_fit) %*% solve(crossprod(model.matrix(PDR_WOR_fit))) * residuals(PDR_WOR_fit)




nuisance_mat_b1 <- with(PDR_WOR_df, mat_merge(mat_merge(mat_merge(cbind(1, 0, 0, 0, W1, W0, X1, X0) * (1 - A0) * (c(1 + exp(q0_Z %*% t0)) - (1 - A1) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1))),
                                                             cbind(1, 0, 1, 0, W1, W0, X1, X0) * A0 * (c(1 + exp(q0_Z %*% t0)) - (1 - A1) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)))), 
                                                   cbind(1, 1, 0, 0, W1, W0, X1, X0) * (1 - A0) * (c(1 + exp(q0_Z %*% t0)) - (1 - A1) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)))), 
                                         cbind(1, 1, 1, 1, W1, W0, X1, X0) * A0 * (c(1 + exp(q0_Z %*% t0)) - (1 - A1) * c(1 + exp(q0_Z %*% t0) + exp(q1_Z0 %*% t1) + exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)))))

nuisance_mat_b0 <- with(PDR_WOR_df,mat_merge(mat_merge(mat_merge(cbind(1, 0, W0, X0) * (1 - (1 - A0) * c(1 + exp(q0_Z %*% t0))),
                                                            cbind(1, 1, W0, X0) * (1 - A0 * c(1 + exp(q0_Z %*% t0)))), 
                                                  cbind(1, 0, W0, X0) * (1 - (1 - A0) * c(1 + exp(q0_Z %*% t0)))), 
                                        cbind(1, 1, W0, X0) * (1 - A0 * c(1 + exp(q0_Z %*% t0)))))


nuisance_mat_t0 <- with(PDR_WOR_df, mat_merge(mat_merge(mat_merge(c(Y - cbind(1, 0, 0, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z + (1 - A0) * c(cbind(1, 0, 0, 0, W1, W0, X1, X0) %*% b1 - cbind(1, 0, W0, X0) %*% b0_0) * c(exp(q0_Z %*% t0)) * q0_Z,
                                                             c(Y - cbind(1, 0, 1, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z + A0 * c(cbind(1, 0, 1, 0, W1, W0, X1, X0) %*% b1 - cbind(1, 1, W0, X0) %*% b0_0) * c(exp(q0_Z %*% t0)) * q0_Z), 
                                                   c(Y - cbind(1, 1, 0, 0, W1, W0, X1, X0) %*% b1) * A1 * (1 - A0) * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z + (1 - A0) * c(cbind(1, 1, 0, 0, W1, W0, X1, X0) %*% b1 - cbind(1, 0, W0, X0) %*% b0_1) * c(exp(q0_Z %*% t0)) * q0_Z), 
                                         c(Y - cbind(1, 1, 1, 1, W1, W0, X1, X0) %*% b1) * A1 * A0 * c(exp(q0_Z %*% t0)) * c(1 + exp(q1_Z1 %*% t1)) * q0_Z + A0 * c(cbind(1, 1, 1, 1, W1, W0, X1, X0) %*% b1 - cbind(1, 1, W0, X0) %*% b0_1) * c(exp(q0_Z %*% t0)) * q0_Z))


nuisance_mat_t1 <- with(PDR_WOR_df, mat_merge(mat_merge(mat_merge(c(Y - cbind(1, 0, 0, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) * q1_Z1),
                                                             c(Y - cbind(1, 0, 1, 0, W1, W0, X1, X0) %*% b1) * (1 - A1) * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) * q1_Z1)), 
                                                   c(Y - cbind(1, 1, 0, 0, W1, W0, X1, X0) %*% b1) * A1 * (1 - A0) * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) * q1_Z1)), 
                                         c(Y - cbind(1, 1, 1, 1, W1, W0, X1, X0) %*% b1) * A1 * A0 * (c(exp(q1_Z0 %*% t1)) * q1_Z0 + c(exp(q0_Z %*% t0) * exp(q1_Z1 %*% t1)) * q1_Z1)))




PDR_WOR_IF <- PDR_WOR_IF + t(solve(crossprod(model.matrix(PDR_WOR_fit))) %*% t(model.matrix(PDR_WOR_fit)) %*% nuisance_mat_b1 %*% t(nuisance_b1_IF) + 
                       solve(crossprod(model.matrix(PDR_WOR_fit))) %*% t(model.matrix(PDR_WOR_fit)) %*% nuisance_mat_b0 %*% t(nuisance_b0_IF))

PDR_WOR_IF <- PDR_WOR_IF + t(solve(crossprod(model.matrix(PDR_WOR_fit))) %*% t(model.matrix(PDR_WOR_fit)) %*% nuisance_mat_t1 %*% t(nuisance_t1_IF) + 
                       solve(crossprod(model.matrix(PDR_WOR_fit))) %*% t(model.matrix(PDR_WOR_fit)) %*% nuisance_mat_t0 %*% t(nuisance_t0_IF))



PDR_WOR_IF <- PDR_WOR_IF[1:N, ] + PDR_WOR_IF[N + 1:N, ] + PDR_WOR_IF[2 * N + 1:N, ] + PDR_WOR_IF[3 * N + 1:N, ]

PDR_WOR_sd_est <- sqrt(colSums(PDR_WOR_IF^2))


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


