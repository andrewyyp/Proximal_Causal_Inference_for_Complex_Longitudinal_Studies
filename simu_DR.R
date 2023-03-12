DR_df <- df
ipw_stage0 <- glm(A0 ~ X0 + W0 + Z0, family = "binomial", data = DR_df)
q0_Z <- with(DR_df, cbind(1, X0, W0, Z0))
ipw_stage1 <- glm(A1 ~ A0 + X1 + W1 + Z1 + X0 + W0 + Z0, family = "binomial", data = DR_df)
q1_Z <- with(DR_df, cbind(1, A0, X1, W1, Z1, X0, W0, Z0))

p_of_A0 <- fitted(ipw_stage0)
p_of_A1 <- fitted(ipw_stage1)
t0_IF <- iid(ipw_stage0)
t1_IF <- iid(ipw_stage1)

standard_or_df <- df
or_stage1 <- lm(Y ~ A1 + A0 + A1 * A0 + X1 + W1 + Z1 + X0 + W0 + Z0, data = standard_or_df)
b1 <- or_stage1$coefficients
b1_IF <- iid(or_stage1)

new_or_df <- standard_or_df
new_or_df$A1 <- 1
or_stage0_a11 <- lm(predict(or_stage1, newdata = new_or_df) ~ A0 + X0 + W0 + Z0, data = standard_or_df)
b0_1 <- or_stage0_a11$coefficients
b0_1_IF <- iid(or_stage0_a11) + t(solve(crossprod(model.matrix(or_stage0_a11))) %*% t(model.matrix(or_stage0_a11)) %*% model.matrix(or_stage1) %*% t(b1_IF))

new_or_df$A1 <- 0
or_stage0_a10 <- lm(predict(or_stage1, newdata = new_or_df) ~ A0 + X0 + W0 + Z0, data = standard_or_df)
b0_0 <- or_stage0_a10$coefficients
b0_0_IF <- iid(or_stage0_a10) + t(solve(crossprod(model.matrix(or_stage0_a11))) %*% t(model.matrix(or_stage0_a10)) %*% model.matrix(or_stage1) %*% t(b1_IF))

DR_work_df <- data.frame(Y = with(DR_df, c((Y - cbind(1, 0, 0, X1, W1, Z1, X0, W0, Z0, 0) %*% b1) * (1 - A1) * (1 - A0) / (1 - p_of_A1) / (1 - p_of_A0) +
                                          (cbind(1, 0, 0, X1, W1, Z1, X0, W0, Z0, 0) %*% b1 - cbind(1, 0, X0, W0, Z0) %*% b0_0) * (1 - A0) / (1 - p_of_A0) +
                                          cbind(1, 0, X0, W0, Z0) %*% b0_0,
                                          (Y - cbind(1, 0, 1, X1, W1, Z1, X0, W0, Z0, 0) %*% b1) * (1 - A1) * A0 / (1 - p_of_A1) / p_of_A0 +
                                            (cbind(1, 0, 1, X1, W1, Z1, X0, W0, Z0, 0) %*% b1 - cbind(1, 1, X0, W0, Z0) %*% b0_0) * A0 / p_of_A0 +
                                            cbind(1, 1, X0, W0, Z0) %*% b0_0,
                                          (Y - cbind(1, 1, 0, X1, W1, Z1, X0, W0, Z0, 0) %*% b1) * A1 * (1 - A0) / p_of_A1 / (1 - p_of_A0) +
                                            (cbind(1, 1, 0,  X1, W1, Z1, X0, W0, Z0, 0) %*% b1 - cbind(1, 0, X0, W0, Z0) %*% b0_1) * (1 - A0) / (1 - p_of_A0) +
                                            cbind(1, 0, X0, W0, Z0) %*% b0_0,
                                          (Y - cbind(1, 1, 1,  X1, W1, Z1, X0, W0, Z0, 1) %*% b1) * A1 * A0 / p_of_A1 / p_of_A0 +
                                            (cbind(1, 1, 1, X1, W1, Z1, X0, W0, Z0, 1) %*% b1 - cbind(1, 1, X0, W0, Z0) %*% b0_1) * A0 / p_of_A0 +
                                            cbind(1, 1, X0, W0, Z0) %*% b0_0)),
                                      A = rep(c(0, 1, 1, 2), each = N))


DR_fit <- lm(Y ~ A, data = DR_work_df)
DR_est <- DR_fit$coefficients

DR_IF <- model.matrix(DR_fit) %*% solve(crossprod(model.matrix(DR_fit))) * residuals(DR_fit)




nuisance_mat_b1 <- with(DR_df, mat_merge(mat_merge(mat_merge(cbind(1, 0, 0, X1, W1, Z1, X0, W0, Z0, 0) * (1 - A0) / (1 - p_of_A0) * (1 - (1 - A1) / (1 - p_of_A1)),
                                                            cbind(1, 0, 1, X1, W1, Z1, X0, W0, Z0, 0) * A0 / p_of_A0 * (1 - (1 - A1) / (1 - p_of_A1))), 
                                                            cbind(1, 1, 0, X1, W1, Z1, X0, W0, Z0, 0) * (1 - A0) / (1 - p_of_A0) * (1 - A1 / p_of_A1)), 
                                                            cbind(1, 1, 1, X1, W1, Z1, X0, W0, Z0, 1) * A0 / p_of_A0 * (1 - A1 / p_of_A1)))

nuisance_b1_IF <- mat_merge(mat_merge(mat_merge(b1_IF, b1_IF), b1_IF), b1_IF)

nuisance_mat_b0 <- with(DR_df, mat_merge(mat_merge(mat_merge(cbind(1, 0, X0, W0, Z0) * (1 - (1 - A0) / (1 - p_of_A0)),
                                                            cbind(1, 1, X0, W0, Z0) * (1 - A0 / p_of_A0)), 
                                                  cbind(1, 0, X0, W0, Z0) * (1 - (1 - A0) / (1 - p_of_A0))), 
                                                  cbind(1, 1, X0, W0, Z0) * (1 - A0 / p_of_A0)))
nuisance_b0_IF <- mat_merge(mat_merge(mat_merge(b0_0_IF, b0_0_IF), b0_1_IF), b0_1_IF)


nuisance_mat_t0 <- with(DR_df, mat_merge(mat_merge(mat_merge(c(Y - cbind(1, 0, 0, X1, W1, Z1, X0, W0, Z0, 0) %*% b1) * (1 - A1) * (1 - A0) / (1 - p_of_A1) / (1 - p_of_A0) / p_of_A0 * -q0_Z + 
                                                               (1 - A0) * c(cbind(1, 0, 0, X1, W1, Z1, X0, W0, Z0, 0) %*% b1 - cbind(1, 0, X0, W0, Z0) %*% b0_0) / (1 - p_of_A0) / p_of_A0 * -q0_Z,
                                                             c(Y - cbind(1, 0, 1, X1, W1, Z1, X0, W0, Z0, 0) %*% b1) * (1 - A1) * A0 / (1 - p_of_A1) / p_of_A0 / (1 - p_of_A0) * q0_Z + 
                                                               A0 * c(cbind(1, 0, 1, X1, W1, Z1, X0, W0, Z0, 0) %*% b1 - cbind(1, 1, X0, W0, Z0) %*% b0_0) / p_of_A0 / (1 - p_of_A0) * q0_Z), 
                                                   c(Y - cbind(1, 1, 0, X1, W1, Z1, X0, W0, Z0, 0) %*% b1) * A1 * (1 - A0) / p_of_A1 / (1 - p_of_A0) / p_of_A0 * -q0_Z + 
                                                     (1 - A0) * c(cbind(1, 1, 0, X1, W1, Z1, X0, W0, Z0, 0) %*% b1 - cbind(1, 0, X0, W0, Z0) %*% b0_1) / (1 - p_of_A0) / p_of_A0 * -q0_Z), 
                                         c(Y - cbind(1, 1, 1, X1, W1, Z1, X0, W0, Z0, 1) %*% b1) * A1 * A0 / p_of_A1 / p_of_A0 / (1 - p_of_A0) * q0_Z + 
                          A0 * c(cbind(1, 1, 1, X1, W1, Z1, X0, W0, Z0, 1) %*% b1 - cbind(1, 1, X0, W0, Z0) %*% b0_1) / p_of_A0 / (1 - p_of_A0) * q0_Z))

nuisance_t0_IF <- mat_merge(mat_merge(mat_merge(t0_IF, t0_IF), t0_IF), t0_IF)

nuisance_mat_t1 <- with(DR_df, mat_merge(mat_merge(mat_merge(c(Y - cbind(1, 0, 0, X1, W1, Z1, X0, W0, Z0, 0) %*% b1) * (1 - A1) * (1 - A0) / (1 - p_of_A1) / p_of_A1 / (1 - p_of_A0) * -q1_Z,
                                                            c(Y - cbind(1, 0, 1, X1, W1, Z1, X0, W0, Z0, 0) %*% b1) * (1 - A1) * A0 / (1 - p_of_A1) / p_of_A1 / p_of_A0 * -q1_Z), 
                                                  c(Y - cbind(1, 1, 0, X1, W1, Z1, X0, W0, Z0, 0) %*% b1) * A1 * (1 - A0) / p_of_A1 / (1 - p_of_A1) / (1 - p_of_A0) * q1_Z), 
                                        c(Y - cbind(1, 1, 1, X1, W1, Z1, X0, W0, Z0, 1) %*% b1) * A1 * A0 / p_of_A1 / (1 - p_of_A1) /p_of_A0 * q1_Z))

nuisance_t1_IF <- mat_merge(mat_merge(mat_merge(t1_IF, t1_IF), t1_IF), t1_IF)




DR_IF <- DR_IF + t(solve(crossprod(model.matrix(DR_fit))) %*% t(model.matrix(DR_fit)) %*% nuisance_mat_b1 %*% t(nuisance_b1_IF) + 
                       solve(crossprod(model.matrix(DR_fit))) %*% t(model.matrix(DR_fit)) %*% nuisance_mat_b0 %*% t(nuisance_b0_IF))

DR_IF <- DR_IF + t(solve(crossprod(model.matrix(DR_fit))) %*% t(model.matrix(DR_fit)) %*% nuisance_mat_t1 %*% t(nuisance_t1_IF) + 
                       solve(crossprod(model.matrix(DR_fit))) %*% t(model.matrix(DR_fit)) %*% nuisance_mat_t0 %*% t(nuisance_t0_IF))



DR_IF <- DR_IF[1:N, ] + DR_IF[N + 1:N, ] + DR_IF[2 * N + 1:N, ] + DR_IF[3 * N + 1:N, ]

DR_sd_est <- sqrt(colSums(DR_IF^2))
