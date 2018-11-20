
params <- c("sigma_b","beta","sigma_u","delta","pi_ranked")

data <- list("Y","X","Z","W","R","N","n","p","d","q","d_star","K"
             ,"zerovec","zerovec_phi","alpha"
             ,"psi_priori","C","psi_phi_priori","C_phi"
             ,"mu_delta","diagonal_matrix_alpha_delta"
             ,"mu_beta","diagonal_matrix_alpha_beta")

cat("model
    {
    
    for( i in 1 : N ) {
    S[i] ~ dcat(pi[1:K])
    
    for(j in 1:n) {
    Y[i,j] ~ dbeta(a1[i,j] ,a2[i,j])
    a1[i,j] <- mu[i,j]*phi[i,j]
    a2[i,j] <- (1-mu[i,j])*phi[i,j]
    logit(mu[i,j]) <- inprod(X[i,j,1:p], beta[S[i],]) + inprod(Z[i,j,1:d],b[i])
    log(phi[i,j]) <- inprod(W[i,j,1:q],delta[S[i],]) + inprod(R[i,j,1:d_star],u[i])
  }
    b[i] ~ dnorm(zerovec,psi)
    u[i] ~ dnorm(zerovec_phi,psi_phi)
    }
    
    pi[1:K] ~ ddirch(alpha[1:K]) 
    psi ~ dgamma(psi_priori, C)
    psi_phi ~ dgamma(psi_phi_priori, C_phi)
    
    sigma_b <- 1/pow(psi,0.5)
    sigma_u <- 1/pow(psi_phi,0.5)
    
    for(s in 1:K) {
    beta[s,1:p]  ~ dmnorm(mu_beta[,s], diagonal_matrix_alpha_beta[,])
    delta[s,1:q] ~ dmnorm(mu_delta[,s], diagonal_matrix_alpha_delta[,])
    pi_ranked[(K+1-s)] <- ranked(pi[],s)
    }
    }", fill=TRUE, file=paste0("beta_mixed_mixture_linear.txt"))

fit_bugs <- bugs(data,inits,params,paste0("beta_mixed_mixture_linear.txt")
                 ,n.thin=nt, n.chains=nc,n.burnin=nb, n.iter=ni
                 ,debug=TRUE, DIC=TRUE, clearWD=TRUE)