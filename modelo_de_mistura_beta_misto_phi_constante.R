
params <- c("sigma_b","beta","phi","pi_ranked")
  
data <- list("Y","X","Z","N","n","p","d","K"
             ,"zerovec","alpha"
             ,"psi_priori","C"
             ,"a_phi","b_phi"
             ,"mu_beta","diagonal_matrix_alpha_beta")

cat("model
    {
    
    for( i in 1 : N ) {
    S[i] ~ dcat(pi[1:K])
    
    for(j in 1:n) {
    Y[i , j] ~ dbeta(a1[i,j] ,a2[i,j])
    a1[i,j] <- mu[i , j]*phi
    a2[i,j] <- (1-mu[i , j])*phi
    logit(mu[i , j]) <-  inprod(X[i, j, 1:p], beta[S[i],]) + inprod(Z[i,j,1:d],b[i])
    }
    
    b[i] ~ dnorm(zerovec,psi)
    }
    
    pi[1:K] ~ ddirch(alpha[1:K]) 
    psi ~ dgamma(psi_priori, C)
    
    sigma_b <- 1/pow(psi,0.5)
    
    phi ~ dunif(a_phi,b_phi)
    
    for(s in 1:K) {
    beta[s,1:p] ~ dmnorm(mu_beta[,s], diagonal_matrix_alpha_beta[,])
    pi_ranked[(K+1-s)] <- ranked(pi[],s)
    }
    }", fill=TRUE, file=paste0("beta_mixed_mixture_linear_d_1_K_phi_1.txt"))

fit_bugs <- bugs(data,inits,params,paste0("beta_mixed_mixture_linear_d_1_K_phi_1.txt")
                 ,n.thin=nt, n.chains=nc,n.burnin=nb, n.iter=ni
                 ,debug=TRUE, DIC=TRUE, clearWD=TRUE)
