# simulation for testing purposes

inv_logit <- function(a){
  exp(a)/(1 + exp(a))
}

nof1.ordinal.simulation <- function(baseline = "baseline", Base.size = 100, Treat.size = 100,
                                    alpha = 0, beta_A = -0.1, beta_B = -0.3, cut = c(0.5,1,1.5,2), ncat = 5){
  
  Treat <- rep(baseline, Base.size)
  Treat <- c(Treat, rep(c("B", "A", "B", "A"), each = Treat.size))
  nobs <- length(Treat)
  
  Y <- epsilon <- mu <- rep(NA, nobs)
  
  q <- matrix(0, nrow = length(Treat), ncol = ncat - 1)
  p <- matrix(0, nrow = nobs, ncol = ncat)
  
  for(i in 1:nobs){
    
    mu[i] <- alpha
    if(Treat[i] == "A"){
      mu[i] <- mu[i] + beta_A
    } else if (Treat[i] == "B"){
      mu[i] <- mu[i] + beta_B
    }
    
    for(r in 1:(ncat-1)){
      q[i,r] <- inv_logit(mu[i] - cut[r])
    }
    
    p[i,1] <- 1 - q[i,1]
    for(r in 2:(ncat-1)){
      p[i,r] <- q[i,r-1] - q[i,r]
    }
    p[i,ncat] <- q[i,(ncat-1)]
    Y[i] <- sample(1:ncat, size = 1,  prob = p[i,])
  }  
  list(Y = Y, Treat = Treat)
}


nof1.binomial.simulation <- function(baseline = "baseline", Base.size = 14, Treat.size = 56,
                                     alpha = 0.5, beta_A = -0.1, beta_B = -0.05){
  
  Treat <- rep(baseline, Base.size)
  Treat <- c(Treat, rep(c("B", "A", "B", "A"), each = Treat.size))
  nobs <- length(Treat)
  
  Y <- epsilon <- mu <- p <- rep(NA, nobs)
    
  for(i in 1:nobs){
  
    mu[i] <- alpha
    if(Treat[i] == "A"){
      mu[i] <- mu[i] + beta_A
    } else if (Treat[i] == "B"){
      mu[i] <- mu[i] + beta_B
    }
    Y[i] <- rbinom(1,1, inv_logit(mu[i]))
  }
  
  Time <- 1:length(Treat)
  nobs <- length(Time)
  
  list(Y = Y, Treat = Treat)
}


  
