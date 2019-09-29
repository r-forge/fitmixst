require(mvtnorm)
require(cubature)
require(tmvtnorm)
require(mnormt)

p <- 3
sigma <- diag(rep(1,p))
nu <- 5
a <- rep(-1000,p)
b <- rep(0,p)
mu  <- rep(1,p)

adaptIntegrate(function(x) x%*%t(x)*dmvt(x,delta=mu,sigma=sigma,df=nu,log=F),
               lowerLimit=a,upperLimit=b,fDim=p)$integral/
                 (pmvt(lower=a,upper=b,delta=mu,sigma=sigma,df=nu,type="shifted")[1])

adaptIntegrate(function(x) x%*%t(x)*dmt(x,mean=mu,S=sigma,df=nu,log=F),lowerLimit=a,upperLimit=b,fDim=p)$integral/
  (sadmvt(lower=a,upper=b,mean=mu,S=sigma,df=nu))


adaptIntegrate(function(x) x*dmvt(x,delta=c(0,0),sigma=sigma[-3,-3],df=nu,log=F),
               lowerLimit=a[-3],upperLimit=b[-3],fDim=p)$integral/
                 (pmvt(lower=a[-3],upper=b[-3],delta=c(0,0),sigma=sigma[-3,-3],df=nu,type="shifted")[1])

adaptIntegrate(function(x) x*dtmvt(x,mean=mu,sigma=sigma,df=nu,log=F,lower=a,upper=b),
               lowerLimit=a,upperLimit=b,fDim=p)



truncm1 <- function (nu, a, b, sigma) {
  p <- length(a)
  qstar <- qstarb <- alpha <- alphab <- delta <- deltab <- vector(length=p)
  
  lambda1 <- (nu-2)/nu;
  lambda2 <- (nu-4)/nu;
  
  alphat <- pmvt(lower=a,upper=b,delta=rep(0,p),sigma=sigma,df=nu);
  
  for(i in 1:p)
  {    
    delta[i] <- (nu-1)/(nu+a[i]^2);
    deltab[i] <- (nu-1)/(nu+b[i]^2);
    
    alpha[i] <- biv.nt.prob(lower=a[-i],upper=b[-i], mean=rep(0,p-1), S=(sigma[-i,-i]/delta[i]),df=nu-1) 
    alpha[i] <-  pmvt(lower=a[-i],upper=b[-i], delta=rep(0,p-1), sigma=(sigma[-i,-i]/delta[i]),df=nu-1)
    
    alpha[i] <- pmvt(upper=b[-i],delta=rep(0,p-1),sigma=sigma[-i,-i]/delta[i],df=nu-1,type="shifted")-
      pmvt(upper=a[-i],delta=rep(0,p-1),sigma=sigma[-i,-i]/delta[i],df=nu-1,type="shifted");
    
    sadmvt(lower=rep(-Inf,p-1),upper=b[-i],mean=rep(0,p-1),S=sigma[-i,-i]/delta[i],df=nu-1)-
      pmt(a[-i],mean=rep(0,p-1),S=sigma[-i,-i]/delta[i],df=nu-1);
    
    #  alphab[i] <- pmt(b[-i],mean=rep(0,p-1),S=sigma[-i,-i]/deltab[i],df=nu-1)-
    #    pmt(a[-i],mean=rep(0,p-1),S=sigma[-i,-i]/deltab[i],df=nu-1);  
    
    #alphab[i] <- sadmvt(lower=a[-i],upper=b[-i], mean=rep(0,p-1), S=(sigma[-i,-i]/deltab[i]),df=nu-1)
    alphab[i] <-  pmvt(lower=a[-i],upper=b[-i], delta=rep(0,p-1), sigma=(sigma[-i,-i]/deltab[i]),df=nu-1)
    #alphat <- biv.nt.prob(lower=c(a,a),upper=c(b,b),mean=c(0,0),S=sigma,df=nu)[1]#-pmt(a,mean=c(0,0),S=sigma,df=nu)
    
    qstar[i] <- lambda1^0.5 * dt(a[i]*lambda1^0.5,df=nu-2) * alpha[i];
    qstarb[i] <- lambda1^0.5 * dt(b[i]*lambda1^0.5,df=nu-2) * alphab[i];
  }
  ew <- (alphat*lambda1)^(-1)* sigma %*% (qstar-qstarb);
  ew
}

truncm1(5,a[-3],b[-3],sigma[-3,-3])
truncm1(5,a,b,sigma)

truncm1(5,rep(-1000,3),a,sigma)




truncm2 <- function (a, b, sigma, nu) {
  
  p <- length(a);
  d <- qstar <- qstarb <- alpha <- alphab <- delta <- deltab <- vector(length=p);
  H <- alphaij <- deltaij <- matrix(0,p,p);
  
  lambda1 <- (nu-2)/nu;
  lambda2 <- (nu-4)/nu;
  
  alphat <- pmvt(lower=a,upper=b,delta=rep(0,p),sigma=sigma,df=nu);
  alphastar <- pmvt(lower=a,upper=b,delta=rep(0,p),sigma=sigma/lambda1,df=nu-2,type="shifted");
  
  for(i in 1:p)
  {    
    delta[i] <- (nu-1)/(nu+a[i]^2);
    deltab[i] <- (nu-1)/(nu+b[i]^2);
    
    alpha[i] <-  pmvt(lower=a[-i],upper=b[-i], delta=rep(0,p-1), sigma=(sigma[-i,-i]/delta[i]),df=nu-1)  
    alphab[i] <-  pmvt(lower=a[-i],upper=b[-i], delta=rep(0,p-1), sigma=(sigma[-i,-i]/deltab[i]),df=nu-1)
    
    
    qstar[i] <- lambda1^0.5 * dt(a[i]*lambda1^0.5,df=nu-2) * alpha[i];
    qstarb[i] <- lambda1^0.5 * dt(b[i]*lambda1^0.5,df=nu-2) * alphab[i];
    
    for(j in 1:p){
      
      deltaij <- function(a,b,s){
        (nu-2)/(nu+(1-s^2)^(-1)*(a^2-2*s*a*b+b^2));
      }
      
      
      hijnew <- function(lower,upper,a,b,sigma,nu,i,j){
        1/lambda2 * dmvt(c(a[i],b[j]),delta=c(0,0),sigma=sigma[c(i,j),c(i,j)]/lambda2,df=nu-4,log=F)*
          alphaijnew(lower[c(-i,-j)],upper[c(-i,-j)],a,b,sigma,nu,i,j);
      }
      
      alphaijnew <- function(lower,upper,a,b,sigma,nu,i,j){
        res <- 1
        if(length(a)>2){
          
          #anpassen in mvtnorm... achtung p==1
          res <- sadmvt(lower=lower,upper=upper,mean=rep(0,p-2),
                        S= sigma[c(-i,-j),c(-i,-j)]/deltaij(a[i],b[j],sigma[i,j]),df=nu-2);
          
        }
        res
      }
      
      if(i!=j){
        H[i,j] <- hijnew(a,b,a,a,sigma,nu,i,j)-hijnew(a,b,a,b,sigma,nu,i,j)-
          hijnew(a,b,b,a,sigma,nu,i,j)+hijnew(a,b,b,b,sigma,nu,i,j);
      }
      
      
    }
    
    d[i] <- a[i]*qstar[i]-b[i]*qstarb[i]-(sigma*H)[i,i]
    
  }
  
  D <- diag(d)
  
  ew2 <- (alphat*lambda1)^(-1) *(alphastar*sigma + sigma%*%(H+D)%*%sigma)
  ew2
}

truncm2(c(0,0),c(5,5),sigma[-3,-3],5)

truncm2(a,b,sigma,5)



#allgemeiner fall fehlt:
#benutze cov2corr für corr matrix
# sigma = sqrt(diag(sigma))


#o'hagan
#instert mu is missing
test1 <- function(a,mu,sigma,nu){
  p <- length(a)
  m1 <- vector(length=p)
  for (i in 1:p){
    astar  <-  (a[-i]-mu[-i])-(a[i]-mu[i])*1/sigma[i,i]*sigma[,i][-i];
    sigmastar <- (nu+1/sigma[i,i]*(a[i]-mu[i])^2)/(nu-1)*(sigma[-i,-i]-1/sigma[i,i]*
      sigma[,i][-i]%*%t(sigma[,i][-i]))
    m1[i] <-(2*pi*sigma[i,i])^(-0.5) * (nu/(nu+1/sigma[i,i]*(a[i]-mu[i])^2))^((nu-1)/2)*
      gamma((nu-1)/2)/gamma(nu/2)*
      sqrt(nu/2)*pmt(astar,mean=rep(0,p-1),S=sigmastar,nu-1)/pmt((a-mu),mean=rep(0,p),S=sigma,df=nu)
  }
  m2 <- mu-sigma%*%m1
  m2
}

#insert mu
test1(b,mu,sigma,nu)

test2 <- function(a,sigma,nu){
  p <- length(a)
  h <- matrix(0,p,p)
  m1 <- vector(length=p)
  tmp <- 1
  tmp2 <- vector(length=p)
  
  for(i in 1:p){
    #mu
    astar  <-  a[-i]-a[-i]*1/sigma[i,i]*sigma[,i][-i];
    sigmastar <- (nu+1/sigma[i,i]*a[i]^2/(nu-1))*(sigma[-i,-i]-1/sigma[i,i]*
      sigma[,i][-i]%*%t(sigma[,i][-i]))
    m1[i] <-(2*pi*sigma[i,i])^(-0.5) * (nu/(nu+1/sigma[i,i]*a[i]^2))^((nu-1)/2)*
      gamma((nu-1)/2)/gamma(nu/2)*
      sqrt(nu/2)*pmt(astar,mean=rep(0,p-1),S=sigmastar,nu-1)/pmt(a,mean=rep(0,p),S=sigma,df=nu)
    
    for(j in 1:p){
      
      if(i!=j){
        
        #sigma
        nustar  <- nu+c(t(a[c(i,j)])%*%solve(sigma[c(i,j),c(i,j)])%*%(a[c(i,j)]))
        astarstar <- (a[c(-i,-j)]) - cbind(sigma[,i][c(-i,-j)],sigma[,j][c(-i,-j)])%*%
          solve(sigma[c(i,j),c(i,j)])%*%a[c(i,j)]        
        if(p>2){
          sigmastarstar <- nustar/(nu-2)*(sigma[c(-i,-j),c(-i,-j)]-cbind(sigma[,i][c(-i,-j)],sigma[,j][c(-i,-j)])%*%
            solve(sigma[c(i,j),c(i,j)])%*%t(cbind(sigma[,i][c(-i,-j)],sigma[,j][c(-i,-j)])))
          tmp <- pmt(astarstar,mean=rep(0,p-2),sigmastarstar,nu-2)
        }
        
        
        h[i,j] <- -1/(2*pi*sqrt(sigma[i,i]*sigma[j,j]-sigma[i,j]^2))*(nu/(nu-2))*(nu/nustar)^(nu/2-1)*tmp/
          pmt(a,mean=rep(0,p),S=sigma,df=nu)
        tmp2[i] <- tmp2[i] + sigma[i,j]*h[i,j]
      }      
    }
  }
  
  for(i in 1:p){
    h[i,i] <- 1/sigma[i,i]*a[i]*m1[i]-1/sigma[i,i]*tmp2[i]
  }
  
  nu/(nu-2)*pmt(a,mean=rep(0,p),S=nu/(nu-2)*sigma,df=nu-2)/pmt(a,mean=rep(0,p),S=sigma,df=nu)*
    sigma -sigma%*%h%*%sigma
  
}

test2(b,sigma,nu)



truncatedt <- function(a,mu,sigma,nu){
  p <- length(a)
  m1 <- vector(length=p)
  h  <- m3 <- matrix(0,p,p)
  tmp <- 1
  tmp2 <- vector(length=p)
  
  for (i in 1:p){
    #evtl falsch in ohagan
    astar  <-  (a[-i]-mu[-i])-(a[i]-mu[i])*1/sigma[i,i]*sigma[,i][-i];
    sigmastar <- (nu+1/sigma[i,i]*(a[i]-mu[i])^2)/(nu-1)*(sigma[-i,-i]-1/sigma[i,i]*
      sigma[,i][-i]%*%t(sigma[,i][-i]))
    m1[i] <-(2*pi*sigma[i,i])^(-0.5) * (nu/(nu+1/sigma[i,i]*(a[i]-mu[i])^2))^((nu-1)/2)*
      gamma((nu-1)/2)/gamma(nu/2)*
      sqrt(nu/2)*pmt(astar,mean=rep(0,p-1),S=sigmastar,nu-1)/pmt((a-mu),mean=rep(0,p),S=sigma,df=nu)
    
    for(j in 1:p){
      if(i!=j){
        
        #sigma
        nustar  <- nu+c(t(a[c(i,j)]-mu[c(i,j)])%*%solve(sigma[c(i,j),c(i,j)])%*%(a[c(i,j)]-mu[c(i,j)]))
        astarstar <- (a[c(-i,-j)]-mu[c(-i,-j)]) - cbind(sigma[,i][c(-i,-j)],sigma[,j][c(-i,-j)])%*%
          solve(sigma[c(i,j),c(i,j)])%*%(a[c(i,j)]-mu[c(i,j)])        
        if(p>2){
          sigmastarstar <- nustar/(nu-2)*(sigma[c(-i,-j),c(-i,-j)]-cbind(sigma[,i][c(-i,-j)],sigma[,j][c(-i,-j)])%*%
            solve(sigma[c(i,j),c(i,j)])%*%t(cbind(sigma[,i][c(-i,-j)],sigma[,j][c(-i,-j)])))
          tmp <- pmt(astarstar,mean=rep(0,p-2),sigmastarstar,nu-2)
        }
        
        
        h[i,j] <- -1/(2*pi*sqrt(sigma[i,i]*sigma[j,j]-sigma[i,j]^2))*(nu/(nu-2))*(nu/nustar)^(nu/2-1)*tmp/
          pmt(a-mu,mean=rep(0,p),S=sigma,df=nu)
        tmp2[i] <- tmp2[i] + sigma[i,j]*h[i,j]
      }  
    }
  }
  
  for(i in 1:p){
    h[i,i] <- 1/sigma[i,i]*(a[i]-mu[i])*m1[i]-1/sigma[i,i]*tmp2[i]
  }
  
  m2 <- mu-sigma%*%m1
  m3  <- -mu %*% t(mu) + mu %*% t(m2) + m2 %*% t(mu) + nu/(nu-2)*pmt(a-mu,mean=rep(0,p),S=nu/(nu-2)*
    sigma,df=nu-2)/pmt(a-mu,mean=rep(0,p),S=sigma,df=nu)*sigma -sigma%*%h%*%sigma
  list(m2,m3)
}

truncatedt(b,mu,sigma,nu)



