context("Fitting mixtures")

test_that("univariate fitting is working correctly",{
  pro = 1; mu=1; Sigma=3; delta=2; nu=10;
  para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
  set.seed(0)
  y <- rmixst(100,para)
  out <- fitmixst(y,g=1,para=para)
  
  expect_that(length(y), equals(100))
  
})

test_that("multivariate fitting is working correctly",{
  pro = 1; mu=list(c(0,0)); Sigma=list(diag(2)); 
  delta=list(c(0,0)); nu=c(10);
  para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
  set.seed(0)
  y <- rmixst(100,para)
  out <- fitmixst(y,g=1,method="self",para=para)
  
  expect_that(nrow(y), equals(100))
  
})