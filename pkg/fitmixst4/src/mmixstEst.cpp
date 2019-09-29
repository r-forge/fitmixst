#include <RcppArmadillo.h>
#include <RcppGSL.h>
#ifndef __APPLE__
#include <omp.h>
#endif
using namespace Rcpp;

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>

#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>

//#define NDEBUG 1
//#define ARMA_NO_DEBUG

extern"C" {
  
  
  extern void mvtdst_(int* n,
  int* nu,
  double* lower,
  double* upper,
  int* infin,
  double* correl,
  double* delta,
  int* maxpts,
  double* abseps,
  double* releps,
  double* error,
  double* value,
  int* inform);
  
}
double pmvnorm(int* n,
int* nu,
double* lower,
double* upper,
int* infin,
double* correl,
double* delta, // non-central parameter
int* maxpts,    // param
double* abseps, // param
double* releps, // param
double* error,  // estimated abs. error. with 99% confidence interval
double* value,     // results store here.
int* inform)    // inform message goes here
{
  // #ifndef __APPLE__
  // #pragma omp critical
  // #endif
  {
    mvtdst_ (n, nu,
    lower, upper, infin, correl, delta,
    maxpts, abseps, releps,
    error, value, inform);
  }
  #ifndef NDEBUG
  printf ("error = %g, value = %g, inform = %d\n", *error, *value, *inform);
  #endif
  switch (*inform) {
    case 0:
    return *value;
    case 1:
    case 2:
    case 3:
    return -1.0;
  };
  
  return *value;
};

/**
 * * * * * * * * @return CDF of multivariate normal P ( X < bound ) where X ~ MVN(0, correlationMatrix)
 * * * * * * * */
 double pmvnorm_P(int n,
 double* bound,
 //double* delta,
 double* correlationMatrix, // (2,1), (3,1), (3,2) .....
 int nu_,
 double* error)
 {
   //int nu_ = 0;
   int maxpts_ = 25000;
   double abseps_ = 0.001;
   double releps_ = 0;
   
   
   double* lower = new double[n];
   int* infin = new int[n];
   double* delta = new double[n];
   
   int i = 0;
   for (i = 0; i < n; ++i) {
     infin[i] = 0; // (-inf, bound]
     lower[i] = 0.0;
     delta[i] = 0.0;
   }
   
   // return values
   double value_ = 0.0;
   int inform_ = 0.0;
   
   double ret = pmvnorm(&n, &nu_, lower, bound, infin, correlationMatrix, delta, &maxpts_, &abseps_, &releps_, error, &value_, &inform_);
   delete[] (lower);
   delete[] (infin);
   delete[] (delta);
   
   return ret;
 }
 
 //#include "/Users/jh/Documents/work/Fitting Skew-t/code/libMvtnorm/mvtnorm.h"
 
 // define interrupt function
 static void chkIntFn1(void *dummy) {
   R_CheckUserInterrupt();
 }
 
 // this will call the above in a top-level context so it wont longjmp-out of your context
 bool checkInterrupt1() {
   return (R_ToplevelExec(chkIntFn1, NULL) == FALSE);
 }
 
 
 NumericVector dmvtc(NumericVector x, NumericVector mu, NumericMatrix sigma, NumericVector df){
   Language call("dmt", x, Named("mean", mu ), Named("S", sigma), Named("df", df), Named("log", 0));
   return call.eval();
 }
 
 //multivariate t distribution
 double dmtc(arma::vec x, arma::vec mu, arma::mat sigma, double df){
   int p = x.size();
   if( sigma.has_nan() || !sigma.is_sympd() ) { throw( std::runtime_error("C++ error: EM algorithm not converging! Error in multivariate t distribution. Sigma not invertible.") );}
   arma::mat SInv = inv_sympd(sigma);
   arma::vec xmu = x-mu;
   double Q = arma::as_scalar(((SInv * xmu) % xmu).t() * arma::ones<arma::vec>(p));
   
   double logDet;
   double b;
   
   arma::log_det(logDet,b, sigma);
   double logPDF = Rf_lgammafn((df + p) / 2.0) - 0.5 * (p* log(M_PI  * df) + logDet) - Rf_lgammafn(df / 2.0) - 0.5 * (df + p) * log(1.0 + Q / df);
   
   //Rprintf("%f",exp(logPDF));
   
   return(exp(logPDF));
 }
 
 
 
 NumericVector pmvtc(arma::vec x, arma::vec mu, arma::mat sigma, int df){
   Language call("pmt", Named("x",wrap(x)), Named("mean", wrap(mu) ), Named("S", wrap(sigma)), Named("df", wrap(df)));
   return call.eval();
 }
 
 //NumericVector pmvtc(NumericVector x, NumericVector mu, NumericMatrix sigma, NumericVector df){
 //Language call("pmt", Named("x",x), Named("mean", mu ), Named("S", sigma), Named("df", df));
 //return call.eval();
 //}
 

#include <mvtnormAPI.h>
 
 arma::vec triangl(const arma::mat& X){
    int n = X.n_cols;
    arma::vec res(n * (n-1) / 2);
    for (int i = 0; i < n; ++i) {
       for (int j = 0; j < i; ++j) {
          res(j + i * (i-1) / 2) = X(i, j);
       }
    }
    return res;
 }
 
 double pmvnorm_cpp(arma::vec& bound,
                    arma::vec& lowertrivec,
                    int nu,
                    double abseps = 1e-3){
    
    int n = bound.n_elem;
    //int nu = 0;
    int maxpts = 25000;     // default in mvtnorm: 25000
    double releps = 0;      // default in mvtnorm: 0
    int rnd = 1;            // Get/PutRNGstate
    
    double* bound_ = bound.memptr();
    double* correlationMatrix = lowertrivec.memptr();
    double* lower = new double[n];
    int* infin = new int[n];
    double* delta = new double[n];
    
    for (int i = 0; i < n; ++i) {
       infin[i] = 0; // (-inf, bound]
       lower[i] = 0.0;
       delta[i] = 0.0;
    }
    
    // return values
    double error;
    double value;
    int inform;
    
    mvtnorm_C_mvtdst(&n, &nu, lower, bound_,
                     infin, correlationMatrix, delta,
                     &maxpts, &abseps, &releps,
                     &error, &value, &inform, &rnd);
                     delete[] (lower);
                     delete[] (infin);
                     delete[] (delta);
                     
                     return value;
 }
 
 double pmtc(arma::vec x, arma::vec mu, arma::mat sigma, int df){
    arma::vec d = arma::sqrt(1 / arma::diagvec(sigma));
    arma::mat D = arma::diagmat(d);
    arma::mat corr = D * sigma * D;
    arma::vec triang = triangl(corr);
    arma::vec bound = (x - mu) % d;
    int nu = df;
    
    double res;
    
#ifndef __APPLE__
#pragma omp critical
#endif
    res = pmvnorm_cpp(bound, triang, nu);
    return res;
 }
 
 //multivariate t distribution
 // double pmtc(arma::vec x, arma::vec mu, arma::mat sigma, int df){
 //   
 //   int p = mu.n_elem;
 //   arma::vec d = arma::sqrt(1 / arma::diagvec(sigma));
 //   arma::mat D = arma::diagmat(d);
 //   arma::mat corr = D * sigma * D;
 //   
 //   int n = (p-1) * p / 2;
 //   
 //   double upper [p];
 //   double corrIn [n];
 //   
 //   double error;
 //   
 //   int k = 0;
 //   for(int i = 0; i < p; i++){
 //     upper[i] = (x(i) - mu(i)) * d(i);
 //     for(int j = 0; j < p; j++){
 //       if(j > i){
 //         corrIn[k] = corr(j,i);
 //         k++;
 //       }
 //     }
 //   }
 //   return pmvnorm_P(p, upper, corrIn, df, &error);
 //   
 // }
 
 //definition of dmixstc function (density function)
 double dmixstc(arma::colvec xx, double pro2, arma::colvec mu_neu, arma::mat Sigma_neu, arma::mat delta_neu, double nu_neu){
   
   nu_neu = round(nu_neu);
   
   double yhat = 0.0;
   int p = mu_neu.size();
   
   arma::mat lambda = Sigma_neu + arma::pow(delta_neu, 2);
   if( lambda.has_nan() || !lambda.is_sympd() ) { throw( std::runtime_error("C++ error: EM algorithm not converging! Error in multivariate t distribution. Lambda not invertible.") );}
   arma::mat lambdaInv = arma::inv(lambda);
   arma::colvec correctedMu = xx - mu_neu;
   
   
   double density = dmtc(xx, mu_neu, lambda, nu_neu);
   //NumericVector density = dmvtc(wrap(arma::trans(xx)), wrap(arma::trans(mu_neu)), wrap(lambda), wrap(nu_neu));
   
   double term1 = std::sqrt((nu_neu + p) / (nu_neu + arma::as_scalar(arma::trans(correctedMu) * lambdaInv * correctedMu)));
   
   arma::mat sigNew = arma::eye(p,p) - arma::trans(delta_neu) * lambdaInv * delta_neu;
   arma::colvec upper = delta_neu * lambdaInv * correctedMu * term1;
   
   //NumericVector upper1 = wrap(arma::conv_to< std::vector<double> >::from(upper));
   
   //NumericVector distribution = pmvtc(upper1, wrap(arma::conv_to< std::vector<double> >::from(arma::zeros(p,1))), wrap(sigNew), wrap(nu_neu + p));
   double distribution = pmtc(upper, arma::zeros(p,1), sigNew, nu_neu + p);
   
   //yhat = pro2 * pow(2,p) * Rcpp::as<double>(density) * Rcpp::as<double>(distribution);
   yhat = pro2 * pow(2,p) * density * distribution;
   
   return yhat;
 }
 
 //first and secon moment for truncated t distribution
 Rcpp::List trunctm1(arma::colvec a, arma::colvec mu, arma::mat sigma, double nu){
   Language call("truncatedt", Named("a", a) , Named("mu", mu ), Named("sigma", sigma), Named("nu", nu));
   return call.eval();
 }
 
 int truncatedt(arma::vec* m2, arma::mat* m3, arma::vec a, arma::vec mu, arma::mat sigma, int nuInt){

   int p = a.size();

   double nu = (double) nuInt;

   *m2 = arma::zeros<arma::vec>(p);
   *m3 = arma::zeros<arma::mat>(p,p);

   arma::vec  m1(p);
   arma::mat  h(p,p);

   arma::vec  tmp2 = arma::zeros<arma::vec>(p);

   arma::vec astar(p - 1);
   arma::mat sigmastar(p - 1, p - 1);
   double nustar;

   arma::vec astarstar = arma::zeros<arma::vec>(p-2);
   arma::mat sigmastarstar = arma::zeros<arma::mat>(p-2, p-2);

   double tmp = 1;

   for(int i = 0; i <p; i++){
     arma::vec atmp = a;
     arma::vec mutmp = mu;
     arma::vec sigmatmp = sigma.col(i);
     arma::mat sigmatmp1 = sigma;


     atmp.shed_row(i);
     mutmp.shed_row(i);
     sigmatmp.shed_row(i);

     sigmatmp1.shed_row(i);
     sigmatmp1.shed_col(i);


     astar = (atmp - mutmp) - (a(i) - mu(i)) / sigma(i,i) * sigmatmp;
     sigmastar = (nu + 1 / sigma(i,i) * pow(a(i) - mu(i), 2)) / (nu - 1) * (sigmatmp1 - 1/sigma(i,i) * sigmatmp * sigmatmp.t());

     /*use lgamma instead of gamma for calculation to avoid numerical instability*/
     m1(i) = pow(2 * M_PI *sigma(i,i), -0.5) * pow((nu / (nu + 1 / sigma(i,i) * pow(a(i) - mu(i), 2))), (nu - 1) / 2.0) * exp(Rf_lgammafn((nu - 1) / 2.0) -
     Rf_lgammafn(nu / 2.0)) *sqrt(nu / 2.0) *
     pmtc(astar, arma::zeros<arma::vec>(p-1), sigmastar, nu - 1) /
     pmtc(a-mu, arma::zeros<arma::vec>(p), sigma, nu);

     for(int j=0; j < p; j++){
       if(i != j){
         arma::uvec indices;
         indices << i << j;

         arma::vec atmp1 = a;
         arma::vec mutmp1 = mu;
         arma::mat sigmatmp2 = sigma.cols(indices);
         arma::mat sigmatmp3 = sigma;

         //wenn element schon entfernt wurde ist die länge anders
         if(i < j){
           atmp1.shed_row(i);
           atmp1.shed_row(j-1);
           mutmp1.shed_row(i);
           mutmp1.shed_row(j-1);
           sigmatmp2.shed_row(i);
           sigmatmp2.shed_row(j-1);
           sigmatmp3.shed_row(i);
           sigmatmp3.shed_row(j-1);
           sigmatmp3.shed_col(i);
           sigmatmp3.shed_col(j-1);
         }

         if(i > j){
           atmp1.shed_row(i);
           atmp1.shed_row(j);
           mutmp1.shed_row(i);
           mutmp1.shed_row(j);
           sigmatmp2.shed_row(i);
           sigmatmp2.shed_row(j);
           sigmatmp3.shed_row(i);
           sigmatmp3.shed_row(j);
           sigmatmp3.shed_col(i);
           sigmatmp3.shed_col(j);
         }
         arma::mat sigmaindices = sigma(indices, indices);
         if(!sigmaindices.is_sympd() ) { throw( std::runtime_error("C++ error: EM algorithm not converging! Moments of truncated t cannot be calculated.") );}
         arma::mat sigmaInv = inv_sympd(sigmaindices);

         nustar = (double) nu + arma::as_scalar((a.elem(indices) - mu.elem(indices)).t() *  sigmaInv * (a.elem(indices) - mu.elem(indices)));


         if(p > 2){
           astarstar = atmp1 - mutmp1 - sigmatmp2 * sigmaInv * (a.elem(indices) - mu.elem(indices));
           sigmastarstar = nustar / (nu - 2) * (sigmatmp3 - sigmatmp2 * sigmaInv * sigmatmp2.t());
           tmp = pmtc(astarstar, arma::zeros<arma::vec>(p-2), sigmastarstar, nu-2);
         }

         h(i,j) = - 1 / (2 * M_PI * sqrt(sigma(i, i) * sigma(j, j) - pow(sigma(i, j), 2))) * (nu / (nu - 2)) * pow((nu / nustar),(nu / 2 -1)) * tmp /
         pmtc(a - mu, arma::zeros<arma::vec>(p), sigma, nu);

         tmp2(i) += sigma(i,j) * h(i,j);


       }
     }

   }

   for(int i = 0; i < p; i++){
     h(i,i) = 1 / sigma(i,i) * (a(i) - mu(i)) * m1(i) - 1 / sigma(i,i) * tmp2(i);
   }

   *m2 = mu - sigma * m1;
   *m3 = -mu * mu.t() + mu * (*m2).t() + (*m2) * mu.t() + nu / (nu - 2) *
   pmtc(a - mu, arma::zeros<arma::vec>(p), nu / (nu - 2) * sigma, nu - 2) /
   pmtc(a - mu, arma::zeros<arma::vec>(p), sigma, nu) *
   sigma - sigma * h * sigma;


   return 0;
 }
 
 
 
 
 

 arma::vec truncatedtm1(arma::vec a, arma::vec mu, arma::mat sigma, int nuInt){
    
   if(sigma.has_nan() || !sigma.is_sympd()) { throw( std::runtime_error("C++ error: truncatedm1, sigma.") );}

   int p = a.size();

   double nu = (double) nuInt;

   arma::vec m2(p);
   //*m3 = arma::zeros<arma::mat>(p,p);

   arma::vec  m1(p);
   arma::mat  h(p,p);

   arma::vec  tmp2 = arma::zeros<arma::vec>(p);

   arma::vec astar(p - 1);
   arma::mat sigmastar(p - 1, p - 1);
   double nustar;

   arma::vec astarstar = arma::zeros<arma::vec>(p-2);
   arma::mat sigmastarstar = arma::zeros<arma::mat>(p-2, p-2);

   double tmp = 1;

   for(int i = 0; i <p; i++){
     arma::vec atmp = a;
     arma::vec mutmp = mu;
     arma::vec sigmatmp = sigma.col(i);
     arma::mat sigmatmp1 = sigma;


     atmp.shed_row(i);
     mutmp.shed_row(i);
     sigmatmp.shed_row(i);

     sigmatmp1.shed_row(i);
     sigmatmp1.shed_col(i);


     astar = (atmp - mutmp) - (a(i) - mu(i)) / sigma(i,i) * sigmatmp;
     sigmastar = (nu + 1 / sigma(i,i) * pow(a(i) - mu(i), 2)) / (nu - 1) * (sigmatmp1 - 1/sigma(i,i) * sigmatmp * sigmatmp.t());
     
     if(astar.has_nan()) { throw( std::runtime_error("C++ error: truncatedm1, astar.") );}
     if(sigmastar.has_nan() || !sigmastar.is_sympd()) { throw( std::runtime_error("C++ error: truncatedm1, sigmastar.") );}

     m1(i) = pow(2 * M_PI *sigma(i,i), -0.5) * pow((nu / (nu + 1 / sigma(i,i) * pow(a(i) - mu(i), 2))), (nu - 1) / 2.0) * exp(Rf_lgammafn((nu - 1) / 2.0) -
     Rf_lgammafn(nu / 2.0)) *sqrt(nu / 2.0) *
     pmtc(astar, arma::zeros<arma::vec>(p-1), sigmastar, nu - 1) /
     pmtc(a-mu, arma::zeros<arma::vec>(p), sigma, nu);

   }
    
   if(mu.has_nan()) { throw( std::runtime_error("C++ error: truncatedm1, mu.") );}
   if(sigma.has_nan()) { throw( std::runtime_error("C++ error: truncatedm1, sigma after.") );}
   if(m1.has_nan()) { throw( std::runtime_error("C++ error: truncatedm1, m1.") );}

   return mu - sigma * m1;

 }


 arma::mat truncatedtm2(arma::vec a, arma::vec mu, arma::mat sigma, int nuInt, arma::vec m2){

   int p = a.size();

   double nu = (double) nuInt;

   //m2 = arma::zeros<arma::vec>(p);
   arma::mat m3(p,p);

   arma::vec  m1(p);
   arma::mat  h(p,p);

   arma::vec  tmp2 = arma::zeros<arma::vec>(p);

   arma::vec astar(p - 1);
   arma::mat sigmastar(p - 1, p - 1);
   double nustar;

   arma::vec astarstar = arma::zeros<arma::vec>(p-2);
   arma::mat sigmastarstar = arma::zeros<arma::mat>(p-2, p-2);
   double tmp = 1;

   for(int i = 0; i <p; i++){

     for(int j=0; j < p; j++){
       if(i != j){
         arma::uvec indices;
         indices << i << j;

         arma::vec atmp1 = a;
         arma::vec mutmp1 = mu;
         arma::mat sigmatmp2 = sigma.cols(indices);
         arma::mat sigmatmp3 = sigma;

         //wenn element schon entfernt wurde ist die länge anders
         if(i < j){
           atmp1.shed_row(i);
           atmp1.shed_row(j-1);
           mutmp1.shed_row(i);
           mutmp1.shed_row(j-1);
           sigmatmp2.shed_row(i);
           sigmatmp2.shed_row(j-1);
           sigmatmp3.shed_row(i);
           sigmatmp3.shed_row(j-1);
           sigmatmp3.shed_col(i);
           sigmatmp3.shed_col(j-1);
         }

         if(i > j){
           atmp1.shed_row(i);
           atmp1.shed_row(j);
           mutmp1.shed_row(i);
           mutmp1.shed_row(j);
           sigmatmp2.shed_row(i);
           sigmatmp2.shed_row(j);
           sigmatmp3.shed_row(i);
           sigmatmp3.shed_row(j);
           sigmatmp3.shed_col(i);
           sigmatmp3.shed_col(j);
         }

         arma::mat sigmaindices = sigma(indices, indices);
         if(!sigmaindices.is_sympd()) { throw( std::runtime_error("C++ error: EM algorithm not converging! Moments of truncated t cannot be calculated.") );}
         arma::mat sigmaInv = inv_sympd(sigmaindices);

         nustar = (double) nu + arma::as_scalar((a.elem(indices) - mu.elem(indices)).t() *  sigmaInv * (a.elem(indices) - mu.elem(indices)));


         if(p > 2){
           astarstar = atmp1 - mutmp1 - sigmatmp2 * sigmaInv * (a.elem(indices) - mu.elem(indices));
           sigmastarstar = nustar / (nu - 2) * (sigmatmp3 - sigmatmp2 * sigmaInv * sigmatmp2.t());
           tmp = pmtc(astarstar, arma::zeros<arma::vec>(p-2), sigmastarstar, nu-2);
         }

         h(i,j) = - 1 / (2 * M_PI * sqrt(sigma(i, i) * sigma(j, j) - pow(sigma(i, j), 2))) * (nu / (nu - 2)) * pow((nu / nustar),(nu / 2 -1)) * tmp /
         pmtc(a - mu, arma::zeros<arma::vec>(p), sigma, nu);

         tmp2(i) += sigma(i,j) * h(i,j);


       }
     }

   }

   for(int i = 0; i < p; i++){
     h(i,i) = 1 / sigma(i,i) * (a(i) - mu(i)) * m1(i) - 1 / sigma(i,i) * tmp2(i);
   }

   //m2 = mu - sigma * m1;
   m3 = -mu * mu.t() + mu * m2.t() + m2 * mu.t() + nu / (nu - 2) *
   pmtc(a - mu, arma::zeros<arma::vec>(p), nu / (nu - 2) * sigma, nu - 2) /
   pmtc(a - mu, arma::zeros<arma::vec>(p), sigma, nu) *
   sigma - sigma * h * sigma;


   return m3;
 }


 
 
 
 
 
 //struct for root search
 struct my_f_params {double k;};
 
 //function definition needed for gsl function
 double nuf1(double nu, void *p){
   struct my_f_params * params = (struct my_f_params *) p;
   double k =(params -> k);
   
   double y=log(nu/2)-Rf_digamma(nu/2)-k;
   return y;
 }
 
 
 
 //root function
 double roots1 (double k)
 {
   int status;
   int iter = 0, max_iter = 1000;
   const gsl_root_fsolver_type *T;
   gsl_root_fsolver *s;
   double r = 0;
   double x_lo = 0.1, x_hi = 10000.0;
   gsl_function F;
   struct my_f_params params = {k};
   
   F.function = &nuf1;
   F.params = &params;
   
   T = gsl_root_fsolver_brent;
   s = gsl_root_fsolver_alloc (T);
   gsl_root_fsolver_set (s, &F, x_lo, x_hi);
   
   do
   {
     iter++;
     status = gsl_root_fsolver_iterate (s);
     r = gsl_root_fsolver_root (s);
     x_lo = gsl_root_fsolver_x_lower (s);
     x_hi = gsl_root_fsolver_x_upper (s);
     status = gsl_root_test_interval (x_lo, x_hi, 0, 0.0001);
     
   }
   while (status == GSL_CONTINUE && iter < max_iter);
   
   gsl_root_fsolver_free (s);
   
   //if (status != GSL_SUCCESS)
   return r;
 }
 
 // declarations
 extern "C" {
   SEXP fitmixstMVC( SEXP y, SEXP g, SEXP itermax, SEXP error,
   SEXP pro, SEXP mu, SEXP Sigma, SEXP delta, SEXP nu, SEXP verbose) ;
 }
 
 // definition
 SEXP fitmixstMVC( SEXP y, SEXP g, SEXP itermax, SEXP error,
 SEXP pro, SEXP mu, SEXP Sigma, SEXP delta, SEXP nu, SEXP verbose)
 {
   BEGIN_RCPP
   
   NumericMatrix yy (y);
   int gg = as<int>((g));
   int iitermax=as<int>((itermax));
   double eerror=as<double>((error));
   NumericVector pro_neu (pro);
   Rcpp::List mu_neu (mu);
   Rcpp::List Sigma_neu (Sigma);
   Rcpp::List delta_neu (delta);
   NumericVector nu_neu (nu);
   int verb = as<int>((verbose));
   
   NumericVector pro2 = pro_neu;
   Rcpp::List mu2 = mu_neu;
   Rcpp::List delta2 = delta_neu;
   Rcpp::List Sigma2 = Sigma_neu;
   NumericVector nu2 = nu_neu;
   
   arma::mat yyArma = Rcpp::as<arma::mat>(yy);
 
   //logLik
   double lik, lik_neu;
   lik = 1;
   lik_neu = 1;
   
   int flag = 1;
   int flag1 = 1;
   int flag2 = 1;
   
   //Aitken acceleration for while condition
   arma::vec liktmp(3);
   double aitkenError = 1.0;
   double likInf, a;
   
   int p = yyArma.n_cols;
   int n = yyArma.n_rows;
   int iter = 0;
   
   
   arma::vec sumpro2(n);
   arma::mat pro4(n,gg);
   
   //add likConvergence Vector for plot
   std::vector <double> likConvergence(iitermax);
   
   //observed fisher
   arma::mat tmpA(p,p);
   arma::mat tmpB(p,p);
   
   int p1 = (p+1) * p /2;
   arma::mat s1pi(n,gg-1);
   arma::cube s1mu(gg,p,n);
   arma::cube s1sigma(gg, p1, n);
   arma::cube s1delta(gg,p,n);
   arma::mat s1nu(n,gg);
   int p2 = gg-1 + gg* (p + p1 + p + 1);
   arma::vec s1(p2);
   arma::mat empcov(p2, p2);
   
   //fabs(lik / lik_neu - 1)
   while(iter < iitermax &&  aitkenError > eerror)
   {
     //interrupt running while
     if (checkInterrupt1()) { break;};
     
     arma::vec pro5 = arma::zeros<arma::vec>(gg);

     sumpro2 = arma::zeros<arma::vec>(n);
     for(int h = 0; h < gg; h++){
       arma::colvec cpMu = Rcpp::as<arma::colvec>(mu2[h]);
       arma::mat cpDelta = arma::diagmat(Rcpp::as<arma::colvec>(delta2[h]));
       arma::mat cpSigma = Rcpp::as<arma::mat> (Sigma2[h]);
       double cpNu = nu2[h];
       double cpPro = pro2[h];
       
       if( cpSigma.has_nan() || !cpSigma.is_sympd()  ) { throw( std::runtime_error("C++ error: EM algorithm not converging! Error in multivariate t density.") );}
       
       #ifndef __APPLE__
       #pragma omp parallel for
       #endif
       for(int j = 0; j < n; j++){
         pro4(j,h) = dmixstc(yyArma.row(j).t(), cpPro, cpMu, cpSigma, cpDelta, cpNu); 
       }
       sumpro2 += pro4.col(h);
     }
     
     for(int h = 0; h < gg; h++){
       pro4.col(h) = pro4.col(h) / sumpro2;
       pro5(h) = sum(pro4.col(h));
     }
     
       for(int i = 0; i < gg; i++){
       
       arma::colvec muArma = Rcpp::as<arma::colvec>(mu2[i]);
       arma::mat SigmaArma = Rcpp::as<arma::mat>(Sigma2[i]);
       arma::colvec deltaArma = Rcpp::as<arma::colvec>(delta2[i]);
       arma::mat DeltaArma = arma::diagmat(deltaArma);
       
       arma::mat Omega, Lambda, OmegaInv;
       
       Omega = SigmaArma + arma::pow(DeltaArma, 2);
       if( Omega.has_nan() || !Omega.is_sympd()  ) { throw( std::runtime_error("C++ error: EM algorithm not converging! Omega not invertible.") );}
       OmegaInv = inv_sympd(Omega);
       Lambda = arma::eye(p,p) - DeltaArma.t() * OmegaInv * DeltaArma;
       
       NumericVector pro3(n);
       arma::mat q(n,p), y_star(n,p), term3(n, p);
       arma::colvec d(n), term1(n), term2(n), tmp(n);
       arma::cube term4(p,p,n);
       
       lik = lik_neu;
       
       lik_neu = 0;
       
       arma::vec t1vec(n), t2vec(n);
       //#pragma omp for ordered schedule(dynamic)

       #ifndef __APPLE__
       #pragma omp parallel for
       #endif
       for(int j = 0; j < n; j++){

         arma::colvec qRow = DeltaArma * OmegaInv * (yyArma.row(j).t() - muArma);
         
         q.row(j) = (arma::mat(qRow)).t();
         
         d[j] = arma::as_scalar((yyArma.row(j).t() - muArma).t() * OmegaInv * (yyArma.row(j).t() - muArma));
         
         y_star.row(j)  = q.row(j)  * sqrt((nu2[i] + p) / (nu2[i] + d[j] ));

         arma::mat nudLambda = ((nu2[i] + d[j]) / (nu2[i] + p + 2)) * Lambda;
         if( Lambda.has_nan() || !nudLambda.is_sympd() ) { throw( std::runtime_error("C++ error: EM algorithm not converging! Scale matrix for distribution function is not invertible.") );}

		   t1vec(j) = pmtc(q.row(j).t() * sqrt((nu2[i] + p + 2) / (nu2[i] + d[j])), arma::zeros(p,1), Lambda, round(nu2[i] + p + 2));
		   t2vec(j) = pmtc(y_star.row(j).t(), arma::zeros(p,1), Lambda, round(nu2[i] + p));
         
         term2[j] = (nu2[i] + p) / (nu2[i] + d[j]) * t1vec(j) / t2vec(j);
         
         term1[j] = term2[j] - log((nu2[i] + d[j]) / 2) -
         (nu2[i] + p) / (nu2[i] + d[j]) + Rf_digamma((nu2[i] + p) / 2);
         
         arma::vec S2;
         arma::mat S3;
         
         arma::mat checkMat = ((nu2[i] + d[j]) / (nu2[i] + p + 2)) * Lambda;

         if( checkMat.has_nan() ||  !checkMat.is_sympd() ) {
         Rprintf("Sigma: ");
         for(int nn = 0; nn < p; nn++){
            for(int kk = 0; kk < p; kk++){
               Rprintf("%f", Rcpp::as<arma::mat>(Sigma_neu[i])(nn,kk));
               if( (nn < p-1) || (kk < p-1) )
                  Rprintf(", ");
            }
         }
         }
         if(  checkMat.has_nan() || !checkMat.is_sympd()) { throw( std::runtime_error("C++ error: EM algorithm not converging! Scale matrix for moments is not invertible.") );}
         // int truncerr = truncatedt(&S2, &S3, arma::zeros(p,1), -arma::conv_to<arma::colvec>::from( q.row(j) ),
         // ((nu2[i] + d[j]) / (nu2[i] + p + 2)) * Lambda, round(nu2[i]) + p + 2);
         
         S2 = truncatedtm1(arma::zeros(p,1), -arma::conv_to<arma::colvec>::from( q.row(j) ),
                           checkMat, round(nu2[i]) + p + 2);

         S3 = truncatedtm2(arma::zeros(p,1), -arma::conv_to<arma::colvec>::from( q.row(j) ),
                           checkMat, round(nu2[i]) + p + 2, S2);

         if(S2.has_nan()) { throw( std::runtime_error("C++ error: EM algorithm not converging! Error in truncatedt, S2.") );}
         if(S3.has_nan()) { throw( std::runtime_error("C++ error: EM algorithm not converging! Error in truncatedt, S3.") );}
         

         term3.row(j) = - term2(j) * S2.t();
         term4.slice(j) = term2(j) * S3;
         
         // arma::mat check12 = S3;
         // Rprintf("\n");
         // Rprintf("slice:%i", j);
         
         #ifndef __APPLE__
         #pragma omp atomic
         #endif
         lik_neu += log(sumpro2[j]);

       }

       if(term3.has_nan()) { throw( std::runtime_error("C++ error: EM algorithm not converging! Error in truncatedt, term3.") );}
       if(term4.has_nan()) { throw( std::runtime_error("C++ error: EM algorithm not converging! Error in truncatedt, term4.") );}
       
       //aitken acceleration
       liktmp(iter % 3) = lik_neu;
       
       if((iter > 1) && (i == 0)){
         a = (liktmp(iter % 3) - liktmp((iter - 1) % 3)) / (liktmp((iter -1)  % 3) - liktmp( (iter -2) % 3));
         likInf = liktmp( (iter -1)  % 3) + 1.0 / (1.0 - a) * (liktmp(iter % 3) - liktmp( (iter - 1) % 3));
         aitkenError = fabs(likInf - liktmp(iter % 3));
       }
       
       for(int h = 0; h < gg; h++){
         pro_neu(h) = pro5(h)/ n;
       }
       
       arma::colvec mutmp = arma::zeros<arma::colvec>(p);
       arma::colvec mutmp1 = arma::zeros<arma::colvec>(p);
       
       double mutmp2 = 0.0;
       
       for(int k = 0; k < n; k++){
         mutmp += pro4(k,i)*term2[k] *yyArma.row(k).t();
         mutmp1 += pro4(k,i)*term3.row(k).t();
         mutmp2 += pro4(k,i) * term2[k];
       }
       
       //for schleife über anzahl der komponenten (i von oben)
       if((mutmp2 == 0) || mutmp2 != mutmp2) { throw( std::runtime_error("C++ error: EM algorithm not converging! Division by zero for recalculating mu.") );}
       muArma = (mutmp - DeltaArma * mutmp1) / mutmp2;
       mu_neu[i] = muArma;

       arma::mat deltmp1 = arma::zeros<arma::mat>(p,p);
       arma::mat deltmp2 = arma::zeros<arma::mat>(p,p);
       
       double sumterm1pro3 = 0;
       double sumterm1pro3nu = 0;
       
       for(int k = 0; k < n; k++){
         deltmp1 += pro4(k,i) * term3.row(k).t() * (yyArma.row(k)  - muArma.t());
         if( term4.slice(k).has_nan() ) { throw( std::runtime_error("C++ error: EM algorithm not converging! term4 slice." + k) );}
         deltmp2 += pro4(k,i) * term4.slice(k);
         if( deltmp2.has_nan() ) { throw( std::runtime_error("C++ error: EM algorithm not converging! Delta tmp has nan.") );}
         
         sumterm1pro3 += term1[k] * pro4(k,i);
         sumterm1pro3nu += pro4(k,i)* (term2[k] - term1[k] -1);
       }

       if( SigmaArma.has_nan() || !SigmaArma.is_sympd() ) { throw( std::runtime_error("C++ error: EM algorithm not converging! Matrix sigma not invertible.") );}
       arma::mat sigInv = inv_sympd(SigmaArma);
       
       arma::mat sigInvDeltmp2 = sigInv % deltmp2;
       if(sigInv.has_nan() ) { throw( std::runtime_error("C++ error: EM algorithm not converging 1! Error in calculation of lambda.") );}
       if(deltmp2.has_nan() ) { throw( std::runtime_error("C++ error: EM algorithm not converging 2! Error in calculation of lambda.") );}
       if(!sigInvDeltmp2.is_sympd() ) {
          sigInvDeltmp2 = symmatu(sigInvDeltmp2);
       }
       
       if(!sigInvDeltmp2.is_sympd() ) { 
          throw( std::runtime_error("C++ error: EM algorithm not converging 4! Error in calculation of lambda.") );}
       deltaArma = inv_sympd(sigInv % deltmp2) * (sigInv % deltmp1) * arma::ones<arma::colvec>(p);
       delta_neu[i] = deltaArma;
       DeltaArma = arma::diagmat(deltaArma); 

       arma::mat sigtmp = arma::zeros<arma::mat>(p,p);
       double sigtmp2 = 0.0;
       arma::rowvec cenmu(p);
       
       for(int k = 0; k < n; k++){
         cenmu = yyArma.row(k)  - muArma.t();
         sigtmp += pro4(k,i) * cenmu.t() * cenmu * term2[k];
         sigtmp2 += pro4(k,i);
       }
       
       sigtmp += DeltaArma * deltmp2 * DeltaArma - DeltaArma * deltmp1 - deltmp1.t() * DeltaArma;
       sigtmp = arma::symmatu(sigtmp);
       
       //check if division by zero
       double sigtmpCheck = arma::as_scalar(sigtmp2);
       if(sigtmpCheck == 0  || sigtmpCheck != sigtmpCheck) { throw( std::runtime_error("C++ error: EM algorithm not converging! Division by zero for recalculating sigma.") );}
       SigmaArma = sigtmp / sigtmp2;

       arma::uword minInd;
       //determinant < eps
       while(!SigmaArma.is_sympd(0.001)) { 
         arma::vec diagSig = SigmaArma.diag();
         diagSig.min(minInd);
         SigmaArma.row(minInd) = arma::zeros(p).t();
         SigmaArma.col(minInd) = arma::zeros(p);
         SigmaArma(minInd, minInd) = 0.01;     
         if(flag){
           Rf_warning("Collapsed cluster detected! Scale matrix was changed to be postive definite (determinant of scale matrix equals 0.000001 now).");
           flag = 0;
         }
       }
       Sigma_neu[i] = SigmaArma;

       double nuSolve;
       nuSolve = sumterm1pro3nu / sigtmp2;
       //nu_neu[i]  = (roots1(nuSolve));
       
       if(nuSolve > 0.001){
         if(nuSolve < 0.271)
         nu_neu[i] = roots1(nuSolve);
         else{
           nu_neu[i] = 4;
           if(flag1)
           Rf_warning("Degrees of freedom getting to small and are set to 4 in this iteration.");
           flag1 = 0;
         }
       }
       else{
         nu_neu[i] = 1000;
         if(flag2)
         Rf_warning("Degrees of freedom getting to large and are set to 1000 in this iteration.");
         flag2 = 0;
       }

       //if while stops in the next step, we calculate here the empcov matrix && fisher 
       if((iter + 1 >= iitermax || aitkenError <= eerror) & (iter > 2) ){
         
        
         // observed fisher
         if( SigmaArma.has_nan() || !SigmaArma.is_sympd() ) { throw( std::runtime_error("C++ error: EM algorithm not converging! New matrix sigma not invertible.") );}
         arma::mat SigmaInvTmp = inv_sympd(SigmaArma);
         
         
         for(int j = 0; j < n; j++)
         {
           cenmu = yyArma.row(j)  - muArma.t();
           
           if(i < (gg - 1)){
             s1pi(j,i) = pro4(j,i) / (pro5(i) / n) - pro4(j,gg-1) / (pro5(gg-1) / n);
           }
           s1mu.slice(j).row(i) = (SigmaInvTmp * (pro4(j,i) * (term2(j) * cenmu.t() - DeltaArma *
           term3.row(j).t()))).t();
           
           tmpA = SigmaInvTmp * (pro4(j,i) * term2(j) *cenmu.t() * cenmu - pro4(j,i) *cenmu.t() * 
           term3.row(j) * DeltaArma - pro4(j,i)* DeltaArma * term3.row(j).t() * cenmu + 
           pro4(j,i) *DeltaArma * term4.slice(j) * DeltaArma) * 
           SigmaInvTmp - SigmaInvTmp;
           
           tmpB = (0.5 * (2* tmpA - arma::diagmat(tmpA)));
           
           int m = 0;
           for(int k = 0; k < p; k++){
             for(int l = 0; l < p; l++){
               if(k >= l){
                 s1sigma.slice(j).row(i)(m) = tmpB(k,l);
                 m++;
               }
             }
           }
           
           s1delta.slice(j).row(i) = ((SigmaInvTmp % (pro4(j,i) * term3.row(j).t() * cenmu)) * arma::zeros<arma::mat>(p,1) - 
           (SigmaInvTmp % term4.slice(j)) * deltaArma).t();
           
           s1nu(j,i) = 0.5 * pro4(j,i) * (log(0.5 * nu_neu(i)) +1 - Rf_digamma(0.5 * nu_neu(i))) + 0.5 * pro4(j,i) * (term1(j) - term2(j));
           
         }
       }
       
     }

     iter++;
     
     likConvergence[iter - 1] = lik_neu;
     if(verb){
       Rprintf("Iteration: %i", iter);
       Rprintf("\nrelative error: %lf", aitkenError);
       Rprintf("\nLoglikelihood: %f", lik_neu);
       Rprintf("\n\n");
       
       Rprintf("\npro: ");
       for(int i = 0; i < gg; i++){
         Rprintf("%f", pro_neu[i]);
         if(i < gg-1)
         Rprintf(", ");
       }
       
       Rprintf("\n\n");
       for(int i = 0; i < gg; i++){
         Rprintf("mu%i: ", i);
         for(int j = 0; j < p; j++){
           Rprintf("%f", Rcpp::as<arma::colvec>(mu_neu[i])(j));
           if(j < p-1)
           Rprintf(", ");
         }
         Rprintf("\n");
       }
       
       Rprintf("\n");
       for(int i = 0; i < gg; i++){
         Rprintf("delta%i: ", i);
         for(int j = 0; j < p; j++){
           Rprintf("%f", Rcpp::as<arma::colvec>(delta_neu[i])(j));
           if(j < p-1)
           Rprintf(", ");
         }
         Rprintf("\n");
       }
       
       Rprintf("\n");
       for(int i = 0; i < gg; i++){
         Rprintf("Sigma%i: ", i);
         for(int j = 0; j < p; j++){
           for(int k = 0; k < p; k++){
             Rprintf("%f", Rcpp::as<arma::mat>(Sigma_neu[i])(j,k));
             if( (j < p-1) || (k < p-1) )
             Rprintf(", ");
           }
         }
         Rprintf("\n");
       }
       
       Rprintf("\nnu: ");
       for(int i = 0; i < gg; i++){
         Rprintf("%f", nu_neu[i]);
         if(i < gg-1)
         Rprintf(", ");
       }
       
       Rprintf("\n\n");
       
     }
     
   }
   
   empcov = arma::zeros<arma::mat>(p2,p2);
   arma::mat empcovRet = arma::zeros<arma::mat>(p2,p2);
   
   //calculate fisher
   int fisher = 1;
   
   if(fisher){
   
   for(int j = 0; j < n; j++){
#ifndef __APPLE__
#pragma omp parallel for
#endif
     for(int i = 0; i < gg; i++){
       if(i < gg - 1){
         s1(i) = s1pi(j, i);
       }
       for(int k = 0; k < p; k++){
         s1(gg-1 + i * p + k) = s1mu(i,k,j);
       }
       
       for(int k = 0; k < p1; k++){
         s1(gg-1 + p * gg + i * p1 + k) = s1sigma(i,k,j);
       }
       
       for(int k = 0; k < p; k++){
         s1(gg-1 + p * gg + p1 * gg + i * p + k) = s1delta(i,k,j);
       }
       
       s1(gg - 1 + 2* p * gg + p1 * gg + i) = s1nu(j,i);
       
     }
     
     empcov += s1 * s1.t();
     
     
   }
   
   }

   //remove empty elements
   likConvergence.erase (likConvergence.begin() + iter,likConvergence.end());
   if(iter == iitermax){
     Rf_warning("Convergence criterion is not met. Stopping because of itermax.");
   }
   
   if(empcov.is_sympd()){
     empcovRet = inv_sympd(empcov);
   }else{
     Rf_warning("Calculating the empricical fisher information failed. Matrix not invertible.");
   }

   return Rcpp::List::create(Rcpp::Named("pro") = pro_neu,
   Rcpp::Named("mu") = mu_neu,
   Rcpp::Named("Sigma") = Sigma_neu,
   Rcpp::Named("delta") = delta_neu,
   Rcpp::Named("nu") = Rcpp::as<Rcpp::List>(nu_neu),
   Rcpp::Named("logLik") = lik_neu,
   Rcpp::Named("iter") = iter,
   Rcpp::Named("empcov") = empcovRet,
   Rcpp::Named("posteriori") = pro4,
   Rcpp::Named("p") = p,
   Rcpp::Named("g") = gg,
   Rcpp::Named("likConvergence") = likConvergence);
   
   END_RCPP
 }
 
 
 
