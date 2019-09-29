
// includes from the plugin
#include <RcppArmadillo.h>
#include "dmmixst.h"


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
  SEXP dmmixstC( SEXP y, SEXP pro, SEXP mu, SEXP Sigma, SEXP delta, SEXP nu, SEXP log) ;
}

// definition

RcppExport SEXP dmmixstC( SEXP y, SEXP pro, SEXP mu, SEXP Sigma, SEXP delta, SEXP nu, SEXP log ){
  BEGIN_RCPP
  
  arma::mat yy = Rcpp::as<arma::mat> (y);
  arma::vec pro_neu = Rcpp::as<arma::vec> (pro);
  Rcpp::List mu_neu (mu);
  Rcpp::List Sigma_neu (Sigma);
  Rcpp::List delta_neu (delta);
  arma::vec nu_neu = Rcpp::as<arma::vec> (nu);
  
  int n = yy.n_rows;
  int g = pro_neu.size();
  
  arma::vec yhat = arma::zeros<arma::vec>(n);
  
  for(int j = 0; j < n; j++){
    
    
    for(int i = 0; i < g; i++){
      yhat(j) += dmixstc(yy.row(j).t(), pro_neu(i), Rcpp::as<arma::vec>(mu_neu(i)), Rcpp::as<arma::mat>(Sigma_neu(i)), 
      arma::diagmat(Rcpp::as<arma::vec>(delta_neu(i))), nu_neu(i));
    }
    
  }
  
  
  if(as<bool>(log)){
    yhat = arma::log(yhat);
  }
  
  
  return wrap(yhat);
  
  END_RCPP
}


