// includes from the plugin

#include <RcppArmadillo.h>
#include <RcppGSL.h>

#include <gsl/gsl_deriv.h>



#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes
//
//definition of dmixstc function (density function)
NumericVector dmixstc1(NumericVector xx, NumericVector pro2, NumericVector mu_neu, NumericVector
                       Sigma_neu, NumericVector delta_neu, NumericVector nu_neu, bool log){
    
    
    int l=xx.size();
    NumericVector yhat (l);
    
    int m=pro2.size();
    int p=1;
    
    
    for(int j=0;j<m;j++){
        yhat+=pro2(j)*2*dt((xx-mu_neu(j))/sqrt(Sigma_neu(j)+delta_neu(j)*delta_neu(j)),
                           nu_neu(j))/sqrt(Sigma_neu(j)+delta_neu(j)*delta_neu(j))*pt(delta_neu(j)*
                                                                                      1/(Sigma_neu(j)+delta_neu(j)*delta_neu(j))* (xx-mu_neu(j)) *
                                                                                      sqrt((nu_neu(j)+p)/(nu_neu(j)+(xx-mu_neu(j))*1/(Sigma_neu(j)+delta_neu(j)*
                                                                                                                                      delta_neu(j))*(xx-mu_neu(j))))/sqrt(1-delta_neu(j)*
                                                                                                                                                                          1/(Sigma_neu(j)+delta_neu(j)*delta_neu(j))*delta_neu(j)),nu_neu(j)+p);
    };
    
    if(log){
        yhat = Rcpp::log(yhat);
    }
    
    return yhat;
    
}

//struct for derivative parameters
struct myparam {double y; NumericVector all; int gg; int index; bool log;};

//my deriv function which is needed for gsl
double myderiv (double x, void *param) {
    struct myparam * params = (struct myparam *) param;
    double y = (params -> y);
    NumericVector all =(params -> all);
    double gg = (params -> gg);
    double index = (params -> index);
    double sumpro = 0.0;
    bool log = (params -> log);
    NumericVector pro(gg), mu(gg), Sigma(gg), delta(gg), nu(gg);
    
    all(index+1) = x;
    
    for(int i=0; i<gg; i++)
    {
        pro(i) = all(i);
        mu(i) = all(gg+i);
        Sigma(i) = all(2*gg+i);
        delta(i) = all(3*gg+i);
        nu(i) = all(4*gg+i);
        sumpro +=pro(i);
    }
    
    pro(0) = 1 + pro(0) - sumpro;
    
    double value=as<double>  (dmixstc1(wrap(y),pro,mu,Sigma,delta,nu,log));
    return value;
}


//take derivative of dmixst at parameter all
double deriv (double y, NumericVector all, int gg, int index, bool log)
{
    gsl_function F;
    double result, abserr;
    
    struct myparam derivparams = {y, all, gg, index, log};
    
    F.function = &myderiv;
    F.params = &derivparams;
    
    gsl_deriv_forward (&F, all(index+1), 1e-8, &result, &abserr);
    
    return result;
}




// declarations
extern "C" {
    SEXP stderrorC( SEXP xs, SEXP alls, SEXP empcov, SEXP ggs, SEXP logs) ;
}

// definition

RcppExport SEXP stderrorC( SEXP xs, SEXP alls, SEXP empcovs,SEXP ggs, SEXP logs){
    BEGIN_RCPP
    
    NumericVector x(xs);
    NumericVector all(alls);
    int gg =as<int>(ggs);
    bool log = as<bool>(logs);
    int l = x.size();
    NumericVector valuederiv(5*gg-1), sedist(l);
    NumericMatrix derivout(gg*5-1,l);
    NumericMatrix empcov(empcovs);
    
    //  for(int i=0; i<gg; i++){
    //  all(i) = pro_neu(i);
    //  all(gg + i) = mu_neu(i);
    //  all(2*gg + i) = Sigma_neu(i);
    //  all(3*gg + i) = delta_neu(i);
    //  all(4*gg + i) = nu_neu(i);
    //}
    
    arma::mat ecov1 = Rcpp::as<arma::mat>(empcov);
    
    for(int k=0; k<l; k++){
        for(int i=0; i<5*gg-1; i++){
            valuederiv(i) = deriv(x(k), all, gg, i, log);
            derivout(i,k) = valuederiv(i);
        }
        
        arma::colvec vderiv(valuederiv.begin(), valuederiv.size(),false);
        sedist(k) =  sqrt(as<double>( wrap(arma::trans(vderiv) * ecov1 * vderiv) ));
    }
    
    
    return sedist;
    
    END_RCPP
}
