
// includes from the plugin
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <iostream>
#include <fstream>
#ifdef __linux
#include <omp.h>
#endif
#include <unistd.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

// define interrupt function
static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
bool checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

//definition of dmixstc function (density function)
double dmixstc(double xx, double pro2,
double mu_neu, double Sigma_neu,
double  delta_neu, double nu_neu){
  
  double yhat = 0.0;
  int p = 1;
  
  double lambda = Sigma_neu + pow(delta_neu, 2);
  double correctedMu = xx - mu_neu;
  double density = Rf_dt(correctedMu / sqrt(lambda), nu_neu, 0);
  double term1 = sqrt((nu_neu + p) / (nu_neu + correctedMu / lambda * correctedMu));
  
  double distribution = Rf_pt(delta_neu / lambda * correctedMu * term1 / sqrt(1-delta_neu / lambda * delta_neu), nu_neu + p, 1, 0);
  
  yhat = pro2 * 2 * density / sqrt(Sigma_neu + delta_neu * delta_neu) * distribution;
  
  return yhat;
}


//first moment for truncated t distribution
double trunctm1(double mu, double sigma, double nu){
  double tmp;
  double m1;
  //if(nu < 50.0){
  tmp = exp(Rf_lgammafn((nu + 1) / 2) - Rf_lgammafn(nu / 2)) / sqrt(PI * nu);
  //}
  //else{
  //  tmp = 1 / sqrt(PI * nu);
  //}
  
  double distribution = Rf_pt((0 - mu) / sigma, nu, 0, 0);
  m1 = mu + (tmp * nu * sigma / (nu - 1)) / (distribution * pow((1 + pow(mu, 2) / (nu * pow(sigma, 2))), ((nu - 1) / 2)));
  return m1;
}

//second moment of truncated t distribution
double trunctm2(double mu, double sigma, double nu){
  double tmp, tmp2;
  //if(nu  < 50.0){
  tmp = exp(Rf_lgammafn((nu + 1) / 2) - Rf_lgammafn(nu / 2)) / sqrt(PI * nu);
  tmp2 = exp(Rf_lgammafn((nu - 1) / 2) - Rf_lgammafn((nu - 2) / 2)) / sqrt(PI * (nu - 2));
  //}
  //else{
  //  tmp = 1 / sqrt(PI * nu);
  //  tmp2 = 1 / sqrt(PI * (nu - 2));
  //}
  
  double m2 = -pow(sigma, 2) * nu - pow(mu, 2) + sigma * nu * tmp / tmp2 * sqrt(nu / (nu - 2)) * sigma *
  (Rf_pt((0 - mu) / (sqrt(nu / (nu - 2)) * sigma), nu - 2, 0, 0)) / (Rf_pt((0 - mu) / sigma, nu, 0, 0))
  + 2 * mu * trunctm1(mu, sigma, nu);
  return m2;
}

//struct for root search
struct my_f_params {double k;};

//function definition needed for gsl function
double nuf(double nu, void *p){
  struct my_f_params * params = (struct my_f_params *) p;
  double k =(params -> k);
  
  double y=log(nu/2)-Rf_digamma(nu/2)-k;
  return y;
}


//searching root for with right side k
double roots (double k)
{
  int status;
  int iter = 0, max_iter = 1000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  double x_lo = 0.1, x_hi = 10000.0;
  gsl_function F;
  struct my_f_params params = {k};
  
  F.function = &nuf;
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
  SEXP fitmixstC( SEXP y, SEXP g, SEXP itermax, SEXP error,
  SEXP pro, SEXP mu, SEXP Sigma, SEXP delta, SEXP nu, SEXP verbose, SEXP bbeta) ;
}

// definition
SEXP fitmixstC( SEXP y, SEXP g, SEXP itermax, SEXP error,
SEXP pro, SEXP mu, SEXP Sigma, SEXP delta, SEXP nu, SEXP verbose, SEXP beta)
{
  BEGIN_RCPP
  
  std::vector<double> yy = Rcpp::as<std::vector< double > >(y);
  int gg = as<int>((g));
  int iitermax=as<int>((itermax));
  double eerror=as<double>(clone(error));
  std::vector<double> pro_neu = Rcpp::as<std::vector< double > >(pro);
  std::vector<double> mu_neu = Rcpp::as<std::vector< double > >(mu);
  std::vector<double> Sigma_neu = Rcpp::as<std::vector< double > >(Sigma);
  std::vector<double> delta_neu = Rcpp::as<std::vector< double > >(delta);
  std::vector<double> nu_neu = Rcpp::as<std::vector< double > > (nu);
  int verb = as<int>(clone(verbose));
  double bbeta = Rcpp::as< double > (beta);  
  
  int p=1;
  int n=yy.size();
  double lik = 1.0;
  double lik_neu = 1.0;
  int flag = 1;
  int flag1 = 1;
  int flag2 = 1;
  
  int iter = 0;
  
  std::vector<double> mu1(gg), Sigma1(gg), delta1(gg), nu1(gg), pro1(gg),
  mu2(gg), Sigma2(gg), delta2(gg), nu2(gg), pro2(gg);
  
  
  //Aitken acceleration for while condition
  arma::vec liktmp(3);
  double aitkenError = 1.0;
  double likInf, a;
  
  //add likConvergence Vector for plot
  std::vector <double> likConvergence(iitermax);
  
  std::vector< std::vector<double> > test1(gg), test2(gg), test3(gg), test4(gg), pro3(gg);
  while( (iter < iitermax) && (aitkenError > eerror) || (iter < 2) )
  {
    //interrupt running while
    if (checkInterrupt()) { break;};
    
    pro2 = pro_neu;
    mu2=mu_neu;
    delta2=delta_neu;
    Sigma2=Sigma_neu;
    nu2=nu_neu;
    
    std::vector<double> k(gg), Omega(gg), Lambda(gg);
    
    std::vector< std::vector<double> > q(gg), d(gg),
    tmp(gg), S2(gg), S3(gg), y_star(gg), sumpro2(gg), sumpro2ann(gg);
    
    std::vector<double> muNom(gg), muDenom(gg), deltaNom(gg), deltaDenom(gg), SigmaNom(gg),
    SigmaDenom(gg), nuNom(gg), nuDenom(gg), sumpro3(gg), sumpro3test1(gg), SigmaNom1(gg),
    SigmaNom2(gg), SigmaNom3(gg), dummy1(n), dummy2(n), dummy3(n), dummy4(n),
    dummy5(n), dummy6(n), dummy7(n), dummy8(n), dummy9(n), dummy10(n), dummy11(n), dummy12(n), dummy13(n);
    
    for(int i = 0; i < gg; i++)
    {
      test1[i] = dummy1;
      test2[i] = dummy2;
      test3[i] = dummy3;
      test4[i] = dummy4;
      pro3[i] = dummy5;
      q[i] = dummy6;
      d[i] = dummy7;
      tmp[i] = dummy8;
      S2[i] = dummy9;
      S3[i] = dummy10;
      y_star[i] = dummy11;
      sumpro2[i] = dummy12;
      sumpro2ann[i] = dummy13;
      
      muNom[i] = 0;
      muDenom[i] = 0;
      sumpro3[i] = 0;
      sumpro3test1[i] = 0;
      deltaNom[i] = 0;
      deltaDenom[i] = 0;
      
      Omega[i] = Sigma2[i] + pow(delta2[i], 2);
      Lambda[i] = 1 - pow(delta2[i], 2) / Omega[i];
      
      lik = lik_neu;
      
      lik_neu = 0;
      
      //#ifdef _WIN64
      //#pragma omp parallel for
      //#elif _WIN32
      //#pragma omp parallel for
      #ifdef __linux
      #pragma omp parallel for
      #endif
      
      for(int j = 0; j < n; j++){
        
        sumpro2[i][j] =0;
        sumpro2ann[i][j] =0;
        for(int h = 0; h < gg; h++){
          sumpro2[i][j] += dmixstc(yy[j], pro2[h], mu2[h], Sigma2[h], delta2[h], nu2[h]);
          sumpro2ann[i][j] += pow(dmixstc(yy[j], pro2[h], mu2[h], Sigma2[h], delta2[h], nu2[h]), bbeta);
        }
        
        #ifdef __linux
        #pragma omp atomic
        #endif
        lik_neu += log(sumpro2[i][j]);
        
        #ifdef __linux
        #pragma omp critical
        #endif
        {
        pro3[i][j] = pow(dmixstc(yy[j], pro2[i], mu2[i], Sigma2[i], delta2[i], nu2[i]), bbeta) / sumpro2ann[i][j];
        }
        
        q[i][j] = delta2[i] * (yy[j] - mu2[i]) / Omega[i];
        
        d[i][j] = (yy[j] - mu2[i]) * (yy[j] - mu2[i]) / Omega[i];
        y_star[i][j]  = q[i][j]  * sqrt((nu2[i] + p) / (nu2[i] + d[i][j] ));
        
        double t1 = Rf_pt(q[i][j] * sqrt((nu2[i] + p + 2) / (nu2[i] + d[i][j])) / sqrt(Lambda[i]), nu2[i] + p + 2, 1, 0);
        double t2 = Rf_pt(y_star[i][j] / sqrt(Lambda[i]), nu2[i] + p, 1, 0);
        test2[i][j] = (nu2[i] + p) / (nu2[i] + d[i][j]) * t1 / t2;
        
        test1[i][j] = (test2[i][j] - log((nu2[i] + d[i][j]) / 2) -
        (nu2[i] + p) / (nu2[i] + d[i][j]) + Rf_digamma((nu2[i] + p) / 2));
        
        
        tmp[i][j] = (nu2[i] + p) / (nu2[i] + d[i][j]) * Rf_pt(q[i][j] /
        sqrt((nu2[i] + d[i][j]) / (nu2[i] + p + 2) * Lambda[i]),
        nu2[i] + p + 2, 1, 0) /
        Rf_pt(y_star[i][j] / sqrt(Lambda[i]), nu2[i] + p, 1, 0);
        
        //moments of truncated t distribution
        S2[i][j] = trunctm1(q[i][j], sqrt(((nu2[i] + d[i][j]) /
        (nu2[i] + p + 2)) * Lambda[i]) , nu2[i] + p + 2);
        
        S3[i][j] = trunctm2(q[i][j], sqrt(((nu2[i] + d[i][j]) /
        (nu2[i] + p + 2)) * Lambda[i]), nu2[i] + p + 2);
        
        
        test3[i][j] = tmp[i][j] * S2[i][j];
        test4[i][j] = tmp[i][j] * S3[i][j];
        
        
        //#ifdef _WIN64
        //#pragma omp critical
        //#elif _WIN32
        //#pragma omp critical
        #ifdef __linux
        #pragma omp critical
        #endif
        {
          muNom[i] += pro3[i][j] * (test2[i][j] * yy[j] - delta2[i] * test3[i][j]);
          muDenom[i] += pro3[i][j] * test2[i][j];
          
          
          deltaDenom[i] += pro3[i][j] * test4[i][j];
          
          sumpro3[i] += pro3[i][j];
          
          sumpro3test1[i] += pro3[i][j] * test1[i][j];
        }
        
      }
      
      //aitken acceleration
      liktmp(iter % 3) = lik_neu;
      
      if(iter > 1 && i == 0){
        a = (liktmp(iter % 3) - liktmp((iter - 1) % 3)) / (liktmp((iter -1)  % 3) - liktmp( (iter -2) % 3));
        likInf = liktmp( (iter -1)  % 3) + 1 / (1 - a) * (liktmp(iter % 3) - liktmp( (iter - 1) % 3));
        aitkenError = fabs(likInf - liktmp(iter % 3));
      }
      
      //updating parameter
      
      if(muDenom[i] == 0) { throw( std::runtime_error("C++ error: EM algorithm not converging! Division by zero for recalculating mu.") );}
      
      mu1[i] = muNom[i] / muDenom[i];
      
      
      deltaNom[i] =0.0;
      for(int l = 0; l < n; l++) {
        deltaNom[i] += pro3[i][l] * (yy[l] - mu1[i]) * test3[i][l];
      }
      
      double d1 = deltaNom[i];
      
      if(deltaNom[i] == 0 || d1 != d1 || deltaDenom[i] == 0) { throw( std::runtime_error("C++ error: EM algorithm not converging! Division by zero for recalculating sigma.") );}
      delta1[i] = deltaNom[i] / deltaDenom[i];
      
      
      for(int j = 0; j < n; j++){
        SigmaNom[i] += pro3[i][j] * (pow(delta1[i], 2) * test4[i][j] - 2 * (yy[j] - mu1[i]) *
        test3[i][j] * delta1[i] + pow((yy[j] - mu1[i]), 2) * test2[i][j]);
      }
      
      if(sumpro3[i] == 0) { throw( std::runtime_error("C++ error: EM algorithm not converging!") );}
      
      Sigma1[i] = SigmaNom[i] / sumpro3[i];
      
      if(Sigma1[i] < 0.000001 ) { 
        Sigma1[i] = 0.000001;     
        if(flag){
          Rf_warning("Collapsed cluster detected! Scale parameter was changed (scale parameter equals 0.000001 now).");
          flag = 0;
        }
      }
      
      
      k[i] = - log(sumpro3[i] / muDenom[i]) -
      sumpro3test1[i] / sumpro3[i];
      
      //Rprintf("\tk: %f", k[i]);
      
      if(k[i] > 0.001){
        if(k[i] < 0.271)
        nu1[i] = roots(k[i]);
        else{
          nu1[i] = 4;
          if(flag1)
          Rf_warning("Degrees of freedom getting to small and are set to 4 in this iteration.");
          flag1 = 0;
        }
      }
      else{
        nu1[i] = 1000;
        if(flag2)
        Rf_warning("Degrees of freedom getting to large and are set to 1000 in this iteration.");
        flag2 = 0;
      }
      
      pro1[i] = sumpro3[i] / n;
      
    }
    //Rprintf("\n");    
    
    
    
    mu_neu = mu1;
    delta_neu = delta1;
    Sigma_neu = Sigma1;
    nu_neu = nu1;
    pro_neu = pro1;
    
    
    //		std::vector<double> liktmp(n);
    //
    //		lik = lik_neu;
    //		lik_neu = 0.0;
    //
    //		for(int i = 0; i < n; i++){
    //			liktmp[i] = 0;
    //			for(int j = 0; j < gg; j++){
    //
    //				liktmp[i] +=  dmixstc(yy[i], pro2[j], mu2[j], Sigma2[j], delta2[j], nu2[j]);
    //
    //			}
    //
    //			lik_neu += log(liktmp[i]);
    //
    //		}
    
    iter++;
    
    likConvergence[iter - 1] = lik_neu;
    if(verb){
      Rprintf("Iteration: %i", iter);
      Rprintf("\nrelative error: %f", aitkenError);
      Rprintf("\nLoglikelihood: %f", lik_neu);
      
      Rprintf("\npro: ");
      for(int i = 0; i < gg; i++){
        Rprintf("%f \t", pro_neu[i]);
      }
      
      Rprintf("\nmu: ");
      for(int i = 0; i < gg; i++){
        Rprintf("%f \t", mu_neu[i]);
      }
      
      Rprintf("\ndelta: ");
      for(int i = 0; i < gg; i++){
        Rprintf("%f \t", delta_neu[i]);
      }
      
      Rprintf("\nSigma: ");
      for(int i = 0; i < gg; i++){
        Rprintf("%f \t", Sigma_neu[i]);
      }
      
      Rprintf("\nnu: ");
      for(int i = 0; i < gg; i++){
        Rprintf("%f \t", nu_neu[i]);
      }
      Rprintf("\n\n");
    }
    
  }
  
  int fisher =1;
  
  NumericMatrix ecov(5*gg-1,5*gg-1);
  
  if(fisher )    {
    //calculate observed fisher information
    NumericMatrix s1pi(n,gg-1), s1mu(n,gg), s1sigma(n,gg), s1delta(n,gg),
    s1nu(n,gg), empcov(5*gg-1,5*gg-1), emp(5*gg-1,5*gg-1);
    
    NumericVector s1(gg-1+gg*4);
    
    //#ifdef _WIN64
    //#pragma omp parallel for
    //#elif _WIN32
    //#pragma omp parallel for
    #ifdef __linux
    #pragma omp parallel for
    #endif
    for (int j = 0; j < n; j++){
      for(int i = 0; i < gg; i++){
        if(i < (gg - 1)){
          s1pi(j,i) = pro3[i][j] / pro_neu[i] - pro3[gg-1][j] / pro_neu[gg-1];
        }
        s1mu(j,i) = pro3[i][j] * (test2[i][j] * (yy[j] - mu_neu[i]) - delta_neu[i] *
        test3[i][j]) /pow(Sigma_neu[i], 2);
        
        s1sigma(j,i) = pro3[i][j]* (pow((yy[j] - mu_neu[i]), 2) * test2[i][j]-
        2.0 * (yy[j] - mu_neu[i]) * test3[i][j] * delta_neu[i] +
        pow(delta_neu[i], 2) * test4[i][j]) / pow(Sigma_neu[i], 3) - pro3[i][j] /
        Sigma_neu[i];
        
        s1delta(j,i) = pro3[i][j] * ((yy[j] - mu_neu[i]) * test3[i][j] -
        delta_neu[i] * test4[i][j]) / pow(Sigma_neu[i],2);
        
        s1nu(j,i) = 0.5 * pro3[i][j] *(log(0.5 * nu_neu[i]) -
        Rf_digamma(0.5 * nu_neu[i]) + 1 + test1[i][j] - test2[i][j]);
        
      }
      
    };
    
    
    
    for (int i = 0; i < n; i++ ){
      for(int l = 0; l < gg; l++){
        if(l < gg - 1){
          s1(l) = s1pi(i, l);
        }
        s1(gg-1+l) = s1mu(i, l);
        s1(2*gg-1+l) = s1sigma(i, l);
        s1(3*gg-1+l) = s1delta(i, l);
        s1(4*gg-1+l) = s1nu(i, l);
        
      }
      
      for (int j = 0; j < 5 * gg - 1; j++){
        for (int k = 0; k < 5 * gg - 1; k++){
          emp(j,k) =  s1(j) * s1(k);
        }
        
      }
      empcov += emp;
    }
    
    
    //invert matrix for standard errors
    arma::mat ecov1 = Rcpp::as<arma::mat>(empcov);
    arma::mat ecov2;
    if(det(ecov1) > 0){
      ecov2 = ecov1.i();
      ecov = wrap(ecov2);
    }else{
      Rf_warning("Calculating the empricical fisher information failed. Matrix not invertible.");
    }
    
  }
  
  
  //remove empty elements
  likConvergence.erase (likConvergence.begin() + iter,likConvergence.end());
  
  if(iter == iitermax){
    Rf_warning("Warning: Convergence criterion is not met. Stopping because of itermax.");
  }
  
  
  return Rcpp::List::create(Rcpp::Named("pro") = pro_neu,
  Rcpp::Named("mu") = mu_neu,
  Rcpp::Named("Sigma") = Sigma_neu,
  Rcpp::Named("delta") = delta_neu,
  Rcpp::Named("nu") = nu_neu,
  Rcpp::Named("logLik") = lik_neu,
  Rcpp::Named("iter") = iter,
  Rcpp::Named("empcov") = ecov,
  Rcpp::Named("posteriori") = pro3,
  Rcpp::Named("p") = 1,
  Rcpp::Named("g") = gg,
  Rcpp::Named("likConvergence") = likConvergence);
  
  END_RCPP
}



