#ifdef INSIDE
#include <Rcpp.h>
#include <RcppGSL.h>
#include <RInside.h>                    // for the embedded R via RInside
//#include "rcpp_hello_world.h"
#include "mixstEst.h"

using namespace Rcpp;
using namespace std;
int main(int argc, char *argv[]) {
    RInside R(argc, argv);              // create an embedded R instance
    
    SEXP y1;
    SEXP g;
    SEXP itermax;
    SEXP error;
    SEXP  pro;
    SEXP mu;
    SEXP Sigma;
    SEXP delta;
    SEXP nu;
    SEXP verbose;
    
    y1 = wrap(Rcpp::rnorm(100));
    g=wrap(1);
    itermax=wrap(100);
    error=wrap(1e-6);
    
    
    pro = wrap(0.5);
    mu = wrap(0);
    Sigma = wrap(1);
    delta = wrap(0);
    nu = wrap(5);
    verbose = wrap(0);
    
    
    
    SEXP out = fitmixstC(y1, g, itermax, error, pro, mu, Sigma, delta, nu, verbose);
    
    Rcpp::List asd(out);
    //asd = wrap(out);
    
    //cout << asd(0);
    
}
// create an embedded R instance
//SEXP s = rcpp_hello_world();
/* Language call("print",s);
 call.eval();
 
 gsl_min(5,5);
 
 }*/
#endif
