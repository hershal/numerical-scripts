/*******************************************************
Function to implement the steepest descent method to 
minimize a function g(x).

Inputs:  
  geval    Function to evaluate g(x), scalar
  dgeval   Function to evaluate dg/dx(x), n vector
  xk       Initial guess for x, n vector  
  maxIter  Maximum number of descent iterations
  tol      Tolerance parameter for ||dg/dx||

  a0       Initial guess for alpha (step size)
  maxSrch  Maximum number of alpha search steps


Outputs: 
  xk       Approx local minimum point of g(x)
  k        Number of descent steps taken


Note 1: This function is incomplete; you'll need to code 
the quadratic interpolation step as indicated below.

Note 2: The functions geval and dgeval are assumed 
to be defined externally (e.g. by calling program).
*******************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std ;

/*** Declare user-defined functions to be used ***/
void geval(vector&, double&) ;  //defined in program8.cpp
void dgeval(vector&, vector&) ; //defined in program8.cpp


int descent(vector& xk, int n, int maxIter, 
                      double tol, int maxSrch, double a0){
  int k=0, kSrch=0 ;  
  double gk, dgkNorm ;
  double a1, a2, a3, y1, y2, y3, ah, gh, ahh, ghh ;
  vector dgk(n), xh(n), xhh(n) ;

  geval(xk,gk) ; //evaluate g(xk)
  dgeval(xk,dgk) ;  //evaluate dg(xk)
  dgkNorm = vecL2Norm(dgk) ; // z_0
  
  while(dgkNorm>=tol && k<maxIter) {
    ah = a0 ;

    // z=1/z_0 \nabla g(x^{(0)})
    // x=x-alpha(minidx).*z;
    xh = xk - scaleVec(ah/dgkNorm,dgk) ;

    // g(2) = gg(x)
    geval(xh,gh) ;

    kSrch = 0 ;
    while( gh>=gk && kSrch<maxSrch ){ 
      ah = ah/2 ;
      xh = xk - scaleVec(ah/dgkNorm,dgk) ;
      geval(xh,gh) ;
      kSrch++ ;
    }
    if(kSrch>=maxSrch){
      cout << "Descent: max iter exceeded in alpha search" << endl;
    }

    a1 = 0 ; y1 = 0 ;
    a2 = 0 ; y2 = 0 ;
    a3 = 0 ; y3 = 0 ;
    ahh = ah ; ghh = gh ; xhh = xh ;

    // ahh= alpha(2) = 0.5
    // ah = alpha(3) = 1
    ahh=0.5;
    ah=1.0;

    // gk = g(1)=gg(x)
    // ghh= g(2)
    // gh = g(3)=gg(x-alpha(3).*z)

    // h1 = y1
    // h2 = y2
    // h3 = y3
    y1=(ghh-gk)/(ah);
    y2=(gh-ghh)/(ahh-ah);
    y3=(y2-y1)/(ahh);

    ahh=0.5*(ah-y1/y3);
    xhh=xk-scaleVec(ahh/dgkNorm,dgk);
    geval(xhh,ghh);
    xhh=xh-scaleVec(ah/dgkNorm, dgk);

    // xhh=
    // ghh=
    /*** compute quadratic interpolation quantities
         (using notation from class notes) a1,y1, 
         a2,y2, a3,y3 and ahh,xhh,ghh here ***/

    // h(alpha_3) < h(alpha_crit)
    if(ahh>=0 && ahh<=ah && ghh<gh){
      xk = xhh ;
      gk = ghh ;
    }
    else {
      xk = xh ;
      gk = gh ;
    }

    dgeval(xk,dgk) ;
    dgkNorm = vecL2Norm(dgk) ; 
    k++ ;
    cout << "Descent: |Gradg_k|_2 = " << dgkNorm << endl ;
  }


  if(dgkNorm < tol) {
    cout << "Descent: solution converged" << endl ;
  } 
  else {
    cout << "Descent: max iterations exceeded" << endl;
  }

  return k;
}

