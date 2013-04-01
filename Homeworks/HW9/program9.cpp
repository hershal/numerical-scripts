/*******************************************************
Program 9. Uses the Newton-RK4 shooting method to find 
an approximate solution of a two-point BVP of the form

              y'' = f(x,y,y'),  a<=x<=b  
              y(a) = alpha,  y(b) = beta

Inputs:  
  feval         Function to evaluate f(x,y,y')
  a,b           Interval params
  alpha,beta    Boundary value params
  t             Initial guess for y'(a)
  N             Number of grid points for RK4
  maxIter,tol   Control params for Newton

Outputs: 
  iter    Number of Newton steps taken
  t       Approx value of y'(a) 
  x       Grid point vector: x(j)=a+jh, j=0...N
  y       Approx soln vector: y(j)=soln at x(j), j=0...N

Note 1: For any given problem, the function feval must 
be changed.  This function computes f for any given x, y 
and y'.  (Below we use the symbol yp in place of y').

Note 2: For any given problem the parameters a, b, alpha,
beta, t, N, maxIter and tol must also be specified.

Note 3: To compile this program use the command (all on 
one line)

  c++ -o program9  matrix.cpp shootRK4.cpp program9.cpp

*******************************************************/
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

/*** Declare external function ***/
int shootRK4(int&, double&, double&, double&, double&, 
                double&, int&, double&, vector&, vector&) ;


/*** Define f(x,y,yp) function ***/
void feval(const double& x, const double& y, 
                              const double& yp, double& f){
  double M=0.2, L=1.3, R=-0.6 ; //define any constants
  double pi=4.0*atan(1.0) ; //the number pi

  f = M*pi + L*y*cos(x) + R*yp*yp ;
}


int main() {
  /*** Define problem parameters ***/
  double tol=1e-6 ;
  int N=20, maxIter=10, iter ;  
  double a=0, b=1, alpha=0, beta=2, t ; 
  vector x(N+1), y(N+1) ; 
  t=2.0 ; // initial guess of slope 

  /*** Call Newton-RK4 method ***/
  cout << setprecision(8) ;
  iter=shootRK4(N,a,b,alpha,beta,t,maxIter,tol,x,y) ;

  /*** Print results to screen ***/
  cout << setprecision(4) ;
  cout << "Number of Newton iterations: " << iter << endl ;
  cout << "Approx solution: t = " << t << endl ;
  cout << "Approx solution: x_j, y_j =  " << endl ;
  for(int j=0; j<N+1; j++){
    cout << "x =" << setw(6) << x(j) ;
    cout << "   " ;
    cout << "y =" << setw(10) << y(j) ;
    cout << "   " ; 
    cout << endl;
  }

  return 0 ; //terminate main program
}

