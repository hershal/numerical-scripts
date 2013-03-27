/*******************************************************
Program 8.  Uses the steepest descent method to find a 
local minimum of a function g(x).

Inputs:  
  geval    Function to evaluate g(x), scalar
  dgeval   Function to evaluate dg/dx(x), n vector
  x        Initial guess for x, n vector  
  maxIter  Maximum number of descent iterations
  tol      Tolerance parameter for ||dg/dx||

  a0       Initial guess for alpha (step size)
  maxSrch  Maximum number of alpha search steps


Outputs: 
  x        Approx local minimum point of g(x)
  g        Approx value of g(x) at local min
  k        Number of descent steps taken


Note 1: For any given problem the functions geval and
dgeval must be changed.  These functions compute the
scalar g and vector dg for any given x.

Note 2: The function file descent.cpp is incomplete; 
you'll need to code the quadratic interpolation step
as indicated in that file.

Note 3: To compile this program use the command (all 
on one line)

  c++ -o program8 matrix.cpp descent.cpp program8.cpp

*******************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std ;

#define a 1.0
#define b 1.5
#define L 0.7

/*** Declare function with steepest descent algorithm ***/
int descent(vector&, int, int, double, int, double) ;
double y(double rr, double aa);
double z(double rr, double aa);
double dzdr(double rr, double aa);
double dzda(double rr, double aa);
double dydr(double rr, double aa);
double dyda(double rr, double aa);


/*** Define g function for problem ***/
void geval(vector& x, double& g){

  double r = x(0);
  double angle = x(1);
  
  double y=pow(L*cos(angle),2)+pow(r-L*sin(angle),2);
  double z=pow(L*cos(angle),2)+pow(r+L*sin(angle),2);

  g=(pow(a,2)/pow(y,2)-(2*a)/y)+(pow(b,2)/pow(z,2)-(2*b)/z);
}


/*** Define dg function (gradient) for problem ***/
void dgeval(vector& x, vector& dg){

  double r = x(0);
  double angle = x(1);

  dg(0) = dydr(r,angle)*(-2.0*pow(a,2)/pow(y(r,angle),3) + 2.0*a/pow(y(r,angle),2) + dzdr(r,angle)*(-2.0*pow(b,2)/pow(z(r,angle),3)+2.0*b/pow(z(r,angle),2)));
  dg(1) = dyda(r,angle)*(-2.0*pow(a,2)/pow(y(r,angle),3) + 2.0*a/pow(y(r,angle),2) + dzda(r,angle)*(-2.0*pow(b,2)/pow(z(r,angle),3)+2.0*b/pow(z(r,angle),2)));

}

double y(double rr, double aa) {
  return pow(L*cos(aa),2)+pow(rr-L*sin(aa),2);
}

double z(double rr, double aa) {
  return pow(L*cos(aa),2)+pow(rr+L*sin(aa),2);
}

double dzdr(double rr, double aa) {
  return 2.0*(rr+L*sin(aa));
}

double dzda(double rr, double aa) {
  return -2.0*(L*sin(aa)*L*cos(aa)*(rr*L*sin(aa)));
}

double dydr(double rr, double aa) {
  return 2.0*(rr-L*sin(aa));
}

double dyda(double rr, double aa) {
  return -2.0*(L*sin(aa)*L*cos(aa)+((rr+L*sin(aa)*L*cos(aa))));
}


int main() {
  /*** Define problem parameters ***/
  int n=2, maxIter=100, iter=0, maxSrch=20 ;
  double tol=1e-6, a0=0.3, g ;
  vector x(n) ;
  x(0) = 1.0 ; x(1) = 0.0 ; //initial guess

  /*** Print problem info ***/
  cout << setprecision(10) ;
  cout << endl ;
  cout << "System size: n = " << n << endl ;
  cout << "Initial guess: x^(0) = " << endl ;
  cout << x << endl ;

  /*** Call steepest descent w/initial guess***/
  iter = descent(x,n,maxIter,tol,maxSrch,a0) ;

  /*** Print results to screen ***/
  geval(x,g) ;
  cout << endl ;
  cout << "Iteration index: k = " << iter << endl ;
  cout << endl ;
  cout << "Approx solution: x = " << endl ;
  cout << x << endl ;
  cout << "Approx g-value: g = " << g << endl ;

  return 0; //terminate main program
}

