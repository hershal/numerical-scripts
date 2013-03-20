/*******************************************************
Program 7.  Uses the Newton method to solve a general
system of equations
                        F(x) = 0.

Inputs:  
  Feval    Function to evaluate F(x), n vector 
  DFeval   Function to evaluate dF/dx(x), nxn matrix
  x        Initial guess, n vector  
  maxIter  Maximum number of iterations
  tol      Tolerance parameter

Outputs: 
  x        Approx soln of F(x)=0, n vector
  iter     Number of Newton steps taken

Note 1: For any given problem the functions Feval and
DFeval must be changed.  These functions compute the
vector F and matrix DF for any given x.

Note 2: To compile this program use the command 
(all on one line)

c++ -o program7 matrix.cpp gauss_elim.cpp 
                                newton.cpp program7.cpp

*******************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std ;


/*** Declare function with Newton algorithm ***/
int newton(vector&, int, int, double) ;


/*** Define F function for problem ***/
void Feval(vector& x, vector& F){
  double a=3, b=1.5, c=1, d=2, h=3.5, r=2.5;
  double pi=4.0*atan(1.0) ; //the number pi

  // For 2b
  // F(0) = pow(x(0),2)+x(1)-37;
  // F(1) = x(0)-pow(x(1),2)-5;
  // F(2) = x(0)+x(1)+x(2)-3;

  // For 7b
  // F(0) = log(pow(x(0),2)+pow(x(1),2))-sin(x(0)*x(1))-log(2)-log(pi);
  // F(1) = exp(x(0)-x(1))+cos(x(0)*x(1));

  // For the programming assignment
  F(0) =
  F(1) = 
  F(2) = 
  F(3) = 
  
}

/*** Define DF function (Jacobian) for problem ***/
void DFeval(vector& x, matrix& DF){
  double a=3, b=1.5, c=1, d=2, h=3.5, r=2.5;
  double pi=4.0*atan(1.0) ; //the number pi

  // For 2b
  // DF(0,0) = 2*x(0);
  // DF(1,0) = 1;
  // DF(2,0) = 0;
  // DF(0,1) = 1;
  // DF(1,1) = -2*x(1);
  // DF(2,1) = 0;
  // DF(0,2) = 1;
  // DF(1,2) = 1;
  // DF(2,2) = 1;

  // For 7b
  // DF(0,0) = (2*x(0))/(pow(x(0),2)+pow(x(1),2)) - x(1)*cos(x(0)*x(1));
  // DF(0,1) = (2*x(1))/(pow(x(0),2)+pow(x(1),2)) - x(0)*cos(x(0)*x(1));
  // DF(1,0) =  exp(x(0) - x(1)) - x(1)*sin(x(0)*x(1));
  // DF(1,1) = -exp(x(0) - x(1)) - x(0)*sin(x(0)*x(1));

}

int main() {
  /*** Define problem parameters ***/
  int n=2, maxIter=100, iter=0 ;  
  double tol=1e-6 ;
  vector x(n) ;
 
  // For 2b
  // x(0) = x(1) = x(2) = 0;

  // for 7b
  x(0) = 2; x(1) = 2;

  /*** Print problem info ***/
  cout << setprecision(10) ;
  cout << endl ;
  cout << "System size: n = " << n << endl ;
  cout << "Initial guess: x^(0) = " << endl ;
  cout << x << endl ;

  /*** Call Newton w/initial guess***/
  iter = newton(x,n,maxIter,tol) ;

  /*** Print results to screen ***/
  cout << endl ;
  cout << "Iteration index: k = " << iter << endl ;
  cout << "Approx solution: x = " << endl ;
  cout << x << endl ;
 
  return 0; //terminate main program
}

