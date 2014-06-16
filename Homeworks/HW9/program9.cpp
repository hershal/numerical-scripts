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

double p(double x,double y);
double q(double x, double y);
double u(double x, double y);
double v(double x, double y);
double w(double x, double y);

/*** Define f(x,y,yp) function ***/
void feval(const double& x, const double& y,
                              const double& yp, double& f){
  double pi=4.0*atan(1.0) ; //the number pi
  f = ((p(x,y)*yp - q(x,y))*(u(x,y) + 2*v(x,y)*yp + w(x,y)*yp*yp))/(1+pow(p(x,y),2) + pow(q(x,y),2));
}

double p(double x,double y) {
  return (-2*x-2)*exp(-pow(x+1,2) - pow(y+1,2)) + 0.5*(-2*x+2)*exp(-pow(x-1,2)-pow(y-1,2));
}
double q(double x, double y) {
  return p(y,x);
}
double u(double x, double y) {
  return -2*exp(-pow(x+1,2)-pow(y+1,2)) + pow(-2*x-2,2)*exp(-pow(x+1,2) - pow(y+1,2))
    - exp(-pow(x-1,2)-pow(y-1,2)) + 0.5*pow(-2*x+2,2)*exp(-pow(x-1,2)-pow(y-1,2));
}
double v(double x, double y) {
  return (-2*x-2)*(-2*y-2)*exp(-pow(x+1,2)-pow(y+1,2))
    + 0.5*(-2*x + 2)*(-2*y+2)*exp(-pow(x-1,2)-pow(y-1,2));
}
double w(double x, double y) {
  return u(y,x);
}

int main() {
  /*** Define problem parameters ***/
  double tol=1e-6 ;
  int N=20, maxIter=10, iter ;
  double a=-3.0, b=3.0, alpha=-2.0, beta=2.0, t ;

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

