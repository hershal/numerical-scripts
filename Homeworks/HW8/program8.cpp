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


/*** Declare function with steepest descent algorithm ***/
int descent(vector&, int, int, double, int, double) ;


/*** Define g function for problem ***/
void geval(vector& x, double& g){
  // double a=1, b=-0.1 ; //define any constants
  // double pi=4.0*atan(1.0) ; //the number pi
  // g = a + b*x(0)*x(1) + x(0)*x(0) + x(1)*x(1) ;

  double a=1.0, b=1.5, L=0.7;
  // vector ref(2);
  // ref(0)=1; ref(0)=0;
  // double angle = arcos(vecDot(ref,x)/(vecL2Norm(ref)*vecL2Norm(x)));


  // double r = sqrt(pow(x(0),2)+pow(x(1),2));
  // double angle = atan2(x(1),x(0));

  double r = x(0);
  double angle = x(1);
  
  double y=pow(L*cos(angle),2)+pow(r-L*sin(angle),2);
  double z=pow(L*cos(angle),2)+pow(r+L*sin(angle),2);

  g=(pow(a,2)/pow(y,2)-(2*a)/y)+(pow(b,2)/pow(z,2)-(2*b)/z);
}


/*** Define dg function (gradient) for problem ***/
void dgeval(vector& x, vector& dg){
  // double a=1, b=-0.1 ; //define any constants
  // double pi=4.0*atan(1.0) ; //the number pi
  // dg(0) = b*x(1) + 2*x(0) ;
  // dg(1) = b*x(0) + 2*x(1) ;


  // double a=1.0, b=1.5, L=0.7;
  double tol = pow(10,-5);

  vector tmp(2);
  double left, right;

  // double r = sqrt(pow(x(0),2)+pow(x(1),2));
  // double angle = atan2(x(1),x(0));

  double r = x(0);
  double angle = x(1);

  double rdenom = (49+100*pow(r,2)-140*r*sin(angle));

  dg(0) = 2000*((200*(-10*r+7*sin(angle)))/pow(rdenom,3)
		+(20*(10*r- 7*sin(angle)))/pow(rdenom,3)
		-(450*(10*r-7*sin(angle)))/pow(rdenom,3)
		+(3*(10*r+7*sin(angle)))/pow(rdenom,2));

  dg(1) = (7/10)*r*cos(angle)*(-(40000)/(pow(rdenom,2))
			     +(4/pow(((49/100)+pow(r,2)-7/5*r*sin(angle)),3))
			     -9000000/(pow(rdenom,3))
			     +60000/pow(rdenom,2));
  

  // tmp(0) = tol; tmp(1) = 0;
  // tmp = tmp+x;
  // geval(tmp,left);
  // geval(x,right);
  // dg(0) = (left-right)/tol;
  // 
  // tmp(0) = 0; tmp(1) = tol;
  // tmp = tmp+x;
  // geval(tmp,left);
  // geval(x,right);
  // dg(1) = (left-right)/tol;
}


int main() {
  /*** Define problem parameters ***/
  int n=2, maxIter=100, iter=0, maxSrch=20 ;
  double tol=1e-6, a0=0.2, g ;
  vector x(n) ;
  x(0) = 1.0 ; x(1) = 0.3 ; //initial guess

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

