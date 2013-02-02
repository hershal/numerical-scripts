/****************************************************
Program 1.  Solves Ax=b using Jacobi method.  

Inputs: A, b, x^(0), maxIter, tol
Outputs: x^(k), iteration count k

Here's how to get started:

1) Copy program1.cpp (this file), jacobi.cpp, 
matrix.cpp and matrix.h into your working directory.

2) Compile (and link) the files by typing
"c++ -o program1 matrix.cpp jacobi.cpp program1.cpp"
at the Linux prompt.

3) Type "program1" to run the program.

4) For any given problem, you'll need to set the 
values of the variables {n,A,b,x,maxIter,tol} below, 
and then compile and run as described above.
****************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;

/*** Declare user-defined functions to be used ***/
int jacobi(matrix&, vector&, vector&, int, double);
int gauss_seidel(matrix&, vector&, vector&, int, double);


/*** Main program ***/
int main() {

  /*** Define and input problem data ***/
  int n=4, maxIter=25, iterGauss, iterJacobi;  
  double tol=1e-3;
  matrix A(n,n);
  vector x(n), xGauss(n), xJacobi(n), b(n);

  xGauss = xJacobi = 0;

  for(int i=0; i<dim(xGauss); i++) {
    xGauss = xJacobi = 88;
  }

  A(0,0)= 4; A(0,1)= 1; A(0,2)=-1; A(0,3)= 1; b(0)=-2;
  A(1,0)= 1; A(1,1)= 4; A(1,2)=-1; A(1,3)=-1; b(1)=-1;
  A(2,0)=-1; A(2,1)=-1; A(2,2)= 5; A(2,3)= 1; b(2)= 0;
  A(3,0)= 1; A(3,1)=-1; A(3,2)= 1; A(3,3)= 3; b(3)= 1;


  /*** Print data to screen ***/
  cout << endl; 
  cout << "Given: A = " << endl;
  cout << A << endl;
  cout << "Given: b = " << endl;
  cout << b << endl;
  cout << "Given: x^(0) = " << endl;
  cout << xGauss << endl;


  /*** Call Jacobi function ***/
  iterJacobi=jacobi(A,b,xJacobi,maxIter,tol);
  

  /*** Print results to screen ***/
  cout << endl; 
  cout << "Iteration index: k = " << iterJacobi << endl;
  cout << endl; 
  cout << "Approximate solution: x^(k) = " << endl;
  cout << xJacobi << endl;

  /*** Call Jacobi function ***/
  iterGauss=gauss_seidel(A,b,xGauss,maxIter,tol);
  
  /*** Print results to screen ***/
  cout << endl; 
  cout << "Iteration index: k = " << iterGauss << endl;
  cout << endl; 
  cout << "Approximate solution: x^(k) = " << endl;
  cout << xGauss << endl;

  
  return 0;
}

