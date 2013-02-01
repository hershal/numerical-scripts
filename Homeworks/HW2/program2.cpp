/*******************************************************
Program 2.  Solves Ax=b using Conjugate Gradient method
with pre-conditioning matrix Cinv.   A is assumed to be
symmetric and positive-definite.  

Inputs: A, b, x^(0), Cinv, maxIter, tol
Outputs: x^(k), iteration count k

Here's how to get started:

1) Copy program2.cpp (this file), conjgrad.cpp, 
matrix.cpp and matrix.h into your working directory.

2) Compile (and link) the files by typing
"c++ -o program2 matrix.cpp conjgrad.cpp program2.cpp"
at the Linux prompt.

3) Type "program2" to run the program.

4) For any given problem, you'll need to set the 
values of the variables {n,A,b,x,Cinv,maxIter,tol} 
below, and then compile and run as described above.

5) If you add another function for making comparisons, 
e.g. jacobi, gauss_seidel or sor, remember to include 
a header for the function in the main program below, 
and to include the filename when compiling as in 
Step 2 above.  Also, remember to re-initialize x before 
running another method.
*******************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


/*** Declare user-defined functions to be used ***/
int conjgrad(matrix&, matrix&, vector&, vector&, int, double);
int sor(matrix&, double, vector&, vector&, int, double);

/*** Main program ***/
int main() {

  /*** Define and input problem data ***/
  int n=100, maxIter=5000, iter, iterSor;
  double tol=1e-5, w=1.6;
  matrix Cinv(n,n), A(n,n);
  vector xCG(n), xSor(n), b(n);

  // A(0,0)= 3; A(0,1)=-1; A(0,2)= 1; b(0)=1;
  // A(1,0)=-1; A(1,1)= 6; A(1,2)= 2; b(1)=0;
  // A(2,0)= 1; A(2,1)= 2; A(2,2)= 7; b(2)=4;

  // A matrix loop-builder
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      if (j==i-1 || j==i+1) {
  	A(i,j) = -1;
      } else if(j==i) {
  	A(i,j) = 2+i/10.0;
      } 
    }
    b(i) = 1+i/20.0;
  }
  
  
  xSor=xCG=0; //initialize  x

  /*** Construct pre-conditioning matrix ***/
  Cinv=0;
  for(int i=0; i<n; i++) {
    if(A(i,i) <= 0) {
      cerr << "CG: bad data -- exiting" << endl;
      exit(EXIT_FAILURE);
    } 
    else { 
      // Cinv(i,i) = 1/sqrt(A(i,i)); 
      Cinv(i,i) = 1; //Use this for no pre-cond
    }
  }

  /*** Call conjugate gradient function ***/
  iter=conjgrad(Cinv,A,b,xCG,maxIter,tol);
  
  /*** Print results to screen ***/
  cout << endl; 
  cout << "Iteration index: k = " << iter << endl;
  //  cout << "Approximate solution: x^(k) = " << endl;
  //  cout << xCG << endl;

  /*** Call conjugate gradient function ***/
  iterSor=sor(A,w,b,xSor,maxIter,tol);
  
  /*** Print results to screen ***/
  cout << endl; 
  cout << "w values: w = " << w << endl;
  cout << "Iteration index: k = " << iterSor << endl;
  //  cout << "Approximate solution: x^(k) = " << endl;
  //  cout << xSor << endl;

  return 0;
}

