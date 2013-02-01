/*******************************************************************
Function to solve Ax=b using the SOR method
Gradient method.  A is assumed to be symmetric, pos-definite.

Inputs: matrix Cinv(n,n), matrix A(n,n), vector b(n), vector x(n),
        integer maxIter, double tol.

Outputs: vector x(n), integer iter (returned)
********************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include <stdio.h>
using namespace std;


int sor(matrix& A, double w, vector& b, vector& x, 
	int maxIter, double tol) {
  
  /*** Check input data ***/
  int n = dim(A,0);

  if(dim(A,1)!=n || dim(b)!=n || dim(x)!=n) {
    cerr << "SOR: bad data -- exiting" << endl; 
    exit(EXIT_FAILURE); 
  }
  if(tol <= 0) {
    cerr << "SOR: bad tol -- exiting" << endl; 
    exit(EXIT_FAILURE);
  }
  if(maxIter <= 0) {
    maxIter = 1;
  }

  /*** Implement algo ***/

  vector r(n);
  int iter = 0;
  double error;
  double sum;

  do {

    for(int i=0; i<n; i++) {

      sum=0;
      for(int j=0; j<n; j++) {
	if (i!=j) {
	  sum += A(i,j)*x(j);
	}
      }
      x(i) = (1-w)*x(i) + ((w/A(i,i))*(b(i)-sum));
    }

    r = b - matVecMult(A,x);
    error = vecMaxNorm(r);

    iter++;
  } while(iter<maxIter && error>=tol);

  /*** Print exit message to screen ***/
  if(error<tol) {
    cout << "SOR: solution converged, |r|_inf = " << error << endl;
  }
  else {
    cout << "SOR: max iterations exceeded" << endl;
  }
  return iter;
}

