/*******************************************************************
Function to find dominant eigenvalue and eigenvector of a square 
matrix A using the symmetric power method.  The eigenvalue and 
eigenvector are assumed to be real, and the matrix A is assumed 
to be diagonalizable.

Inputs: matrix A(n,n), vector x(n), double lambda,
        integer n, integer maxIter, double tol.

Outputs: vector x(n), double lambda, integer iter (returned)
********************************************************************/
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
using namespace std;


int sympower(matrix& A, vector& x, double& lambda, 
                          int n, int maxIter, double tol) {
  
  int iter=0, pindex ;
  double error=1, ynorm ;
  vector y(n), r(n) ;

  y = x ;
  while(iter<maxIter && error>=tol) {

    ynorm = vecL2Norm(y) ;
    x = scaleVec(1/ynorm,y) ;

    y = matVecMult(A,x) ;
    lambda = vecDot(x,y);
    r = scaleVec(1/lambda,y) - x ;
    error = vecL2Norm(r) ;
    iter++ ;
  }

  if(error < tol) {
    cout << "SymPwr: soln converged, |r|_inf = " << error << endl;
  }
  else {
    cout << "SymPwr: max iterations exceeded" << endl;
  }
  return iter;
}


