#ifndef GAUSS_ELIM_HPP
#define GAUSS_ELIM_HPP

#include "matrix.hpp"

int pivot_index(matrix& A, int k, int n);
void row_exchange(matrix& A, vector& b, int k, int p, int n);
void fwd_elimination(matrix& A, vector& b, int k, int n);
void bwd_substitution(matrix& A, vector& b, vector& x, int n);
int gauss_elim(matrix& A, vector& b, vector& x);

#endif // GAUSS_ELIM_HPP
