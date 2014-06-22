#ifndef LINEAR_FINITE_ELEMENT_HPP
#define LINEAR_FINITE_ELEMENT_HPP

void (*BVPeval)(const double&, double&, double&, double&, double&, double&);
void setBVP(void (*)(const double &, double &, double &, double &, double &, double &));

void PHIeval(int N, vector& x, int i, double xval, double& phi, double& dphi);
void FEMeval(int N, vector& x, matrix& A, vector& F);
int linearfem(int N, vector& x, vector& y);

#endif // LINEAR_FINITE_ELEMENT_HPP
