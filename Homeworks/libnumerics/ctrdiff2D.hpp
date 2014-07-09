#ifndef CENTERED_DIFFERENCE_2D_HPP
#define CENTERED_DIFFERENCE_2D_HPP

void (*ctr_PDEeval)(const double&, const double&, double&,
                double&, double&, double&, double&, 
                double&, double&, double&);
void ctr_setPDEeval(void (*)(const double&, const double&, double&,
                         double&, double&, double&, double&, 
                         double&, double&, double&));

void (*ctr_BCeval)(const double&, const double&, double&,
               double&, double&, double&, double&);
void ctr_setBCeval(void (*)(const double&, const double&, double&,
                        double&, double&, double&, double&));

void (*ctr_ICeval)(const double&, const double&, double&, double&);
void ctr_setICeval(void (*)(const double&, const double&, double&, double&));

int ctrdiff2D(int N, int M, 
              double a, double b, double c, double d, 
              vector& x, vector& y, double& t, double& dt, 
              matrix& uold, matrix& u);

#endif // CENTERED_DIFFERENCE_2D_HPP
