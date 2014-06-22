#ifndef LINEAR_CENTERD_DIFFERENCE_2D_HPP
#define LINEAR_CENTERD_DIFFERENCE_2D_HPP

void (*PDEeval)(const double&, const double&, double&, double&, double&, 
		double&, double&, double&);
void setPDEeval(void (*)(const double&, const double&, double&, double&, double&, 
		     double&, double&, double&));

void (*BCeval)(const double&, const double&, double&, double&, double&, double&);
void setBCeval(void (*)(const double&, const double&, double&, double&, double&,
			 double&));


void AGeval(int N, int M, double a, double b, double c, double d, vector& x,
	    vector& y, matrix& A, vector& G);

int linearcd2D(int N, int M, double a, double b, double c, double d,
	       vector& x, vector& y, matrix& u);

#endif // LINEAR_CENTERD_DIFFERENCE_2D_HPP
