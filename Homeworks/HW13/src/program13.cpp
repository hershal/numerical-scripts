/************************************************************
Program 13.  Uses the central-difference method to find an
approximate solution of a hyperbolic IBVP in a rectangular
domain of the form

 utt + s ut = P uxx + Q uyy + p ux + q uy + r u + eta,
                                   a<=x<=b, c<=y<=d, 0<=t<=T

 u(a,y,t) = ga(y,t),   u(b,y,t) = gb(y,t),  c<=y<=d, 0<=t<=T
 u(x,c,t) = gc(x,t),   u(x,d,t) = gd(x,t),  a<=x<=b, 0<=t<=T

 u(x,y,0) = f(x,y),                         a<=x<=b, c<=y<=d
 ut(x,y,0) = gamma(x,y),                    a<=x<=b, c<=y<=d


Inputs:
  PDEeval       Function to evaluate P,Q,p,q,r,s,eta
  BCeval        Function to evaluate ga,gb,gc,gd
  ICeval        Function to evaluate f,gamma
  a,b,c,d       Space domain parameters
  N,M           Number of interior x,y pts (N+2,M+2 total pts)
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  T             Time interval parameter (final time)
  L             Number of time steps

Outputs:
  x             Grid point vector: x(i)=a+i*dx, i=0...N+1
  y             Grid point vector: y(j)=c+j*dy, j=0...M+1
  u             Approx soln at time t=T: u(i,j) is soln
                at x(i),y(j), i=0...N+1, j=0...M+1


Note 1: The function file ctrdiff2D.cpp is incomplete; you'll
need to finish coding the method as indicated in that file.

Note 2: For any given problem, the functions PDEeval, BCeval
and ICeval must be changed.

Note 3: For any given problem, the grid parameters a,b,c,d,T
and N,M,L must be specified.

Note 4: To compile this program use the command (all on one
line)

  c++ -o program13  matrix.cpp  ctrdiff2D.cpp  program13.cpp

Note 5: The program output is written to a file.
************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#include "matrix.hpp"
#include "ctrdiff2D.hpp"

using namespace std;

/*** Define output file ***/
const char myfile[20]="program13.out" ;
ofstream prt(myfile) ;

/*** Define P(x,y,t), Q(x,y,t),
     p(x,y,t), q(x,y,t), r(x,y,t), s(x,y,t), eta(x,y,t) ***/
void ctr_PDEeval_override(const double& x, const double& y, double& t,
                          double& P, double& Q, double& p, double& q,
                          double& r, double& s, double& eta) {

    // P = 0.1*(1 + exp(-x*y*t/10)) ;
    // Q = 0.1 ;
    // p = 0.02*x*y ;
    // q = 0.04*x ;
    // r = 0.01 ;
    // s = 0.05 ;
    // eta = 0.01*t - 0.02*x*sin(y) ;

    // For Programming minilab
    double pi=4.0*atan(1.0);
    double xs = 1.5;
    double ys = 0.0;
    double xx = -1;
    double yy = 0;
    double tt = 0.5;
    double w = sqrt(pow(x-xs,2) + 0.25*pow(y-ys,2));
    eta = 10.0*(exp(-30*pow(x-xx,2)
                    - 30*pow(y-yy,2)
                    - 30*pow(t-tt,2)));
    if(w>=0 && w<=0.5) {
        s = 10*pow(cos(pi*w),2);
    } else {
        s = 0;
    }
    P = Q = 0.3;
    p = q = r = 0;
}

/*** Define ga(y,t), gb(y,t), gc(x,t), gd(x,t) ***/
void ctr_BCeval_override(const double& x, const double& y, double& t,
                         double& ga, double& gb, double& gc, double& gd) {

    // ga = y ;
    // gb = sqrt(y) ;
    // gc = 0 ;
    // gd = x ;

    ga = gb = gc = gd = 0;
}

/*** Define f(x,y), gamma(x,y) ***/
void ctr_ICeval_override(const double& x, const double& y,
                         double& f, double& gamma) {

    // f = 1 ;
    // gamma = x*y*(1-x)*(1-y) ;
    f = gamma = 0;
}

int main() {
    /*** Define problem parameters ***/
    ctr_setPDEeval(&ctr_PDEeval_override);
    ctr_setBCeval(&ctr_BCeval_override);
    ctr_setICeval(&ctr_ICeval_override);

    // int N=9, M=9, L=50, success_flag=0 ;
    // matrix u(N+2,M+2), uold(N+2,M+2) ;
    // vector x(N+2), y(M+2) ;
    // double a=0, b=1, c=0, d=1 ;
    // double dx=(b-a)/(N+1), dy=(d-c)/(M+1) ;
    // double T=5, dt=T/L ;

    //  Minilab Part B
    // // For \Delta t = \frac{1}{10}
    // int L = 70;
    // double T = 7;

    // // For \Delta t = \frac{1}{15}
    // int L = 105;
    // double T = 7;

    // // For \Delta t = \frac{1}{20}
    // int L = 140;
    // double T = 7;

    // // For \Delta t = \frac{1}{25}
    // int L = 175;
    // double T = 7;

    // // For \Delta t = \frac{1}{30}
    // int L = 210;
    // double T = 7;

    //  Minilab Part C
    // // For [0,T] = [0, 3.5]
    // double T = 3.5;
    // int L = 105;

    // // For [0,T] = [0, 5]
    // double T = 5;
    // int L = 150;

    // For [0,T] = [0, 7]
    double T = 7;
    int L = 210;

    int N=119 , M=119 , success_flag=0;
    matrix u(N+2,M+2), uold(N+2,M+2) ;
    vector x(N+2), y(M+2) ;
    double a=-3, b=3, c=-3, d=3;
    double dx=(b-a)/(N+1), dy=(d-c)/(M+1) ;
    double dt=T/L ;

    /*** Construct xy-grid ***/
    for(int i=0; i<=N+1; i++){
        x(i) = a + i*dx ;
    }
    for(int j=0; j<=M+1; j++){
        y(j) = c + j*dy ;
    }

    /*** Load initial values for u ***/
    double t=0 ;
    double gLeft, gRight, gBottom, gTop, f, gamma ;
    for(int i=0; i<=N+1; i++){ //actual i-index on grid
        for(int j=0; j<=M+1; j++){ //actual j-index on grid
            ctr_BCeval_override(x(i),y(j),t,gLeft,gRight,gBottom,gTop) ;
            ctr_ICeval_override(x(i),y(j),f,gamma) ;
            if( j==M+1 ){ u(i,j) = gTop ; }
            if( i==N+1 ){ u(i,j) = gRight ; }
            if( i==0 ){ u(i,j) = gLeft ; }
            if( j==0 ){ u(i,j) = gBottom ; }
            if((i>0)&&(i<N+1)&&(j>0)&&(j<M+1)){ u(i,j) = f ; }
        }
    }
    uold = u ; //initialize uold

    /*** Call ctr-diff method at each time step (overwrites u,uold) ***/
    for(int n=0; n<L; n++){
        success_flag = ctrdiff2D(N,M,a,b,c,d,x,y,t,dt,uold,u) ;
        t = t + dt ;
    }

    /*** Print results at final time to output file ***/
    prt.setf(ios::fixed) ;
    prt << setprecision(5) ;
    cout << "Ctr-Diff-2D: output written to " << myfile << endl ;
    prt << "Ctr-Diff-2D results" << endl ;
    prt << "Number of interior x-grid pts: N = " << N << endl ;
    prt << "Number of interior y-grid pts: M = " << M << endl ;
    prt << "Number of time steps: L = " << L << endl ;
    prt << "Final time of simulation: t = " << t << endl ;
    prt << "Approximate solution at time t: x_i, y_j, u_ij" << endl ;
    for(int i=0; i<=N+1; i++){
        for(int j=0; j<=M+1; j++){
            prt << setw(8) << x(i) ;
            prt << "   " ;
            prt << setw(8) << y(j) ;
            prt << "   " ;
            prt << setw(8) << u(i,j) ;
            prt << endl;
        }
    }

    return 0 ; //terminate main program
}
