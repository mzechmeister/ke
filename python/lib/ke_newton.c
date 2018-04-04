#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>

/*
   gcc -c  -Wall -O2 -fPIC ke_newton.c; gcc -o ke_newton.so -shared ke_newton.o
*/


static const double 
half =  5.00000000000000000000e-01, /* 0x3FE00000, 0x00000000 */
S1  = -1.66666666666666324348e-01, /* 0xBFC55555, 0x55555549 */
S2  =  8.33333333332248946124e-03, /* 0x3F811111, 0x1110F8A6 */
S3  = -1.98412698298579493134e-04, /* 0xBF2A01A0, 0x19C161D5 */
S4  =  2.75573137070700676789e-06, /* 0x3EC71DE3, 0x57B1FE7D */
S5  = -2.50507602534068634195e-08, /* 0xBE5AE5E6, 0x8A2B9CEB */
S6  =  1.58969099521155010221e-10; /* 0x3DE5D93A, 0x5ACFD57C */

static const double 
one =  1.00000000000000000000e+00, /* 0x3FF00000, 0x00000000 */
C1  =  4.16666666666666019037e-02, /* 0x3FA55555, 0x5555554C */
C2  = -1.38888888888741095749e-03, /* 0xBF56C16C, 0x16C15177 */
C3  =  2.48015872894767294178e-05, /* 0x3EFA01A0, 0x19CB1590 */
C4  = -2.75573143513906633035e-07, /* 0xBE927E4F, 0x809C52AD */
C5  =  2.08757232129817482790e-09, /* 0x3E21EE9E, 0xBDB4B1C4 */
C6  = -1.13596475577881948265e-11; /* 0xBDA8FAE9, 0xBE8838D4 */

/* simple version of sin and cos for speed comparison */

double my_sin(double x) {
   long n; 
   double z;
   n = round(1/M_PI*x);
   x -= n*M_PI;
   if (n&1) x = -x;
   z = x * x;
   return x+ x*z*(S1+z*(S2+z*(S3+z*(S4+z*(S5+z*S6)))));
}

double my_cos(double x) {
   long n; 
   double z, cosx;
   n = round(1/M_PI*x);
   x -= n*M_PI;
   z = x * x;
   cosx = one - z*(0.5 - z*(C1+z*(C2+z*(C3+z*(C4+z*(C5+z*C6))))));
   return (n&1)? -cosx:cosx;
}

double _E0(double M, double e) {
   /* start guess for E0 for range reduction */
   return (M>2*M_PI*round(1/(2*M_PI)*M))? M+0.85*e : M-0.85*e;
}

double _E(double M, double e, double eps) {
   /* Danby start guess */
   double dE, En=_E0(M,e);
   int nmax=20;
   do {
      En -= (dE = (En-e*sin(En)-M) / (1-e*cos(En)));
   } while (fabs(dE)>eps && --nmax);
   return En;
}

double _E_my(double M, double e, double eps) {
   /* Danby start guess*/
   double dE, En=_E0(M,e);
   int nmax=20;
   do {
      dE = (En-e*my_sin(En)-M) / (1-e*my_cos(En));
      En -= dE;
   } while (fabs(dE)>eps && --nmax);
   return En;
}

double _E_N(double M, double e, int N) {
   /* Danby start guess */
   // fixed number of iterations to check how many cordic rotations this corresponds to.
   double dE, En=_E0(M,e);
   int n;
   for (n=0; n<N; ++n) {
      En -= (dE = (En-e*sin(En)-M) / (1-e*cos(En)));
   };
   return En;
}

double _E_myN(double M, double e, int N) {
   /* Danby start guess*/
   double dE, En=_E0(M,e);
   int n;
   for (n=0; n<N; ++n) {
      dE = (En-e*my_sin(En)-M) / (1-e*my_cos(En));
      En -= dE;
   };
   return En;
}

int vectorise(double (*func)(), double *M, double *En, double e, double eps, long nM) {
   // vectorised version of _E
   long i;
   for (i=0; i<nM; ++i) {
      En[i] = func(M[i], e, eps);
   }
   return 0;
}

int v_sin(double *x, double *y, long nM) {
   long i;
   for (i=0; i<nM; ++i) {
      y[i] = my_sin(x[i]);
   }
   return 0;
}

int v_cos(double *x, double *y, long nM) {
   long i;
   for (i=0; i<nM; ++i) {
      y[i] = my_cos(x[i]);
   }
   return 0;
}



