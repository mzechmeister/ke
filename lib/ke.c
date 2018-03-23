#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
   gcc -c  -Wall -O2 -fPIC ke.c; gcc -o ke.so -shared ke.o
   gcc -c  -Wall -O2 -ansi -pedantic -fPIC ke.c; gcc -o ke.so -shared ke.o; ./ke
   gcc -lm -o kepcordic kepcordic.c; ./kepcordic 
*/

#define TABLE_SIZE 60
static double a[TABLE_SIZE+1] = {
   0x1.921fb54442d18p+0 , 0x1.921fb54442d18p-1 , 0x1.921fb54442d18p-2 , 0x1.921fb54442d18p-3 ,
   0x1.921fb54442d18p-4 , 0x1.921fb54442d18p-5 , 0x1.921fb54442d18p-6 , 0x1.921fb54442d18p-7 ,
   0x1.921fb54442d18p-8 , 0x1.921fb54442d18p-9 , 0x1.921fb54442d18p-10, 0x1.921fb54442d18p-11,
   0x1.921fb54442d18p-12, 0x1.921fb54442d18p-13, 0x1.921fb54442d18p-14, 0x1.921fb54442d18p-15,
   0x1.921fb54442d18p-16, 0x1.921fb54442d18p-17, 0x1.921fb54442d18p-18, 0x1.921fb54442d18p-19,
   0x1.921fb54442d18p-20, 0x1.921fb54442d18p-21, 0x1.921fb54442d18p-22, 0x1.921fb54442d18p-23,
   0x1.921fb54442d18p-24, 0x1.921fb54442d18p-25, 0x1.921fb54442d18p-26, 0x1.921fb54442d18p-27,
   0x1.921fb54442d18p-28, 0x1.921fb54442d18p-29, 0x1.921fb54442d18p-30, 0x1.921fb54442d18p-31,
   0x1.921fb54442d18p-32, 0x1.921fb54442d18p-33, 0x1.921fb54442d18p-34, 0x1.921fb54442d18p-35,
   0x1.921fb54442d18p-36, 0x1.921fb54442d18p-37, 0x1.921fb54442d18p-38, 0x1.921fb54442d18p-39,
   0x1.921fb54442d18p-40, 0x1.921fb54442d18p-41, 0x1.921fb54442d18p-42, 0x1.921fb54442d18p-43,
   0x1.921fb54442d18p-44, 0x1.921fb54442d18p-45, 0x1.921fb54442d18p-46, 0x1.921fb54442d18p-47,
   0x1.921fb54442d18p-48, 0x1.921fb54442d18p-49, 0x1.921fb54442d18p-50, 0x1.921fb54442d18p-51,
   0x1.921fb54442d18p-52, 0x1.921fb54442d18p-53, 0x1.921fb54442d18p-54, 0x1.921fb54442d18p-55,
   0x1.921fb54442d18p-56, 0x1.921fb54442d18p-57, 0x1.921fb54442d18p-58, 0x1.921fb54442d18p-59,
   };
static double cosa[TABLE_SIZE+1] = {
   0x1.1a62633145c07p-54, 0x1.6a09e667f3bcdp-1 , 0x1.d906bcf328d46p-1 , 0x1.f6297cff75cb0p-1 ,
   0x1.fd88da3d12526p-1 , 0x1.ff621e3796d7ep-1 , 0x1.ffd886084cd0dp-1 , 0x1.fff62169b92dbp-1 ,
   0x1.fffd8858e8a92p-1 , 0x1.ffff621621d02p-1 , 0x1.ffffd88586ee6p-1 , 0x1.fffff62161a34p-1 ,
   0x1.fffffd8858675p-1 , 0x1.ffffff621619cp-1 , 0x1.ffffffd885867p-1 , 0x1.fffffff62161ap-1 ,
   0x1.fffffffd88586p-1 , 0x1.ffffffff62162p-1 , 0x1.ffffffffd8858p-1 , 0x1.fffffffff6216p-1 ,
   0x1.fffffffffd886p-1 , 0x1.ffffffffff621p-1 , 0x1.ffffffffffd88p-1 , 0x1.fffffffffff62p-1 ,
   0x1.fffffffffffd9p-1 , 0x1.ffffffffffff6p-1 , 0x1.ffffffffffffep-1 , 0x1.fffffffffffffp-1 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   };
static double sina[TABLE_SIZE+1] = {
   0x1.0000000000000p+0 , 0x1.6a09e667f3bccp-1 , 0x1.87de2a6aea963p-2 , 0x1.8f8b83c69a60ap-3 ,
   0x1.917a6bc29b42cp-4 , 0x1.91f65f10dd814p-5 , 0x1.92155f7a3667ep-6 , 0x1.921d1fcdec784p-7 ,
   0x1.921f0fe670071p-8 , 0x1.921f8becca4bap-9 , 0x1.921faaee6472dp-10, 0x1.921fb2aecb360p-11,
   0x1.921fb49ee4ea6p-12, 0x1.921fb51aeb57bp-13, 0x1.921fb539ecf31p-14, 0x1.921fb541ad59ep-15,
   0x1.921fb5439d73ap-16, 0x1.921fb544197a0p-17, 0x1.921fb544387bap-18, 0x1.921fb544403c1p-19,
   0x1.921fb544422c2p-20, 0x1.921fb54442a83p-21, 0x1.921fb54442c73p-22, 0x1.921fb54442cefp-23,
   0x1.921fb54442d0ep-24, 0x1.921fb54442d15p-25, 0x1.921fb54442d17p-26, 0x1.921fb54442d18p-27,
   0x1.921fb54442d18p-28, 0x1.921fb54442d18p-29, 0x1.921fb54442d18p-30, 0x1.921fb54442d18p-31,
   0x1.921fb54442d18p-32, 0x1.921fb54442d18p-33, 0x1.921fb54442d18p-34, 0x1.921fb54442d18p-35,
   0x1.921fb54442d18p-36, 0x1.921fb54442d18p-37, 0x1.921fb54442d18p-38, 0x1.921fb54442d18p-39,
   0x1.921fb54442d18p-40, 0x1.921fb54442d18p-41, 0x1.921fb54442d18p-42, 0x1.921fb54442d18p-43,
   0x1.921fb54442d18p-44, 0x1.921fb54442d18p-45, 0x1.921fb54442d18p-46, 0x1.921fb54442d18p-47,
   0x1.921fb54442d18p-48, 0x1.921fb54442d18p-49, 0x1.921fb54442d18p-50, 0x1.921fb54442d18p-51,
   0x1.921fb54442d18p-52, 0x1.921fb54442d18p-53, 0x1.921fb54442d18p-54, 0x1.921fb54442d18p-55,
   0x1.921fb54442d18p-56, 0x1.921fb54442d18p-57, 0x1.921fb54442d18p-58, 0x1.921fb54442d18p-59,
   };
static double tanatr[TABLE_SIZE+1] = {
   0x1.0000000000000p+0 , 0x1.0000000000000p-1 , 0x1.0000000000000p-2 , 0x1.0000000000000p-3 ,
   0x1.0000000000000p-4 , 0x1.0000000000000p-5 , 0x1.0000000000000p-6 , 0x1.0000000000000p-7 ,
   0x1.0000000000000p-8 , 0x1.0000000000000p-9 , 0x1.0000000000000p-10, 0x1.0000000000000p-11,
   0x1.0000000000000p-12, 0x1.0000000000000p-13, 0x1.0000000000000p-14, 0x1.0000000000000p-15,
   0x1.0000000000000p-16, 0x1.0000000000000p-17, 0x1.0000000000000p-18, 0x1.0000000000000p-19,
   0x1.0000000000000p-20, 0x1.0000000000000p-21, 0x1.0000000000000p-22, 0x1.0000000000000p-23,
   0x1.0000000000000p-24, 0x1.0000000000000p-25, 0x1.0000000000000p-26, 0x1.0000000000000p-27,
   0x1.0000000000000p-28, 0x1.0000000000000p-29, 0x1.0000000000000p-30, 0x1.0000000000000p-31,
   0x1.0000000000000p-32, 0x1.0000000000000p-33, 0x1.0000000000000p-34, 0x1.0000000000000p-35,
   0x1.0000000000000p-36, 0x1.0000000000000p-37, 0x1.0000000000000p-38, 0x1.0000000000000p-39,
   0x1.0000000000000p-40, 0x1.0000000000000p-41, 0x1.0000000000000p-42, 0x1.0000000000000p-43,
   0x1.0000000000000p-44, 0x1.0000000000000p-45, 0x1.0000000000000p-46, 0x1.0000000000000p-47,
   0x1.0000000000000p-48, 0x1.0000000000000p-49, 0x1.0000000000000p-50, 0x1.0000000000000p-51,
   0x1.0000000000000p-52, 0x1.0000000000000p-53, 0x1.0000000000000p-54, 0x1.0000000000000p-55,
   0x1.0000000000000p-56, 0x1.0000000000000p-57, 0x1.0000000000000p-58, 0x1.0000000000000p-59,
   };
static double atr[TABLE_SIZE+1] = {
   0x1.921fb54442d18p-1 , 0x1.dac670561bb4fp-2 , 0x1.f5b75f92c80ddp-3 , 0x1.fd5ba9aac2f6ep-4 ,
   0x1.ff55bb72cfdeap-5 , 0x1.ffd55bba97625p-6 , 0x1.fff555bbb729bp-7 , 0x1.fffd555bbba97p-8 ,
   0x1.ffff5555bbbb7p-9 , 0x1.ffffd5555bbbcp-10, 0x1.fffff55555bbcp-11, 0x1.fffffd55555bcp-12,
   0x1.ffffff555555cp-13, 0x1.ffffffd555556p-14, 0x1.fffffff555555p-15, 0x1.fffffffd55555p-16,
   0x1.ffffffff55555p-17, 0x1.ffffffffd5555p-18, 0x1.fffffffff5555p-19, 0x1.fffffffffd555p-20,
   0x1.ffffffffff555p-21, 0x1.ffffffffffd55p-22, 0x1.fffffffffff55p-23, 0x1.fffffffffffd5p-24,
   0x1.ffffffffffff5p-25, 0x1.ffffffffffffdp-26, 0x1.fffffffffffffp-27, 0x1.0000000000000p-27,
   0x1.0000000000000p-28, 0x1.0000000000000p-29, 0x1.0000000000000p-30, 0x1.0000000000000p-31,
   0x1.0000000000000p-32, 0x1.0000000000000p-33, 0x1.0000000000000p-34, 0x1.0000000000000p-35,
   0x1.0000000000000p-36, 0x1.0000000000000p-37, 0x1.0000000000000p-38, 0x1.0000000000000p-39,
   0x1.0000000000000p-40, 0x1.0000000000000p-41, 0x1.0000000000000p-42, 0x1.0000000000000p-43,
   0x1.0000000000000p-44, 0x1.0000000000000p-45, 0x1.0000000000000p-46, 0x1.0000000000000p-47,
   0x1.0000000000000p-48, 0x1.0000000000000p-49, 0x1.0000000000000p-50, 0x1.0000000000000p-51,
   0x1.0000000000000p-52, 0x1.0000000000000p-53, 0x1.0000000000000p-54, 0x1.0000000000000p-55,
   0x1.0000000000000p-56, 0x1.0000000000000p-57, 0x1.0000000000000p-58, 0x1.0000000000000p-59,
   };
static double sinatr[TABLE_SIZE+1] = {
   0x1.6a09e667f3bccp-1 , 0x1.c9f25c5bfedd9p-2 , 0x1.f0b6848d2af1cp-3 , 0x1.fc0bd88a0f1d9p-4 ,
   0x1.ff00bf608b827p-5 , 0x1.ffc00bfd808bep-6 , 0x1.fff000bff6009p-7 , 0x1.fffc000bffd80p-8 ,
   0x1.ffff0000bfff6p-9 , 0x1.ffffc0000c000p-10, 0x1.fffff00000c00p-11, 0x1.fffffc00000c0p-12,
   0x1.ffffff000000cp-13, 0x1.ffffffc000001p-14, 0x1.fffffff000000p-15, 0x1.fffffffc00000p-16,
   0x1.ffffffff00000p-17, 0x1.ffffffffc0000p-18, 0x1.fffffffff0000p-19, 0x1.fffffffffc000p-20,
   0x1.ffffffffff000p-21, 0x1.ffffffffffc00p-22, 0x1.fffffffffff00p-23, 0x1.fffffffffffc0p-24,
   0x1.ffffffffffff0p-25, 0x1.ffffffffffffcp-26, 0x1.fffffffffffffp-27, 0x1.0000000000000p-27,
   0x1.0000000000000p-28, 0x1.0000000000000p-29, 0x1.0000000000000p-30, 0x1.0000000000000p-31,
   0x1.0000000000000p-32, 0x1.0000000000000p-33, 0x1.0000000000000p-34, 0x1.0000000000000p-35,
   0x1.0000000000000p-36, 0x1.0000000000000p-37, 0x1.0000000000000p-38, 0x1.0000000000000p-39,
   0x1.0000000000000p-40, 0x1.0000000000000p-41, 0x1.0000000000000p-42, 0x1.0000000000000p-43,
   0x1.0000000000000p-44, 0x1.0000000000000p-45, 0x1.0000000000000p-46, 0x1.0000000000000p-47,
   0x1.0000000000000p-48, 0x1.0000000000000p-49, 0x1.0000000000000p-50, 0x1.0000000000000p-51,
   0x1.0000000000000p-52, 0x1.0000000000000p-53, 0x1.0000000000000p-54, 0x1.0000000000000p-55,
   0x1.0000000000000p-56, 0x1.0000000000000p-57, 0x1.0000000000000p-58, 0x1.0000000000000p-59,
   };
static double cosatr[TABLE_SIZE+1] = {
   0x1.6a09e667f3bcdp-1 , 0x1.c9f25c5bfedd9p-1 , 0x1.f0b6848d2af1cp-1 , 0x1.fc0bd88a0f1d9p-1 ,
   0x1.ff00bf608b827p-1 , 0x1.ffc00bfd808bep-1 , 0x1.fff000bff6009p-1 , 0x1.fffc000bffd80p-1 ,
   0x1.ffff0000bfff6p-1 , 0x1.ffffc0000c000p-1 , 0x1.fffff00000c00p-1 , 0x1.fffffc00000c0p-1 ,
   0x1.ffffff000000cp-1 , 0x1.ffffffc000001p-1 , 0x1.fffffff000000p-1 , 0x1.fffffffc00000p-1 ,
   0x1.ffffffff00000p-1 , 0x1.ffffffffc0000p-1 , 0x1.fffffffff0000p-1 , 0x1.fffffffffc000p-1 ,
   0x1.ffffffffff000p-1 , 0x1.ffffffffffc00p-1 , 0x1.fffffffffff00p-1 , 0x1.fffffffffffc0p-1 ,
   0x1.ffffffffffff0p-1 , 0x1.ffffffffffffcp-1 , 0x1.fffffffffffffp-1 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 , 0x1.0000000000000p+0 ,
   };
static double K[TABLE_SIZE+1] = {
   0x1.0000000000000p+0 , 0x1.6a09e667f3bccp-1 , 0x1.43d136248490ep-1 , 0x1.3a261ba6d7a36p-1 ,
   0x1.37b9141deb3fep-1 , 0x1.371dac182eef5p-1 , 0x1.36f6cfabd961fp-1 , 0x1.36ed1869f27e9p-1 ,
   0x1.36eaaa970b20fp-1 , 0x1.36ea0f222a6d1p-1 , 0x1.36e9e844efd24p-1 , 0x1.36e9de8da104bp-1 ,
   0x1.36e9dc1fcd4eep-1 , 0x1.36e9db8458614p-1 , 0x1.36e9db5d7b25dp-1 , 0x1.36e9db53c3d6fp-1 ,
   0x1.36e9db5156034p-1 , 0x1.36e9db50ba8e5p-1 , 0x1.36e9db5093b11p-1 , 0x1.36e9db5089f9cp-1 ,
   0x1.36e9db50878bfp-1 , 0x1.36e9db5086f08p-1 , 0x1.36e9db5086c9ap-1 , 0x1.36e9db5086bffp-1 ,
   0x1.36e9db5086bd8p-1 , 0x1.36e9db5086bcep-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 ,
   0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 ,
   0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 ,
   0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 ,
   0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 ,
   0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 ,
   0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 ,
   0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 ,
   0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 , 0x1.36e9db5086bccp-1 ,
   };

int n;

double _E0(double M) {
   /* start guess for E0 for range reduction */
   return 2*M_PI * round(1/(2*M_PI)*M);
}

double _E(double M, double e, int N) {
   // half angle rotations, one-sided (zero-one).
   double En=_E0(M), cosE=1, sinE=0, Et, st;
   if (0<(Et=M-En) && Et<M_PI) {
      for (n=0; n<N; ++n) {
         st = cosE*sina[n] + sinE*cosa[n];
         Et = En + a[n];
         if (Et-e*st < M) {
            En = Et;
            cosE = cosE*cosa[n] - sinE*sina[n];
            sinE = st;
         };
      }
   }
   else {
      for (n=0; n<N; ++n) {
         st = -cosE*sina[n] + sinE*cosa[n];
         Et = En - a[n];
         if (Et-e*st > M) {
            En = Et;
            cosE = cosE*cosa[n] + sinE*sina[n];
            sinE = st;
         };
      }
   }
   return En;
}

double _Epn(double M, double e, int N) {
   // Calculate eccentric anomaly cordic-like and with half-angle, two-sided (positive and negative rotations).
   double En=_E0(M), cosE=1, sinE=0, an, sinan, cosEnew;
   for (n=0; n<N; ++n) {
      an = a[n]; //2<<N;
      sinan = sina[n];
      /* printf("%d %f %f %f %f  %f %fb\n", n, En, an,cosE, sinE,  cosa[n], sinan  );*/
      if (En-e*sinE > M) {
         an = -an; 
         sinan = -sinan;
      };
      En += an;
      cosEnew = cosE*cosa[n] - sinE*sinan;
      sinE    = cosE*sinan   + sinE*cosa[n];
      cosE = cosEnew;
   }
   return En;
}

double _Eatr(double M, double e, int N) {
   // Calculate eccentric anomaly cordic-like and with arc tangens radix, two-sided.
   // almost as fast _E, faster than _Epn 
   double En=_E0(M), cosE=0, sinE, cosEnew, tanan;
   sinE = En>M? -1 : 1;
   En += sinE * (M_PI/2);
   for (n=0; n<N; ++n) {
      if (En-e*K[n]*sinE > M) {
         tanan = -tanatr[n];
         En -= atr[n];
      } else {
         tanan = tanatr[n];
         En += atr[n];
      }
      cosEnew = cosE - sinE*tanan;
      sinE += cosE*tanan;
      cosE = cosEnew;
   }
   return En;
}

double _Eatrdiv(double M, double e, int N) {
   // Calculate eccentric anomaly cordic-like and with arc tangens radix, two-sided.
   // test whether division by power of two compiles to a faster program. => No.
   double En=_E0(M), cosE=0, sinE, cosEnew;
   long s;
   sinE = En>M? -1 : 1;
   En += sinE * (M_PI/2);
   for (n=0; n<N; ++n) {
      if (En-e*K[n]*sinE > M) {
         s = -1<<n;
         En -= atr[n];
      } else {
         s = 1<<n;
         En += atr[n];
      }
      cosEnew = cosE - sinE/s;
      sinE += cosE/s;
      cosE = cosEnew;
   }
   return En;
}

double _Eldexp(double M, double e, int N) {
   // Calculate eccentric anomaly cordic-like and with arc tangens radix, two-sided.
   // test whether ldexp can do fast division with power of two. => No.
   double En=_E0(M), cosE=0, sinE, cosEnew;
   sinE = En>M? -1 : 1;
   En += sinE * (M_PI/2);
   for (n=0; n<N; ++n) {
      if (En-e*K[n]*sinE > M) {
         En -= atr[n];
         cosEnew = cosE + ldexp(sinE, -n);
         sinE -= ldexp(cosE, -n);
      } else {
         En += atr[n];
         cosEnew = cosE - ldexp(sinE, -n);
         sinE += ldexp(cosE, -n);
      }
      cosE = cosEnew;
   }
   return En;
}

double _Eatrone(double M, double e, int N) {
   // atr, one-sided (zero-one).
   // same as _E but a => atr
   double En=_E0(M), cosE=0, sinE, Et, st;
   sinE = En>M? -1 : 1;
   En += sinE * (M_PI/2);
   if (0<(Et=M-En) && Et<M_PI) {
      for (n=0; n<N; ++n) {
//         st = cosE*sinatr[n] + sinE*cosatr[n];
         st = (cosE*tanatr[n] + sinE)*cosatr[n];
         Et = En + atr[n];
         if (Et-e*st < M) {
            En = Et;
//             cosE = cosE*cosatr[n] - sinE*sinatr[n];
            cosE = (cosE - sinE*tanatr[n])*cosatr[n];
            sinE = st;
         };
      }
   }
   else {
      for (n=0; n<N; ++n) {
         st = -cosE*sinatr[n] + sinE*cosatr[n];
         Et = En - atr[n];
         if (Et-e*st > M) {
            En = Et;
            cosE = cosE*cosatr[n] + sinE*sinatr[n];
            sinE = st;
         };
      }
   }
   return En;
}

struct Ecossin{
   double E, cosE, sinE;
};

struct Ecossin _Ecs(double M, double e, int N) {
   // half angle zero-one rotations
   // same as _E, but returns the triple
   double En=_E0(M), cosE=1, sinE=0, Et, st;
   struct Ecossin ke;
   if (0<(Et=M-En) && Et<M_PI) {
      for (n=0; n<N; ++n) {
         st = cosE*sina[n] + sinE*cosa[n];
         Et = En + a[n];
         /* printf("%d %f %f %f %f  %f %fb\n", n, En, an,cosE, sinE,  cosa[n], sinan  );*/
         if (Et-e*st < M) {
            En = Et;
            cosE = cosE*cosa[n] - sinE*sina[n];
            sinE = st;
         };
      }
   }
   else {
      for (n=0; n<N; ++n) {
         st = -cosE*sina[n] + sinE*cosa[n];
         Et = En - a[n];
         /* printf("%d %f %f %f %f  %f %fb\n", n, En, an,cosE, sinE,  cosa[n], sinan  );*/
         if (Et-e*st > M) {
            En = Et;
            cosE = cosE*cosa[n] + sinE*sina[n];
            sinE = st;
         };
      }
   }
   ke.E = En;
   ke.cosE = cosE;
   ke.sinE = sinE;
   return ke;
}

double _E1N(double M, double e, int N) {
   // half angle zero-one rotations
   // one Newton iteration 
   double an;
   struct Ecossin ke=_Ecs(M, e, N);

   an = (M-(ke.E-e*ke.sinE)) / (1-e*ke.cosE);
   return ke.E + an;
}

double _EN(double M, double e, int N, double eps) {
   // half angle zero-one rotations
   // multiple iterations
   int nmax=20;
   double dE, cosE;
   struct Ecossin ke=_Ecs(M, e, N);

   do {
      cosE = ke.cosE;
      dE = (M-(ke.E-e*ke.sinE)) / (1-e*cosE);
      ke.E += dE;
      ke.cosE += dE*ke.sinE;
      ke.sinE += dE*cosE;
   } while (fabs(dE)>eps && --nmax);
   return ke.E;
}

int v_Efunc(double (*func)(), double *M, double *En, double e, int N, long nM) {
   // vectorised version of _E
   long i;
//    printf("%d\n", **func);
//   printf("The address of the function is =%p\n",func);
//     printf("The address of the function pointer is =%p\n",&func);
   for (i=0; i<nM; ++i) {
      En[i] = (*func)(M[i], e, N);
   }
   return 0;
}

int v_E(double *M, double *En, double e, int N, long nM) {
   // vectorised version of _E
   long i;
   for (i=0; i<nM; ++i) {
      En[i] = _E(M[i], e, N);
   }
   return 0;
}

int v_Ecs(double *M, double *En, double *cosE, double *sinE,  double e, int N, long nM) {
   // vectorised version of _E
   long i;
   struct Ecossin ke;
   for (i=0; i<nM; ++i) {
      ke = _Ecs(M[i], e, N);
      En[i] = ke.E;
      cosE[i] = ke.cosE;
      sinE[i] = ke.sinE;
   }
   return 0;
}


int main() {
    int (*funcptr)() = main;
    unsigned char *p = (unsigned char *)&funcptr;
    size_t i;
        printf("%p \n", funcptr);
       printf("%p \n", _E);

    for (i = 0; i < sizeof funcptr; i++)
    {
        printf("%02x ", p[i]);
    }
    putchar('\n');

    return 0;
}
