#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
gcc ke_dbl_f.c -lm -O2 -fPIC -pie -o ke_dbl_f.so -Wl,-E  &&  ./ke_dbl_f.so 2. 1. 10
*/

#define Nmax 89
#define K 0x1.799b34c7fac93p-2 

// atan(2**-k)
static double f_a[Nmax] = {
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
0x1.0000000000000p-60, 0x1.0000000000000p-61, 
};

// shift values (double iterations for k<27)
static short k_n[Nmax] = {
 0,  0,  1,  1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,
 8,  8,  9,  9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15,
16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23,
24, 24, 25, 25, 26, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52,
53, 54, 55, 56, 57, 58, 59, 60, 61,
};


double f_Ecs(double M, double e, double *E, double *ecosE, double *esinE, int N){
    // floating point 45 ms
    // convert 1/2**n into a double multiplier
// f64_Ecs_mulcpsgn,
    double ec=K*e, es=0, ect, En=0, t;
    long long mul;
    short k, n;
    for (n=0; n<N; n++){
        k = k_n[n];
        t = (En-es) - M;
        En -= copysign(f_a[k], t);
        mul = ((1LL<<(11-1)) - 1 - k) << 52;
        t = copysign(*(double *) &mul, t);
//        printf("%d %lld %g %g %g\n", n, mul, t, t, En);
        ect = ec;
        ec += es * t;
        es -= ect * t;
    }
    *ecosE = ec;
    *esinE = es;
    *E = En;
    return En;
}



int v_Efunc(double (*func)(), double *M, double e, double *E, double *ecosE, double *esinE, long n, int N) {
   // vectorised version for *_Ecs
   // cannot be vectorised by gcc (?) https://gcc.gnu.org/projects/tree-ssa/vectorization.html
   while (n--) {
      (*func)(*M++, e, E++, ecosE++, esinE++, N);
   }
   return 0;
}

int main(int argc, char *argv[]) {
    int N;
    double x, e;
    x = atof(argv[1]);
    e = atof(argv[2]);
    N = atoi(argv[3]);

    double u, hh, M=x-e*sin(x);
    printf("E:               %.20g\n", x);
    printf("i64_Ecs          %.20g\n", f_Ecs(M, e, &hh, &hh, &hh, N));

    printf("v_Ecs       %.20g\n", u);

    return 0;
}

