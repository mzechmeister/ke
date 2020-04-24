#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
gcc ke_dbl_i.c -lm -O2 -fPIC -pie -o ke_dbl_i.so -Wl,-E  &&  ./ke_dbl_i.so 2. 1. 56

This is the same as ke_dbl, but loop end is N instead of Nmax
*/

#define R (1LL<<61)
#define dR (1./R)
#define Nmax 89
#define KR 0x1.799b34c7fac93p+59

// atan(2**-k)
static long long i64_a[Nmax] = {
0x1921fb54442d1800, 0xed63382b0dda780, 0x7d6dd7e4b203740, 0x3fab7535585edc0,
0x1ff55bb72cfdea0, 0xffeaaddd4bb128, 0x7ffd556eedca6c, 0x3fffaaab77752e,
0x1ffff5555bbbb7, 0xffffeaaaaddde, 0x7ffffd55556ef, 0x3fffffaaaaab7,
0x1ffffff555555, 0xffffffeaaaaa, 0x7ffffffd5555, 0x3fffffffaaaa,
0x1ffffffff555, 0xffffffffeaa, 0x7ffffffffd5, 0x3fffffffffa,
0x1ffffffffff, 0xffffffffff, 0x7fffffffff, 0x3fffffffff,
0x1fffffffff, 0xfffffffff, 0x7ffffffff, 0x400000000,
0x200000000, 0x100000000, 0x80000000, 0x40000000,
0x20000000, 0x10000000, 0x8000000, 0x4000000,
0x2000000, 0x1000000, 0x800000, 0x400000,
0x200000, 0x100000, 0x80000, 0x40000,
0x20000, 0x10000, 0x8000, 0x4000,
0x2000, 0x1000, 0x800, 0x400,
0x200, 0x100, 0x80, 0x40,
0x20, 0x10, 0x8, 0x4,
0x2, 0x1,
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



double i_Ecs(double M, double e, double *E, double *ecosE, double *esinE, int N){
    // simple conversion dbl-fx + bit shift
    // s^s+, y before s
    // 19.2 ms
    long long x=KR*e, y=0, t=remainder(M,2*M_PI)*R, Y;
    int k, n, s;
    for (n=0; n<N; n++){
        k = k_n[n];
        Y = y;
        s = t+y >> 63;
        t -= s ^ s + i64_a[k];
        y += s ^ s + (x>>k);
        x -= s ^ s + (Y>>k);
    }
    *ecosE = x * dR;
    *esinE = y * dR;
    *E = M + *esinE;
    return M + y*dR;
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
    double x, e;
    int N;
    x = atof(argv[1]);
    e = atof(argv[2]);
    N = atoi(argv[3]);

    double hh, M=x-e*sin(x);
    printf("E:               %.20g\n", x);
    printf("i64_Ecs          %.20g\n", i_Ecs(M, e, &hh, &hh, &hh, N));

    return 0;
}

