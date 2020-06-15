#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
gcc dbl_setup_i.c -lm -O2 -o dbl_setup_i  &&  ./dbl_setup_i 61 61
*/

void make_tables(int kmax, int B){
    // prepare tables and scale correction
    // B: binary point
    int k, n, N, k_n[2*(kmax+1)];
    long long R=1LL<<B, i64_a[2*(kmax+1)];
    double sec2a, K=1., dR=1./R;

    for (k=0, n=0; k<=kmax; k++){
        i64_a[k] = round(atan(1./(1LL<<k)) * R);
        k_n[n++] = k;
        sec2a = 1 + 1./(1LL<<k)/(1LL<<k);
        if (sec2a>1) {
            K /= sec2a;
            k_n[n++] = k;
        }
    }
    N = n;

    // print parameters and tables

    printf("#define R (1LL<<%d)\n", B);
    printf("#define dR (1./R)\n");
    printf("#define N %d\n", N);
    printf("#define KR %-21.13a\n", K*R);
    printf("// R: %lld = 0x%llx = 2^%d\n", R, R, B);
    printf("// K: %lf\n", K);
    printf("\n");

    printf("atan(a_k)\n");
    for (k=0; k<=kmax; ++k) printf("0x%llx,%s", i64_a[k], (k+1)%4?" ":"\n");

    printf("\n\nshift values k_n\n");
    for (n=0; n<N; ++n) printf("%2d,%s", k_n[n], (n+1)%16?" ":"\n");
    printf("\n");

}


int main(int argc, char *argv[]) {
    int kmax, bp;
    kmax = atoi(argv[1]);   // largest shift value
    bp = atoi(argv[2]);     // location of binary point
    make_tables(kmax, bp);

    return 0;
}
