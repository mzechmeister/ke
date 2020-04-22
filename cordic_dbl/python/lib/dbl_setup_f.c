#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
gcc dbl_setup_f.c -lm -O2 -o dbl_setup_f  &&  ./dbl_setup_f 61
*/

void make_tables(int kmax){
    // prepare tables and scale correction
    int k, n, Nmax, k_n[2*(kmax+1)];
    double sec2a, K=1., f64_a[2*(kmax+1)];

    for (k=0, n=0; k<=kmax; k++){
        f64_a[k] = atan(1./(1LL<<k));
        k_n[n++] = k;
        sec2a = 1 + 1./(1LL<<k)/(1LL<<k);
        if (sec2a>1) {
            K /= sec2a;
            k_n[n++] = k;
        }
    }
    Nmax = n;

    // print parameters and tables

    printf("#define Nmax %d\n", Nmax);
    printf("#define K %-21.13a\n", K);
    printf("\n");

    printf("atan(a_k)\n");
    for (k=0; k<=kmax; ++k) printf("%-21.13a,%s", f64_a[k], (k+1)%4?" ":"\n");

    printf("\n\nshift values k_n\n");
    for (n=0; n<Nmax; ++n) printf("%2d,%s", k_n[n], (n+1)%16?" ":"\n");
    printf("\n");
}


int main(int argc, char *argv[]) {
    int kmax, bp;
    kmax = atoi(argv[1]);   // largest shift value
    make_tables(kmax);

    return 0;
}
