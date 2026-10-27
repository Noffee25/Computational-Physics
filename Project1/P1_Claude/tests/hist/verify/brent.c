/* verifier: Brent cycle detection for y <- 4y(1-y) in IEEE double,
   y0 = frac(seed*0.6180339887498949). Prints seed, tail mu, cycle lam, entry state */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static inline double f(double y) { return 4.0 * y * (1.0 - y); }
int main(int argc, char **argv) {
    for (int a = 1; a < argc; a++) {
        long long seed = atoll(argv[a]);
        double y0 = fmod((double)seed * 0.6180339887498949, 1.0);
        if (y0 <= 0.0 || y0 >= 1.0) y0 = 0.3141592653589793;
        long long power = 1, lam = 1;
        double t = y0, h = f(y0);
        while (t != h) {
            if (power == lam) { t = h; power *= 2; lam = 0; }
            h = f(h); lam++;
            if (power > 8000000000LL) break;
        }
        t = y0; h = y0;
        for (long long i = 0; i < lam; i++) h = f(h);
        long long mu = 0;
        while (t != h) { t = f(t); h = f(h); mu++; }
        printf("%lld %lld %lld %.17g\n", seed, mu, lam, t);
        fflush(stdout);
    }
    return 0;
}
