/* hsim.c : independent re-implementation (verifier) of the deep-histogram
 * test.  Generators are coded from their definitions, NOT from rnglib:
 *   randu : s <- 65539 s mod 2^31,  x = s/2^31
 *   randm : s <- 6539  s mod 2^31,  x = s/2^31
 *   minstd: s <- 16807 s mod (2^31-1), x = s/(2^31-1)
 *   rand  : MINSTD state s, x = floor((s-1)/2^9)/2^22   (claimed identity)
 *   logis : y <- 4y(1-y) (double), x = (2/pi) asin(sqrt(y)),
 *           y0 = frac(seed*0.6180339887498949)
 * Optional single precision storage: x <- (double)(float)x, 1.0 -> largest float < 1.
 * Bins: idx = (int)(x*B), clamped to [0,B-1].
 * usage: hsim GEN SEED B N NSEG SINGLE
 *   runs NSEG consecutive segments of N draws (histogram reset each time) and
 *   prints one line per segment: seg D z nempty min max
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage\n"); return 1; }
    const char *gen = argv[1];
    long long seed = atoll(argv[2]);
    long long B = atoll(argv[3]);
    long long N = atoll(argv[4]);
    long long nseg = atoll(argv[5]);
    int single = atoi(argv[6]);
    int g = 0;
    if (!strcmp(gen, "randu")) g = 1; else if (!strcmp(gen, "randm")) g = 2;
    else if (!strcmp(gen, "minstd")) g = 4; else if (!strcmp(gen, "rand")) g = 3;
    else if (!strcmp(gen, "logis")) g = 9; else if (!strcmp(gen, "sm64")) g = 10; else { fprintf(stderr, "gen?\n"); return 1; }
    const uint64_t M31 = 2147483648ULL, PM = 2147483647ULL;
    uint64_t s;
    double y = 0.0;
    if (g == 1 || g == 2) s = (uint64_t)(((seed % (long long)M31) + (long long)M31) % (long long)M31);
    else { s = (uint64_t)(((seed % (long long)PM) + (long long)PM) % (long long)PM); if (s == 0) s = 1; }
    if (g == 9) { y = fmod((double)seed * 0.6180339887498949, 1.0); if (y <= 0.0 || y >= 1.0) y = 0.3141592653589793; }
    uint64_t sm = (uint64_t)seed; sm = (sm ^ (sm >> 33)) * 0xFF51AFD7ED558CCDULL; sm = (sm ^ (sm >> 33)) * 0xC4CEB9FE1A85EC53ULL; sm ^= sm >> 33; /* hashed seed: distinct streams */
    uint32_t *h = (uint32_t *)malloc(sizeof(uint32_t) * B);
    double dB = (double)B;
    const double ONEM = 1.0 - 2.220446049250313e-16;
    const double ONEMSP = 0.99999994039535522;
    for (long long sg = 0; sg < nseg; sg++) {
        memset(h, 0, sizeof(uint32_t) * B);
        for (long long i = 0; i < N; i++) {
            double x;
            switch (g) {
            case 1: s = (65539ULL * s) & (M31 - 1); x = (double)s / 2147483648.0; break;
            case 2: s = (6539ULL * s) & (M31 - 1); x = (double)s / 2147483648.0; break;
            case 4: s = (16807ULL * s) % PM; x = (double)s / 2147483647.0; break;
            case 3: s = (16807ULL * s) % PM; x = (double)((s - 1) >> 9) / 4194304.0; break;
            case 10: { /* splitmix64 reference generator, 53-bit output */
                uint64_t zz = (sm += 0x9E3779B97F4A7C15ULL);
                zz = (zz ^ (zz >> 30)) * 0xBF58476D1CE4E5B9ULL;
                zz = (zz ^ (zz >> 27)) * 0x94D049BB133111EBULL;
                zz ^= zz >> 31;
                x = (double)(zz >> 11) / 9007199254740992.0; break; }
            default:
                y = 4.0 * y * (1.0 - y);
                x = (2.0 / 3.141592653589793) * asin(sqrt(y));
            }
            if (x >= 1.0) x = ONEM;
            if (x < 0.0) x = 0.0;
            if (single) { float f = (float)x; x = (double)f; if (x >= 1.0) x = ONEMSP; }
            long long idx = (long long)(x * dB);
            if (idx >= B) idx = B - 1;
            if (idx < 0) idx = 0;
            h[idx]++;
        }
        double lam = (double)N / dB, s2 = 0.0;
        long long ne = 0; uint32_t mn = h[0], mx = h[0];
        for (long long b = 0; b < B; b++) {
            double d = (double)h[b] - lam; s2 += d * d;
            if (h[b] == 0) ne++;
            if (h[b] < mn) mn = h[b];
            if (h[b] > mx) mx = h[b];
        }
        double D = s2 / (double)(B - 1) / lam;
        double z = (s2 / lam - (double)(B - 1)) / sqrt(2.0 * (double)(B - 1) * (1.0 - 1.0 / (double)N));
        printf("%lld %.10f %.6f %lld %u %u\n", sg, D, z, ne, mn, mx);
        fflush(stdout);
    }
    return 0;
}
