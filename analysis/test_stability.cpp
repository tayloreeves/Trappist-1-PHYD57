// tests for stability.hpp using a fake TRAPPIST-1-like system
// g++ -O2 -std=c++17 test_stability.cpp -o test_stability && ./test_stability

#include "stability.hpp"
#include <cstdio>
#include <cmath>

int failures = 0;

void expect(const char* name, StabilityResult got, int code, int i, int j) {
    bool ok = got.code == code && got.i == i && got.j == j;
    printf("[%s] %-40s %s (i=%d, j=%d)\n", ok ? "PASS" : "FAIL",
           name, stability_name(got.code), got.i, got.j);
    if (!ok) failures++;
}

// planets on circular orbits, star at the origin (rough values, testing only)
void make_system(double* x, double* y, double* z, double* m) {
    double a[8]  = {0, 1.000, 1.369, 1.930, 2.535, 3.335, 4.058, 5.363};
    double mE[8] = {0, 1.37, 1.31, 0.39, 0.69, 1.04, 1.32, 0.33};
    double earth_over_star = 3.003e-6 / 0.0898;
    x[0] = y[0] = z[0] = 0; m[0] = 1;
    for (int i = 1; i < 8; i++) {
        double phi = 0.9 * i;
        x[i] = a[i] * cos(phi);
        y[i] = a[i] * sin(phi);
        z[i] = 0;
        m[i] = mE[i] * earth_over_star;
    }
}

int main() {
    double x[8], y[8], z[8], m[8];

    make_system(x, y, z, m);
    StabilityLimits L = make_limits(x, y, z);
    expect("normal system", check_stability(x, y, z, m, L), STABLE, -1, -1);

    make_system(x, y, z, m);
    x[2] = x[1] + 0.001; y[2] = y[1];
    expect("planet c right next to b", check_stability(x, y, z, m, L), CLOSE_ENCOUNTER, 1, 2);

    make_system(x, y, z, m);
    x[7] *= 2.1; y[7] *= 2.1;
    expect("planet h too far out", check_stability(x, y, z, m, L), ESCAPE, 7, -1);

    make_system(x, y, z, m);
    x[1] *= 0.9; y[1] *= 0.9;
    expect("planet b too far in", check_stability(x, y, z, m, L), INFALL, 1, -1);

    // shifting everything shouldn't matter (barycentric frame)
    make_system(x, y, z, m);
    for (int i = 0; i < 8; i++) { x[i] += 0.3; y[i] -= 0.2; z[i] += 0.1; }
    expect("whole system shifted", check_stability(x, y, z, m, L), STABLE, -1, -1);

    // 2.5 Roche radii apart is still ok
    make_system(x, y, z, m);
    double rH1 = roche_radius(1.0, m[1], m[0]);
    x[2] = x[1] + 2.5 * rH1 * cos(0.9); y[2] = y[1] + 2.5 * rH1 * sin(0.9);
    expect("planets 2.5 Roche radii apart", check_stability(x, y, z, m, L), STABLE, -1, -1);

    // without the 95% margin, round-off alone flags planet b
    make_system(x, y, z, m);
    StabilityLimits Lstrict = make_limits(x, y, z, 2.0, 1.0);
    for (int i = 0; i < 8; i++) { x[i] += 0.3; y[i] -= 0.2; z[i] += 0.1; }
    expect("no margin -> false infall", check_stability(x, y, z, m, Lstrict), INFALL, 1, -1);

    printf("\nRoche radius of b: %.4e (%.1f%% of a_b)\n", rH1, 100 * rH1);
    printf(failures == 0 ? "All tests passed.\n" : "SOME TESTS FAILED.\n");
    return failures;
}
