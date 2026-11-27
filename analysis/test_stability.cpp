// tests for stability.hpp with real TRAPPIST-1 values, in AU, yr, Msun (G = 4 pi^2)
// g++ -O2 -std=c++17 test_stability.cpp -o test_stability && ./test_stability

#include "stability.hpp"
#include <cstdio>
#include <cstring>

constexpr int ND = 3;
const double G = 4.0 * STAB_PI * STAB_PI;

// Agol et al. 2021
const double M_STAR = 0.0898;                       // Msun
const double R_STAR = 0.1192 * 0.00465047;          // Rsun -> AU
const double A[7]  = {0.01154, 0.01580, 0.02227, 0.02925, 0.03849, 0.04683, 0.06189};
const double ME[7] = {1.374, 1.308, 0.388, 0.692, 1.039, 1.321, 0.326};
const double RE[7] = {1.116, 1.097, 0.788, 0.920, 1.045, 1.129, 0.755};   // Earth radii
const double MEARTH = 3.0034896e-6, REARTH = 4.26352e-5;                  // Msun, AU

double x[NB*ND], v[NB*ND], m[NB], R[NB];
int failures = 0;

// circular orbits in the xy plane, star at rest at the origin
void make_system() {
    memset(x, 0, sizeof x);
    memset(v, 0, sizeof v);
    m[0] = M_STAR;  R[0] = R_STAR;
    for (int i = 1; i < NB; i++) {
        double phi = 0.9 * i, a = A[i-1];
        m[i] = ME[i-1] * MEARTH;
        R[i] = RE[i-1] * REARTH;
        double vc = sqrt(G * (M_STAR + m[i]) / a);
        x[ND*i] = a * cos(phi);   x[ND*i + 1] = a * sin(phi);
        v[ND*i] = -vc * sin(phi); v[ND*i + 1] = vc * cos(phi);
    }
}

// put planet j at distance d from planet i, pointing away from the star
void place_near(int i, int j, double d) {
    double r = sqrt(x[ND*i]*x[ND*i] + x[ND*i+1]*x[ND*i+1]);
    for (int k = 0; k < 2; k++) x[ND*j + k] = x[ND*i + k] * (1 + d / r);
}

void expect(const char* name, StabResult r, int code, int i, int j, bool want_flag) {
    bool flag = r.hill_flags != 0;
    bool ok = r.code == code && r.i == i && r.j == j && flag == want_flag;
    printf("[%s] %-38s %-17s i=%2d j=%2d hill_flag=%d\n", ok ? "PASS" : "FAIL",
           name, stab_name(r.code), r.i, r.j, (int)flag);
    if (!ok) failures++;
}

int main() {
    make_system();
    StabLimits L = stab_setup(x, v, m, R, ND, G);
    double rH_b = hill_radius(A[0], m[1], m[0]);
    double rL_b = roche_lobe_radius(A[0], m[1], m[0]);
    double mH_bc = mutual_hill_radius(A[0], A[0], m[1], m[2], m[0]);
    double rlim_star_b = roche_limit(R[0], L.rho[0], L.rho[1]);
    double rlim_bc = roche_limit(R[1], L.rho[1], L.rho[2]);

    printf("planet b: Hill radius %.3e AU, Roche lobe %.3e AU (%.2f x Hill)\n", rH_b, rL_b, rL_b / rH_b);
    printf("b-c mutual Hill radius (b at a_b) %.3e AU; 2 Roche lobes %.3e AU\n", mH_bc, 2 * rL_b);
    printf("Roche limit of star for b: %.3e AU = %.2f a_b = %.1f R_star\n",
           rlim_star_b, rlim_star_b / A[0], rlim_star_b / R[0]);
    printf("Roche limit b-c: %.3e AU = %.2f (R_b+R_c)\n", rlim_bc, rlim_bc / (R[1] + R[2]));
    printf("infall radius %.5f AU = %.3f a_b; escape radius %.4f AU\n\n",
           L.r_infall, L.r_infall / A[0], L.r_escape);

    make_system();
    expect("normal system", stab_check(x, m, L), STABLE, -1, -1, false);

    make_system();
    place_near(1, 2, 1.5 * rL_b);
    expect("c within 2 Roche lobes of b", stab_check(x, m, L), CLOSE_ENCOUNTER, 1, 2, false);

    // with the assignment's close-encounter stop off, Hill entry is only a flag
    StabLimits Lflag = L;
    Lflag.use_close_encounter = false;
    make_system();
    place_near(1, 2, 0.8 * mH_bc);
    expect("c inside b's Hill radius (flag only)", stab_check(x, m, Lflag), STABLE, -1, -1, true);

    make_system();
    place_near(1, 2, 0.9 * (R[1] + R[2]));
    expect("b and c touching", stab_check(x, m, Lflag), COLLISION, 1, 2, false);

    make_system();
    place_near(1, 2, 0.5 * (rlim_bc + R[1] + R[2]));
    expect("c inside b-c Roche limit", stab_check(x, m, Lflag), TIDAL_DISRUPTION, 1, 2, false);

    StabLimits Lnoinf = L;
    Lnoinf.use_infall = false;
    make_system();
    x[ND*1] *= 0.9 * rlim_star_b / A[0];  x[ND*1 + 1] *= 0.9 * rlim_star_b / A[0];
    expect("b inside star's Roche limit", stab_check(x, m, Lnoinf), TIDAL_DISRUPTION, 0, 1, false);

    make_system();
    x[ND*1] *= 0.9;  x[ND*1 + 1] *= 0.9;
    expect("b moved 10% inward (infall)", stab_check(x, m, L), INFALL, 0, 1, false);

    make_system();
    x[ND*7] *= 2.1;  x[ND*7 + 1] *= 2.1;
    expect("h beyond 2x its start", stab_check(x, m, L), ESCAPE, 0, 7, false);

    make_system();
    for (int i = 0; i < NB; i++) { x[ND*i] += 0.01; x[ND*i + 1] -= 0.02; x[ND*i + 2] += 0.005; }
    expect("whole system shifted in 3D", stab_check(x, m, L), STABLE, -1, -1, false);

    make_system();
    x[ND*7 + 2] = 0.5 * A[6];   // h pushed out of the plane, still within limits
    expect("h pushed out of plane (3D)", stab_check(x, m, L), STABLE, -1, -1, false);

    printf("\npair_index check: (1,2)=%d (1,7)=%d (2,3)=%d (6,7)=%d\n",
           pair_index(1, 2), pair_index(1, 7), pair_index(2, 3), pair_index(6, 7));
    printf(failures == 0 ? "All tests passed.\n" : "SOME TESTS FAILED.\n");
    return failures;
}
