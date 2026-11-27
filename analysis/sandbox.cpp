// sandbox.cpp - small Forest-Ruth integrator for trying edge cases
// uses stability.hpp, runs copies in parallel with OpenMP
//
// g++ -O3 -std=c++17 -fopenmp sandbox.cpp -o sandbox
// ./sandbox convergence
// ./sandbox massscan <years> <copies>
// ./sandbox spacing <years> <copies>
//
// code units: G = 1, star mass = 1, a_b ~ 1  (P_b = 2pi, 1 yr = 1519 time units)
// orbits are flat (z = 0) since TRAPPIST-1 is nearly coplanar

#include "stability.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <string>
#include <omp.h>

const int N = NB;
const double TWO_PI = 2.0 * M_PI;

// Agol et al. 2021
const double M_STAR_SUN = 0.0898;
const double M_EARTH_SUN = 3.0034896e-6;
const double P_DAYS[7]  = {1.510826, 2.421937, 4.049219, 6.101013, 9.207540, 12.352446, 18.772866};
const double M_EARTH[7] = {1.374, 1.308, 0.388, 0.692, 1.039, 1.321, 0.326};
const double ECOSW[7] = {-0.00215, 0.00055, -0.00496, 0.00433, -0.00840, 0.00380, -0.00365};
const double ESINW[7] = {0.00217, 0.00001, 0.00267, -0.00461, -0.00051, 0.00128, -0.00002};
const double R_EARTH[7] = {1.116, 1.097, 0.788, 0.920, 1.045, 1.129, 0.755};

// radii in code units (lengths in units of a_b = 0.01154 AU)
const double A_B_AU = 0.01154;
const double R_STAR_CODE = 0.1192 * 0.00465047 / A_B_AU;
const double R_EARTH_CODE = 4.26352e-5 / A_B_AU;

const double YEAR = 365.25 / (P_DAYS[0] / TWO_PI);

struct System { double x[N], y[N], z[N], vx[N], vy[N], vz[N], m[N]; };

struct Setup {
    double mass_factor = 1.0;   // multiply all planet masses
    int stretch_planet = -1;    // planet whose period gets changed
    double period_scale = 1.0;
};

// random phase and orientation, velocities from the ellipse equation
System make_system(const Setup& S, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> U(0.0, TWO_PI);
    System s{};
    s.m[0] = 1.0;
    for (int i = 1; i < N; i++) {
        double P = P_DAYS[i-1] / P_DAYS[0];
        if (i == S.stretch_planet) P *= S.period_scale;
        s.m[i] = S.mass_factor * M_EARTH[i-1] * M_EARTH_SUN / M_STAR_SUN;
        double GM = 1.0 + s.m[i];
        double a = cbrt(GM * P * P);   // Kepler's 3rd law
        double e = hypot(ECOSW[i-1], ESINW[i-1]);
        double w = U(rng), f = U(rng);

        double p = a * (1 - e*e);
        double r = p / (1 + e * cos(f));
        double vr = sqrt(GM / p) * e * sin(f);
        double vt = sqrt(GM * p) / r;   // L / r
        double th = w + f;

        s.x[i] = r * cos(th);
        s.y[i] = r * sin(th);
        s.vx[i] = vr * cos(th) - vt * sin(th);
        s.vy[i] = vr * sin(th) + vt * cos(th);
    }

    // move to the centre-of-mass frame
    double M = 0, cx = 0, cy = 0, px = 0, py = 0;
    for (int i = 0; i < N; i++) {
        M += s.m[i];
        cx += s.m[i] * s.x[i];  cy += s.m[i] * s.y[i];
        px += s.m[i] * s.vx[i]; py += s.m[i] * s.vy[i];
    }
    for (int i = 0; i < N; i++) {
        s.x[i] -= cx / M;  s.y[i] -= cy / M;
        s.vx[i] -= px / M; s.vy[i] -= py / M;
    }
    return s;
}

// gravity, each pair once (Newton's 3rd law)
void accel(const System& s, double* ax, double* ay, double* az) {
    for (int i = 0; i < N; i++) ax[i] = ay[i] = az[i] = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            double dx = s.x[j] - s.x[i], dy = s.y[j] - s.y[i], dz = s.z[j] - s.z[i];
            double r2 = dx*dx + dy*dy + dz*dz;
            double inv3 = 1.0 / (r2 * sqrt(r2));
            ax[i] += s.m[j] * dx * inv3;  ay[i] += s.m[j] * dy * inv3;  az[i] += s.m[j] * dz * inv3;
            ax[j] -= s.m[i] * dx * inv3;  ay[j] -= s.m[i] * dy * inv3;  az[j] -= s.m[i] * dz * inv3;
        }
}

double energy(const System& s) {
    double E = 0;
    for (int i = 0; i < N; i++) {
        E += 0.5 * s.m[i] * (s.vx[i]*s.vx[i] + s.vy[i]*s.vy[i] + s.vz[i]*s.vz[i]);
        for (int j = i + 1; j < N; j++) {
            double dx = s.x[j] - s.x[i], dy = s.y[j] - s.y[i], dz = s.z[j] - s.z[i];
            E -= s.m[i] * s.m[j] / sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    return E;
}

// Forest-Ruth coefficients (lecture 9), computed exactly instead of rounded
struct FR { double c[4], d[3]; };
FR forest_ruth() {
    double cr = cbrt(2.0);
    double w1 = 1.0 / (2.0 - cr), w0 = -cr / (2.0 - cr);
    return {{w1/2, (w0+w1)/2, (w0+w1)/2, w1/2}, {w1, w0, w1}};
}

// drift c1, kick d1, drift c2, kick d2, drift c3, kick d3, drift c4
void step(System& s, double dt, const FR& k) {
    double ax[N], ay[N], az[N];
    for (int sub = 0; sub < 4; sub++) {
        for (int i = 0; i < N; i++) {
            s.x[i] += k.c[sub] * dt * s.vx[i];
            s.y[i] += k.c[sub] * dt * s.vy[i];
            s.z[i] += k.c[sub] * dt * s.vz[i];
        }
        if (sub == 3) break;
        accel(s, ax, ay, az);
        for (int i = 0; i < N; i++) {
            s.vx[i] += k.d[sub] * dt * ax[i];
            s.vy[i] += k.d[sub] * dt * ay[i];
            s.vz[i] += k.d[sub] * dt * az[i];
        }
    }
}

struct RunResult { int code, i, j; double t_years, max_dE; };

// run one copy until t_end or until it goes unstable
// copy into the integrator layout pos[3*i+k] that stability.hpp uses
void pack(const System& s, double* pos, double* vel) {
    for (int i = 0; i < N; i++) {
        pos[3*i] = s.x[i];   pos[3*i+1] = s.y[i];   pos[3*i+2] = s.z[i];
        vel[3*i] = s.vx[i];  vel[3*i+1] = s.vy[i];  vel[3*i+2] = s.vz[i];
    }
}

RunResult run_copy(System s, double dt, double t_end_years, int check_every) {
    FR k = forest_ruth();
    double pos[3*N], vel[3*N], R[N];
    R[0] = R_STAR_CODE;
    for (int i = 1; i < N; i++) R[i] = R_EARTH[i-1] * R_EARTH_CODE;
    pack(s, pos, vel);
    StabLimits L = stab_setup(pos, vel, s.m, R, 3, 1.0);   // G = 1 in code units
    double E0 = energy(s), max_dE = 0;
    long nsteps = (long)(t_end_years * YEAR / dt);
    for (long n = 1; n <= nsteps; n++) {
        step(s, dt, k);
        if (n % check_every == 0) {
            pack(s, pos, vel);
            StabResult r = stab_check(pos, s.m, L);
            if (r.code != STABLE) return {r.code, r.i, r.j, n * dt / YEAR, max_dE};
        }
        if (n % 10000 == 0) {
            double dE = fabs((energy(s) - E0) / E0);
            if (dE > max_dE) max_dE = dE;
        }
    }
    return {STABLE, -1, -1, t_end_years, max_dE};
}

// halving dt should make the energy error ~16x smaller (4th order)
void convergence() {
    printf("dt          dt/P_b   max|dE/E|    ratio\n");
    double prev = 0;
    for (int fr : {25, 50, 100, 200}) {
        std::mt19937_64 rng(12345);
        System s = make_system(Setup{}, rng);
        double dt = TWO_PI / fr;
        FR k = forest_ruth();
        double E0 = energy(s), maxdE = 0;
        long nsteps = (long)(10 * YEAR / dt);
        for (long n = 1; n <= nsteps; n++) {
            step(s, dt, k);
            if (n % 50 == 0) maxdE = std::max(maxdE, fabs((energy(s) - E0) / E0));
        }
        printf("%-11.5f 1/%-6d %-12.3e %s\n", dt, fr, maxdE,
               prev > 0 ? std::to_string(prev / maxdE).c_str() : "-");
        prev = maxdE;
    }
}

// run every setup x copy in parallel, print a summary and save a csv
void scan(const char* label, const std::vector<Setup>& setups, const std::vector<double>& values,
          double years, int ncopies, const char* csv) {
    int nset = setups.size(), total = nset * ncopies;
    std::vector<RunResult> res(total);
    double dt = TWO_PI / 100.0;   // P_b / 100
    double t0 = omp_get_wtime();

    #pragma omp parallel for schedule(dynamic)
    for (int idx = 0; idx < total; idx++) {
        int si = idx / ncopies, c = idx % ncopies;
        std::mt19937_64 rng(1000 + c);   // same start for copy c in every setup
        System s = make_system(setups[si], rng);
        res[idx] = run_copy(s, dt, years, 4);
    }
    double wall = omp_get_wtime() - t0;
    double simyears = 0;
    for (auto& r : res) simyears += r.t_years;

    FILE* f = fopen(csv, "w");
    fprintf(f, "%s,copy_id,unstable,t_instability,reason,planet_i,planet_j,max_dE\n", label);
    printf("\n%-14s unstable   median t [yr]   max|dE/E|\n", label);
    for (int si = 0; si < nset; si++) {
        std::vector<double> ts;
        double mde = 0;
        int nu = 0;
        for (int c = 0; c < ncopies; c++) {
            RunResult& r = res[si * ncopies + c];
            fprintf(f, "%g,%d,%d,%.6e,%s,%d,%d,%.3e\n", values[si], c, r.code != STABLE,
                    r.t_years, stab_name(r.code), r.i, r.j, r.max_dE);
            if (r.code != STABLE) { nu++; ts.push_back(r.t_years); }
            else mde = std::max(mde, r.max_dE);
        }
        std::sort(ts.begin(), ts.end());
        char med[32] = "-";
        if (!ts.empty()) snprintf(med, 32, "%.3g", ts[ts.size() / 2]);
        printf("%-14g %2d/%-7d %-15s %.2e\n", values[si], nu, ncopies, med, mde);
    }
    fclose(f);
    printf("%.1f s on %d threads, %.0f simulated yr/s\n", wall, omp_get_max_threads(), simyears / wall);
}

int main(int argc, char** argv) {
    if (argc < 2) { printf("modes: convergence | massscan <yr> <n> | spacing <yr> <n>\n"); return 1; }
    std::string mode = argv[1];
    double years = argc > 2 ? atof(argv[2]) : 100;
    int ncopies = argc > 3 ? atoi(argv[3]) : 8;

    if (mode == "convergence") {
        convergence();
    } else if (mode == "massscan") {
        std::vector<double> f = {1, 3, 10, 30, 100};
        std::vector<Setup> S;
        for (double v : f) { Setup s; s.mass_factor = v; S.push_back(s); }
        scan("mass_factor", S, f, years, ncopies, "massscan.csv");
    } else if (mode == "spacing") {
        // move planet h off its 3:2 resonance with g
        std::vector<double> f = {0.90, 0.95, 0.98, 1.00, 1.02, 1.05};
        std::vector<Setup> S;
        for (double v : f) { Setup s; s.stretch_planet = 7; s.period_scale = v; S.push_back(s); }
        scan("h_period_scale", S, f, years, ncopies, "spacing.csv");
    }
    return 0;
}
