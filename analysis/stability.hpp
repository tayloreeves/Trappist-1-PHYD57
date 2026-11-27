/* stability.hpp - stability checks for one TRAPPIST-1 copy
 *
 * C++ (also compiles with nvcc for the GPU).
 * Same layout as the integrator: x[nd*i+k] = coordinate k of body i (nd = 2 or 3).
 * Body 0 = star, 1-7 = planets b-h. Any consistent units, G is passed in
 * (Taylor's units: AU, yr, Msun, G = 4 pi^2).
 *
 * Checks, every few steps:
 *   stop the run:
 *     COLLISION         bodies touch: d < R_i + R_j (planet-planet or planet-star)
 *     TIDAL_DISRUPTION  inside the Roche limit of the star or another planet
 *     CLOSE_ENCOUNTER   within 2 Roche-lobe radii (assignment criterion, can be switched off)
 *     ESCAPE            further than 2x planet h's starting distance (assignment)
 *     INFALL            closer than planet b's starting periapsis minus 2 Hill radii
 *                       (assignment criterion, can be switched off)
 *   don't stop, just flag:
 *     hill_flags        bit set for each planet pair inside its mutual Hill radius
 */
#pragma once
#include <cmath>

constexpr double STAB_PI = 3.14159265358979323846;

#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

constexpr int NB = 8;       // star + 7 planets
constexpr int NPAIR = 21;   // planet-planet pairs

enum { STABLE = 0, COLLISION, TIDAL_DISRUPTION, CLOSE_ENCOUNTER, ESCAPE, INFALL };

struct StabLimits {
    int nd;                     // 2 or 3 dimensions
    double G;
    double R[NB];               // physical radii
    double rho[NB];             // mean densities, from mass and radius
    double r_escape;            // 2 x starting distance of planet h
    double r_infall;            // planet b's starting periapsis - 2 Hill radii of b
    double encounter_factor;    // 2 Roche-lobe radii (assignment)
    double hill_factor;         // flag when d < hill_factor x mutual Hill radius
    bool use_close_encounter;   // true = stop at 2 Roche-lobe radii, false = only flag Hill entries
    bool use_infall;            // true = assignment infall rule on, false = rely on Roche limit only
};

struct StabResult {
    int code;                   // STABLE or the reason the run stopped
    int i, j;                   // bodies involved (-1 if none)
    unsigned int hill_flags;    // bit pair_index(i,j) set if that pair is inside its Hill radius
};

// ---- radii ----

// Hill radius: where the planet's gravity beats the star's tidal pull
HD inline double hill_radius(double r, double m, double M) {
    return r * cbrt(m / (3.0 * M));
}

// mutual Hill radius of two planets at distances r1, r2 from the star
HD inline double mutual_hill_radius(double r1, double r2, double m1, double m2, double M) {
    return cbrt((m1 + m2) / (3.0 * M)) * 0.5 * (r1 + r2);
}

// Roche LOBE radius (Eggleton 1983): the region a body's gravity controls in a
// two-body orbit. ~0.71 x Hill radius for planets. This is what the assignment's
// "2 Roche lobe radii" refers to - not the same as the Roche LIMIT.
HD inline double roche_lobe_radius(double r, double m, double M) {
    double q = m / M, q13 = cbrt(q), q23 = q13 * q13;
    return r * 0.49 * q23 / (0.6 * q23 + log(1.0 + q13));
}

// Roche LIMIT (fluid body): distance from a primary inside which tides rip the
// secondary apart. d = 2.44 R_primary (rho_primary / rho_secondary)^(1/3)
HD inline double roche_limit(double R_primary, double rho_primary, double rho_secondary) {
    return 2.44 * R_primary * cbrt(rho_primary / rho_secondary);
}

// index 0..20 for planet pair (i, j), 1 <= i < j <= 7
HD inline int pair_index(int i, int j) {
    return (i - 1) * (14 - i) / 2 + (j - i - 1);
}

// ---- setup, once per copy from the starting state ----
HD inline StabLimits stab_setup(const double* x, const double* v, const double* m,
                                       const double* R, int nd, double G) {
    StabLimits L;
    L.nd = nd;
    L.G = G;
    L.encounter_factor = 2.0;
    L.hill_factor = 1.0;
    L.use_close_encounter = true;
    L.use_infall = true;

    for (int i = 0; i < NB; i++) {
        L.R[i] = R[i];
        L.rho[i] = m[i] / (4.0 / 3.0 * STAB_PI * R[i] * R[i] * R[i]);
    }

    // escape: 2 x planet h's starting distance from the star
    double r2 = 0;
    for (int k = 0; k < nd; k++) {
        double d = x[nd*7 + k] - x[k];
        r2 += d * d;
    }
    L.r_escape = 2.0 * sqrt(r2);

    // infall: planet b's periapsis a(1-e) at the start, minus 2 of its Hill radii.
    // b naturally swings in to its periapsis every orbit, so that can't count as
    // infall; moving in by more than 2 Hill radii past it means its orbit changed.
    double r = 0, vv = 0, rv = 0;
    for (int k = 0; k < nd; k++) {
        double dr = x[nd*1 + k] - x[k], dv = v[nd*1 + k] - v[k];
        r += dr * dr;  vv += dv * dv;  rv += dr * dv;
    }
    r = sqrt(r);
    double mu = G * (m[0] + m[1]);
    double a = 1.0 / (2.0 / r - vv / mu);                 // vis-viva
    double e2 = 0;                                        // eccentricity vector
    for (int k = 0; k < nd; k++) {
        double dr = x[nd*1 + k] - x[k], dv = v[nd*1 + k] - v[k];
        double ek = ((vv - mu / r) * dr - rv * dv) / mu;
        e2 += ek * ek;
    }
    L.r_infall = a * (1.0 - sqrt(e2)) - 2.0 * hill_radius(a, m[1], m[0]);
    return L;
}

// ---- the check, every few steps ----
HD inline StabResult stab_check(const double* x, const double* m, const StabLimits& L) {
    StabResult res = {STABLE, -1, -1, 0u};
    const int nd = L.nd;
    double rs[NB];   // distance of each planet from the star

    // planet vs star
    for (int i = 1; i < NB; i++) {
        double r2 = 0;
        for (int k = 0; k < nd; k++) {
            double d = x[nd*i + k] - x[k];
            r2 += d * d;
        }
        double r = sqrt(r2);
        rs[i] = r;
        int code = STABLE;
        if (r < L.R[0] + L.R[i])                                   code = COLLISION;
        else if (r < roche_limit(L.R[0], L.rho[0], L.rho[i]))     code = TIDAL_DISRUPTION;
        else if (r > L.r_escape)                                    code = ESCAPE;
        else if (L.use_infall && r < L.r_infall)                   code = INFALL;
        if (code != STABLE) { res.code = code; res.i = 0; res.j = i; return res; }
    }

    // planet vs planet
    for (int i = 1; i < NB; i++) {
        for (int j = i + 1; j < NB; j++) {
            double d2 = 0;
            for (int k = 0; k < nd; k++) {
                double d = x[nd*i + k] - x[nd*j + k];
                d2 += d * d;
            }
            double d = sqrt(d2);

            double rl_ij = roche_limit(L.R[i], L.rho[i], L.rho[j]);   // j torn apart by i
            double rl_ji = roche_limit(L.R[j], L.rho[j], L.rho[i]);   // i torn apart by j
            int code = STABLE;
            if (d < L.R[i] + L.R[j])                     code = COLLISION;
            else if (d < (rl_ij > rl_ji ? rl_ij : rl_ji))  code = TIDAL_DISRUPTION;
            else if (L.use_close_encounter) {
                double li = roche_lobe_radius(rs[i], m[i], m[0]);
                double lj = roche_lobe_radius(rs[j], m[j], m[0]);
                if (d < L.encounter_factor * (li > lj ? li : lj)) code = CLOSE_ENCOUNTER;
            }
            if (code != STABLE) { res.code = code; res.i = i; res.j = j; return res; }

            // red flag only: inside the mutual Hill radius
            if (d < L.hill_factor * mutual_hill_radius(rs[i], rs[j], m[i], m[j], m[0]))
                res.hill_flags |= 1u << pair_index(i, j);
        }
    }
    return res;
}

inline const char* stab_name(int code) {
    switch (code) {
        case STABLE:           return "stable";
        case COLLISION:        return "collision";
        case TIDAL_DISRUPTION: return "tidal_disruption";
        case CLOSE_ENCOUNTER:  return "close_encounter";
        case ESCAPE:           return "escape";
        case INFALL:           return "infall";
        default:               return "unknown";
    }
}
