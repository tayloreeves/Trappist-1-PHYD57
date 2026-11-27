// stability.hpp - checks if a TRAPPIST-1 copy has gone unstable
// index 0 = star, 1-7 = planets b-h. Works in code units.

#pragma once
#include <cmath>

#ifdef __CUDACC__
#define HD __host__ __device__   // so it also compiles for the GPU
#else
#define HD
#endif

constexpr int NBODY = 8;

enum StabilityCode { STABLE = 0, CLOSE_ENCOUNTER = 1, ESCAPE = 2, INFALL = 3 };

struct StabilityLimits {
    double r_max;          // 2x starting distance of planet h
    double r_min;          // 0.95x starting distance of planet b
    double roche_factor;   // 2 Roche radii
};

struct StabilityResult {
    int code;
    int i, j;   // planets involved (-1 if none)
};

// Hill radius used as the Roche radius (check with L12)
HD inline double roche_radius(double r, double m, double M) {
    return r * cbrt(m / (3.0 * M));
}

// set limits from the starting positions of a copy
// planet b starts right on the infall line, so we use 95% of it
HD inline StabilityLimits make_limits(const double* x, const double* y, const double* z,
                                      double roche_factor = 2.0, double infall_margin = 0.95) {
    double rmin = 1e300, rmax = 0.0;
    for (int i = 1; i < NBODY; ++i) {
        double dx = x[i] - x[0], dy = y[i] - y[0], dz = z[i] - z[0];
        double r = sqrt(dx*dx + dy*dy + dz*dz);
        if (r < rmin) rmin = r;
        if (r > rmax) rmax = r;
    }
    StabilityLimits L;
    L.r_max = 2.0 * rmax;
    L.r_min = infall_margin * rmin;
    L.roche_factor = roche_factor;
    return L;
}

// returns the first problem found, or STABLE
HD inline StabilityResult check_stability(const double* x, const double* y, const double* z,
                                          const double* m, const StabilityLimits& L) {
    double rx[NBODY], ry[NBODY], rz[NBODY], rH[NBODY];

    // escape / infall (positions relative to the star)
    for (int i = 1; i < NBODY; ++i) {
        rx[i] = x[i] - x[0];
        ry[i] = y[i] - y[0];
        rz[i] = z[i] - z[0];
        double r = sqrt(rx[i]*rx[i] + ry[i]*ry[i] + rz[i]*rz[i]);
        if (r > L.r_max) return {ESCAPE, i, -1};
        if (r < L.r_min) return {INFALL, i, -1};
        rH[i] = roche_radius(r, m[i], m[0]);
    }

    // close encounters between planet pairs (uses the bigger Roche radius)
    for (int i = 1; i < NBODY; ++i) {
        for (int j = i + 1; j < NBODY; ++j) {
            double dx = rx[i] - rx[j], dy = ry[i] - ry[j], dz = rz[i] - rz[j];
            double d2 = dx*dx + dy*dy + dz*dz;
            double lim = L.roche_factor * (rH[i] > rH[j] ? rH[i] : rH[j]);
            if (d2 < lim * lim) return {CLOSE_ENCOUNTER, i, j};
        }
    }
    return {STABLE, -1, -1};
}

inline const char* stability_name(int code) {
    switch (code) {
        case STABLE:          return "stable";
        case CLOSE_ENCOUNTER: return "close_encounter";
        case ESCAPE:          return "escape";
        case INFALL:          return "infall";
        default:              return "unknown";
    }
}
