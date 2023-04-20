/* NOTE: Literal translation of Kepler solver using Universal formulation using Stumpff functions
 * from ABIE integrator_wisdom_holman.c. Please refer to that file for more information */
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <float.h>

#define NMAX 20
#define TOL 1e-10
#define TOL_ENERGY 0.0
#define MAX_KEPLER_ITERATION  100
#define real double


 real vector_norm(const real *vec, size_t N) {
     real sum = 0.0;
     for (size_t i = 0; i < N; i++) sum += (vec[i] * vec[i]);
     return sqrt(sum);
 }
 real dot(const real *vec1, const real *vec2) {
     return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
 }

 real cross_norm(const real *vec1, const real *vec2) {
     real c0 = vec1[1] * vec2[2] - vec1[2] * vec2[1];
     real c1 = vec1[2] * vec2[0] - vec1[0] * vec2[2];
     real c2 = vec1[0] * vec2[1] - vec1[1] * vec2[0];
     return sqrt(c0 * c0 + c1 * c1 + c2 * c2);
 }

 real sign(real x) {
     return (x > 0) - (x < 0);
 }

void compute_c2c3(real psi, real *cs) {
    real c2, c3;
    if (psi > 1.e-10){
        c2 = (1.0 - cos(sqrt(psi))) / psi;
        c3 = (sqrt(psi) - sin(sqrt(psi))) / sqrt(pow(psi, 3.0));

    } else {
        if (psi < -1.e-6){
            c2 = (1 - cosh(sqrt(-psi))) / psi;
            c3 = (sinh(sqrt(-psi)) - sqrt(-psi)) / sqrt(-pow(psi, 3.0));
        } else {
           c2 = 0.5;
           c3 = 1.0 / 6.0;
        }
    }
    cs[0] = c2; cs[1] = c3;
    return;
}



void propagate_kepler_wh(real x_star[NMAX][3], real v_star[NMAX][3], real gm, real dt, size_t N, size_t particle_id) {

    real vr0[3];
    real vv0[3];
    for (size_t d = 0; d < 3; d++) vr0[d] = x_star[particle_id][d];
    for (size_t d = 0; d < 3; d++) vv0[d] = v_star[particle_id][d];

    // Compute the magnitude of the initial position and velocity vectors:
    real r0 = vector_norm(vr0, 3);
    real v0 = vector_norm(vv0, 3);

    // Precompute sqrt(gm):
    real sqrtgm = sqrt(gm);

    // Initial value of the Keplerian energy:
    real xi = 0.5 * v0 * v0  - gm / r0;

    // Semimajor axis:
    real sma = -gm / (2 * xi);
    real alpha = 1.0 / sma;
    real chi0;
    if (alpha > TOL_ENERGY){
        // Elliptic orbits:
        chi0 = sqrtgm * dt * alpha;
    } else if (alpha < TOL_ENERGY) {
        // Hyperbolic orbits:
        chi0 = sign(dt) * (sqrt(-sma) * log(-2.0 * gm * alpha * dt / ( dot(vr0, vv0) + sqrt(-gm * sma) * (1.0 - r0 * alpha))));
    } else {
        // Parabolic orbits:
        real cn = cross_norm(vr0, vv0);
        real p = cn * cn / gm;
        real s = 0.5 * atan( 1.0 / (3.0 * sqrt(gm / pow(p, 3.0)) * dt));
        real w = atan(pow(tan(s), 1.0 / 3.0));
        chi0 = sqrt(p) * 2 / tan(2 * w);
    }

    // Solve Kepler's equation:
    real cs[2];
    real chi = 0.0, r = 0.0, c2 = 0.0, c3 = 0.0, psi = 0.0;
    for (size_t j = 0; j < MAX_KEPLER_ITERATION; j++) {
        // Compute universal variable:
        psi = chi0 * chi0 * alpha;

        // Compute C2 and C3:
        compute_c2c3(psi, cs); c2 = cs[0]; c3 = cs[1];

        // Propagate radial distance:
        r = chi0 * chi0 * c2 + dot(vr0, vv0) / sqrtgm * chi0 * (1.0 - psi * c3) + r0 * (1 - psi * c2);

        // Auxiliary variable for f and g functions:
        chi = chi0 + (sqrtgm * dt - pow(chi0, 3.0) * c3 - dot(vr0, vv0) / sqrtgm * pow(chi0, 2.0) * c2 - r0 * chi0 * (1.0 - psi * c3)) / r;

        // Convergence:
        if (fabs(chi - chi0) < TOL){
            break;
        }

        chi0 = chi;
    }

    if (fabs(chi - chi0) > TOL) {
#ifdef LONGDOUBLE
        printf("WARNING: failed to solver Kepler's equation, error = %23.15Lg\n", fabs(chi - chi0));
#else
        printf("WARNING: failed to solver Kepler's equation, error = %23.15g\n", fabs(chi - chi0));
#endif
    }

    // Compute f and g functions, together with their derivatives:
    real f  = 1.0 - chi * chi / r0 * c2;
    real g  = dt - pow(chi, 3.0) / sqrtgm * c3;
    real dg = 1.0 - chi * chi / r * c2;
    real df = sqrtgm / (r * r0) * chi * (psi * c3 - 1.0);

    // Propagate states:
    for (size_t d = 0; d< 3; d++) {
        x_star[particle_id][d] = f * vr0[d] + g * vv0[d];
        v_star[particle_id][d] = df * vr0[d] + dg * vv0[d];
    }

    return;
}


