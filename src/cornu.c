#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

#include "cornu.h"
#include "fresnel.h"
#include "integration.h"
#include "newton.h"

#define COMPARISON_EPSILON (1e-4) // TODO: replace

#define FRESNEL_EPSILON (1e-8)
#define NEWTON_EPSILON (1e-6)

static bool numeric_equal(double a, double b)
{
  return fabs(a - b) / (fabs(a) + fabs(b)) < COMPARISON_EPSILON;
}

// Shorthand names for Fresnel integrals
static inline double S(double x)
{
  return fresnel_sine_fleckner(x, FRESNEL_EPSILON);
}

static inline double C(double x)
{
  return fresnel_cosine_fleckner(x, FRESNEL_EPSILON);
}

static double dS(double x)
{
  if(fabs(x) < 1e-10) {
    return 0;
  }

  return 1/sqrt(2*M_PI) * sin(x) / sqrt(x);
}

static double dC(double x)
{
  if(fabs(x) < 1e-10) {
    return 0;
  }

  return 1/sqrt(2*M_PI) * cos(x) / sqrt(x);
}


/**
 * Function to solve for C-shaped case.
 *
 * @see newton_bisection
 */
static double fdf(void *user_data, double theta, double *out_df)
{
  double f, df, phi1, phi2, Ssum, Stheta, Csum, Ctheta;
  const cornu_spiral *cornu;

  cornu = (const cornu_spiral *)user_data;
  phi1 = cornu->phi1;
  phi2 = cornu->phi2;

  // Initialize
  f = df = 0;

  // Add terms with S integral
  Ssum = S(theta + phi1 + phi2);
  Stheta = S(theta);

  f += (Ssum - Stheta) * cos(theta + phi1);
  df += (dS(theta + phi1 + phi2) - dS(theta)) * cos(theta + phi1);
  df -= (Ssum - Stheta) * sin(theta + phi1);

  // Add terms with C integral
  Csum = C(theta + phi1 + phi2);
  Ctheta = C(theta);

  f -= (Csum - Ctheta) * sin(theta + phi1);
  df -= (dC(theta + phi1 + phi2) - dC(theta)) * sin(theta + phi1);
  df -= (Csum - Ctheta) * cos(theta + phi1);

  // TODO
  f *= sqrt(2*M_PI);
  df *= sqrt(2*M_PI);

  *out_df = df;
  return f;
}

static void arc_shaped(cornu_spiral *cornu, double D, double phi1, double phi2);

static void Cshaped(cornu_spiral *cornu,
  double D, double phi1, double phi2)
{
  double lambda, theta0, theta, a;

  //printf("C-shaped  phi1 = %g, phi2 = %g\n", phi1, phi2);

  // See: Theorem 2 of Walton & Meek 2009
  lambda = (1 - cos(phi1)) / (1 - cos(phi2));
  theta0 = (lambda * lambda) / (1 - lambda * lambda) * (phi1 + phi2);
  //printf("lambda = %lg\ntheta0 = %lg\n", lambda, theta0);

  if(theta0 > fresnel_max_fleckner()) {
    arc_shaped(cornu, D, phi1, phi2);
    return;
  }

  cornu->theta_min = theta = newton_bisection(0, theta0, NEWTON_EPSILON, fdf, cornu);
  cornu->theta_max = theta + phi1 + phi2;

  /*
  printf("theta = %lg\n", theta);
  printf("theta + phi1 + phi2 = %lg\n", theta + phi1 + phi2);

  printf("S(theta) = %lg\n", S(theta));
  printf("S(theta + phi1 + phi2) = %lg\n", S(theta + phi1 + phi2));
  printf("C(theta) = %lg\n", C(theta));
  printf("C(theta + phi1 + phi2) = %lg\n", C(theta + phi1 + phi2));
  */

  a = (S(theta + phi1 + phi2) - S(theta)) * sin(theta + phi1);
  a += (C(theta + phi1 + phi2) - C(theta)) * cos(theta + phi1);
  assert(fabs(a) > 1e-8);
  a = D / a;
  cornu->a = a;
}

/**
 * Function to solve for S-shaped case.
 *
 * @see newton_bisection
 */
static double gdg(void *user_data, double omega, double *out_dg)
{
  double g, dg, phi1, phi2, Ssum, Somega, Csum, Comega;
  const cornu_spiral *cornu;

  cornu = (const cornu_spiral *)user_data;
  phi1 = cornu->phi1;
  phi2 = cornu->phi2;

  // Initialize
  g = dg = 0;

  // Add terms with S integral
  Ssum = S(omega + phi1 + phi2);
  Somega = S(omega);

  g += (Ssum + Somega) * cos(omega + phi1);
  dg += (dS(omega + phi1 + phi2) - dS(omega)) * cos(omega + phi1);
  dg -= (Ssum + Somega) * sin(omega + phi1);

  // Add terms with C integral
  Csum = C(omega + phi1 + phi2);
  Comega = C(omega);

  g -= (Csum + Comega) * sin(omega + phi1);
  dg -= (dC(omega + phi1 + phi2) + dC(omega)) * sin(omega + phi1);
  dg -= (Csum + Comega) * cos(omega + phi1);

  // TODO
  g *= sqrt(2*M_PI);
  dg *= sqrt(2*M_PI);

  *out_dg = dg;
  return g;
}

static void Sshaped(cornu_spiral *cornu,
  double D, double phi1, double phi2)
{
  double omega, omega_min, omega_max, a, phis[2];

  //printf("S-shaped  phi1 = %g, phi2 = %g\n", phi1, phi2);

  // Find bounds
  if(phi1 > 0) {
    omega_min = 0;
  }
  else {
    omega_min = -phi1;
  }
  omega_max = M_PI/2 - phi1;

  // See: Theorem 3 & 4 of Walton & Meek 2009
  omega = newton_bisection(omega_min, omega_max, NEWTON_EPSILON, gdg, cornu);

  cornu->theta_min = -omega;
  cornu->theta_max = omega + phi1 + phi2;

  a = (S(omega + phi1 + phi2) + S(omega)) * sin(omega + phi1);
  a += (C(omega + phi1 + phi2) + C(omega)) * cos(omega + phi1);
  a = D / a;
  cornu->a = a;
}

static void rotate(double angle, const double *in, double *out)
{
  double s, c;

  s = sin(angle);
  c = cos(angle);

  out[0] = c * in[0] - s * in[1];
  out[1] = s * in[0] + c * in[1];
}

static void reflect_vector(const double across[2], double v[2])
{
  double dot, proj[2];
  int i;

  dot = 0;
  for(i = 0; i < 2; i++) {
    dot += across[i] * v[i];
  }
  
  for(i = 0; i < 2; i++) {
    proj[i] = dot * across[i];
    v[i] = 2 * proj[i] - v[i];
  }
}

/**
 * Relect a 2D point p across the line p0 + t * m
 */
static void reflect_point(const double p0[2], const double m[2], double p[2])
{
  double dp[2];
  int i;

  for(i = 0; i < 2; i++) {
    dp[i] = p[i] - p0[i];
  }
  reflect_vector(m, dp);
  for(i = 0; i < 2; i++) {
    p[i] = p0[i] + dp[i];
  }
}

static void normalize(double p[2])
{
  double mag;
  int i;
  
  mag = 0;
  for(i = 0; i < 2; i++) {
    mag += p[i] * p[i];
  }
  mag = sqrt(mag);

  for(i = 0; i < 2; i++) {
    p[i] /= mag;
  }
}

/**
 * Calculate the offset from the inflection point for a given parameter theta
 */
static void cornu_spiral_offset(const cornu_spiral *cornu, double theta, double Q[2], double epsilon)
{
  double Ctheta, Stheta, a, N0[2], Dvec[2];
  int i;

  //if(cornu->reflect) {
  //  theta = cornu->theta-theta;
  //}

  a = cornu->a;
  N0[0] = -cornu->T0[1];
  N0[1] = cornu->T0[0];

  epsilon = FRESNEL_EPSILON;
  Ctheta = fresnel_cosine_fleckner(fabs(theta), epsilon);
  Stheta = fresnel_sine_fleckner(fabs(theta), epsilon);

  for(i = 0; i < 2; i++) {
    Q[i] = a * Ctheta * cornu->T0[i] + a * Stheta * N0[i];
  }

  if(theta < 0) {
    for(i = 0; i < 2; i++) {
      Q[i] = -Q[i];
    }
  }

  if(cornu->reflect) {
    Dvec[0] = cornu->P2[0] - cornu->P1[0];
    Dvec[1] = cornu->P2[1] - cornu->P1[1];
    normalize(Dvec);

    reflect_point(cornu->P1, Dvec, Q);
  }
}

static void arc_shaped(cornu_spiral *cornu,
  double D, double phi1, double phi2)
{
  double d1[2], d2[2], r, c[2], swap;
  int d;

  //fprintf(stderr, "Arc case\n");
  cornu->arc = 1;

  // Calculate circle parameters
  cornu->arc_r = fabs(r = D/(2 * sin(phi1)));
  cornu->arc_center[0] = c[0] = cornu->P1[0] - r * cornu->T1[1];
  cornu->arc_center[1] = c[1] = cornu->P1[1] + r * cornu->T1[0];

  // Find angle limits
  for(d = 0; d < 2; d++) {
    d1[d] = cornu->P1[d] - c[d];
    d2[d] = cornu->P2[d] - c[d];
  }

  cornu->theta_min = atan2(d1[1], d1[0]);
  cornu->theta_max = atan2(d2[1], d2[0]);

  if(cornu->theta_max < cornu->theta_min) {
    swap = cornu->theta_min;
    cornu->theta_min = cornu->theta_max;
    cornu->theta_max = swap;
  }

  if(cornu->theta_max < 0) {
    cornu->theta_min += 2 * M_PI;
    cornu->theta_max += 2 * M_PI;
  }
  assert(cornu->theta_max > 0);
}

static void cornu_spiral_phi(cornu_spiral *cornu,
  double D, double phi1, double phi2)
{
  int i;
  double swap, h, pimult, P[2], T[2], Dvec[2];

  //printf("D = %lg\nphi1 = %lg\nphi2 = %lg\n", D, phi1, phi2);

  // Check for degenerate cases
  cornu->line = 0;
  cornu->arc = 0;

  //printf("arc criterion: %lg\n", fabs(phi1 - phi2));

  if(fabs(fabs(phi1) - fabs(phi2)) < COMPARISON_EPSILON) {
    if(fabs(phi1) < COMPARISON_EPSILON) {
      cornu->line = 1;
      return;
    }

    if(fabs(phi1 - phi2) < COMPARISON_EPSILON) {
      arc_shaped(cornu, D, phi1, phi2);
      return;
    }

    //fprintf(stderr, "%s:%d Bad degenerate case (TODO?) phi1 = %g, phi2 = %g\n", __FILE__, __LINE__, phi1, phi2);
  }

  cornu->reflect = 0;
  cornu->reverse = 0;

  // Check for transformations
  // TODO: weirdness when both happens...
  if(fabs(phi1) > fabs(phi2)) {
    cornu->reverse = 1;

    swap = phi1;
    phi1 = -phi2;
    phi2 = -swap;

    for(i = 0; i < 2; i++) {
      swap = cornu->P1[i];
      cornu->P1[i] = cornu->P2[i];
      cornu->P2[i] = swap;

      swap = cornu->T1[i];
      cornu->T1[i] = -cornu->T2[i];
      cornu->T2[i] = -swap;
    }

    cornu->Tangle1 = atan2(cornu->T1[1], cornu->T1[0]);
    cornu->Tangle2 = atan2(cornu->T2[1], cornu->T2[0]);
  }

  if(phi2 > M_PI) {
    phi1 -= 2*M_PI;
    phi2 -= 2*M_PI;
    cornu->reflect = !cornu->reflect;
  }
  if(phi2 < 0) {
    phi1 = -phi1;
    phi2 = -phi2;
    cornu->reflect = !cornu->reflect;
  }

  assert(fabs(phi1) < fabs(phi2) + 1e-8);
  assert(0 < phi2 && phi2 <= M_PI);

  // Figure out C-shaped vs S-shaped case
  h = S(phi1 + phi2) * cos(phi1) - C(phi1 + phi2) * sin(phi1);

  //printf("D = %lg\nphi1 = %lg = %lg deg\nphi2 = %lg = %lg deg\nh = %lg\n",
  //       D, phi1, phi1 * 180/M_PI, phi2, phi2 * 180/M_PI, h);
  //if(cornu->reverse) {
  //  printf("reverse\n");
  //}
  //if(cornu->reflect) {
  //  printf("reflect\n");
  //}

  cornu->phi1 = phi1;
  cornu->phi2 = phi2;

  if(0 < phi1 && h <= 0) {
    cornu->shape = 'C';
    Cshaped(cornu, D, phi1, phi2);
    if(cornu->arc) {
      return;
    }
  }
  else {
    cornu->shape = 'S';
    Sshaped(cornu, D, phi1, phi2);
  }

  // Calculate inflection point parameters
  P[0] = cornu->P1[0];
  P[1] = cornu->P1[1];

  T[0] = cornu->T1[0];
  T[1] = cornu->T1[1];

  if(cornu->reflect) {
    Dvec[0] = cornu->P2[0] - cornu->P1[0];
    Dvec[1] = cornu->P2[1] - cornu->P1[1];
    normalize(Dvec);

    reflect_vector(Dvec, T);
  }

  if(cornu->shape == 'C') {
    rotate(-cornu->theta_min, T, cornu->T0);
  }
  else {
    rotate(cornu->theta_min, T, cornu->T0);
  }

  cornu_spiral_offset(cornu, cornu->theta_min, cornu->P0, cornu->a * 1e-4);
  for(i = 0; i < 2; i++) {
    cornu->P0[i] = P[i] - cornu->P0[i];
  }
}

void cornu_spiral_tangents(cornu_spiral *cornu, const double P1[2], const double T1[2], const double P2[2], const double T2[2])
{
  int i;
  double a1, a2;

  for(i = 0; i < 2; i++) {
    cornu->P1[i] = P1[i];
    cornu->T1[i] = T1[i];

    cornu->P2[i] = P2[i];
    cornu->T1[i] = T1[i];
  }

  cornu->Tangle1 = a1 = atan2(T1[1], T1[0]);
  cornu->Tangle2 = a2 = atan2(T2[1], T2[0]);

  cornu_spiral_angles(cornu, P1, a1, P2, a2);
}

static double phirange(double phi)
{
  phi = fmod(phi, 2*M_PI);

  if(phi > M_PI) {
    phi -= 2 * M_PI;
  }
  if(phi < -M_PI) {
    phi += 2 * M_PI;
  }

  return phi;
}

void cornu_spiral_angles(cornu_spiral *cornu, const double P1[2], double a1, const double P2[2], double a2)
{
  int i;
  double Dx, Dy, D, Dangle, phi1, phi2;

  for(i = 0; i < 2; i++) {
    cornu->P1[i] = P1[i];
    cornu->P2[i] = P2[i];
  }

  cornu->T1[0] = cos(a1);
  cornu->T1[1] = sin(a1);
  cornu->Tangle1 = a1;

  cornu->T2[0] = cos(a2);
  cornu->T2[1] = sin(a2);
  cornu->Tangle2 = a2;

  // Calculate input invariants
  Dx = P2[0] - P1[0];
  Dy = P2[1] - P1[1];

  D = hypot(Dy, Dx);
  Dangle = atan2(Dy, Dx);

  phi1 = phirange(Dangle - a1);
  phi2 = phirange(a2 - Dangle);

  //printf("D = %lg\nphi1 = %lg\nphi2 = %lg\n", D, phi1, phi2);

  cornu_spiral_phi(cornu, D, phi1, phi2);
}

void cornu_spiral_points(const cornu_spiral *cornu, double epsilon, int N, double *p)
{
  int i, k;
  double theta, dtheta, P0[2];

  dtheta = (cornu->theta_max - cornu->theta_min) / (N - 1);

  // Special cases
  if(cornu->line) {
    for(i = 0; i < N; i++) {
      theta = ((double) i) / (N-1);
      for(k = 0; k < 2; k++) {
        p[2*i + k] = (1 - theta) * cornu->P1[k] + theta * cornu->P2[k];
      }
    }
    return;
  }
  if(cornu->arc) {
    //fprintf(stderr, "Plotting arc case\n");
    for(i = 0; i < N; i++) {
      theta = i * dtheta + cornu->theta_min;
      p[2*i] = cornu->arc_center[0] + cornu->arc_r * cos(theta);
      p[2*i+1] = cornu->arc_center[1] + cornu->arc_r * sin(theta);
    }
    return;
  }

  // Real cornu spiral
  for(i = 0; i < N; i++) {
    theta = i * dtheta + cornu->theta_min;

    cornu_spiral_offset(cornu, theta, p + 2*i, epsilon);
    for(k = 0; k < 2; k++) {
      p[2*i + k] += cornu->P0[k];
    }
  }
}
