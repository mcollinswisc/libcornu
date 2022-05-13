#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "fresnel.h"

#include "integration.h"

int fresnel_counter = 0;

static double fS(void *user_data, double u)
{
  if(u < 1e-8) {
    return 0;
  }

  return sin(u) / sqrt(u);
}

static double fC(void *user_data, double u)
{
  if(u < 1e-8) {
    return 0;
  }

  return cos(u) / sqrt(u);
}

double fresnel_sine_integral(double theta, int N)
{
  double S;

  S = integrate_composite_simpsons(fS, 0, theta, N, NULL);

  return S / sqrt(2*M_PI);
}

double fresnel_cosine_integral(double theta, int N)
{
  double C;

  C = integrate_composite_simpsons(fC, 0, theta, N, NULL);

  return C / sqrt(2*M_PI);
}

double factorial(int n)
{
  double f;
  int i;

  f = 1;
  for(i = 1; i < n; i++) {
    f *= (double)n;
  }

  return f;
}

double fresnel_sine_series(double x, int N, double epsilon)
{
  double sgn, fact, x2, xpow, S, term;
  int n;

  if(x < 1e-10) {
    return 0;
  }

  epsilon *= sqrt((2*M_PI)/x) / x;

  term = 0;
  sgn = 1;
  xpow = 1;
  x2 = x * x;
  fact = 1;

  S = 0;
  for(n = 0; n < N; n++) {
    term = sgn * xpow / (fact * (2*n + 1.5));
    S += term;

    xpow *= x2;
    fact *= (double)(2*n + 3) * (double)(2*n + 2);
    sgn = -sgn;

    //printf("S = %g\n", S * x * sqrt(x / (2*M_PI)));

    if(fabs(term) < epsilon) {
      break;
    }
  }
  //printf("fresnelS: iters = %d, delta = |%g| < %g, x =  %g\n", n, term, epsilon, x);
  if(fabs(term) > epsilon) {
    fprintf(stderr, "Error: fresnel integral S(%g) failed to converge in %d iterations.", x, n);
  }

  S *= x * sqrt(x / (2*M_PI));

  fresnel_counter++;
  return S;
}

double fresnel_cosine_series(double x, int N, double epsilon)
{
  double sgn, fact, x2, xpow, C, term;
  int n;

  if(x < 1e-10) {
    return 0;
  }

  epsilon *= sqrt((2*M_PI)/x);

  term = 0;
  sgn = 1;
  xpow = 1;
  x2 = x * x;
  fact = 1;

  C = 0;
  for(n = 0; n < N; n++) {
    term = sgn * xpow / (fact * (2*n + 0.5));
    C += term;

    xpow *= x2;
    fact *= (double)(2*n + 2) * (double)(2*n + 1);
    sgn = -sgn;

    //printf("C = %g\n", C);

    if(fabs(term) < epsilon) {
      break;
    }
  }
  //printf("fresnelC: iters = %d, delta = |%g| < %g, x =  %g\n", n, term, epsilon, x);
  if(fabs(term) > epsilon) {
    fprintf(stderr, "Error: fresnel integral C(%g) failed to converge in %d iterations.", x, n);
  }

  C *= sqrt(x / (2*M_PI));

  fresnel_counter++;
  return C;
}
