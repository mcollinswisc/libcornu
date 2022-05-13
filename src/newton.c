#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "newton.h"

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

double clamped_newton(
  double x, double xmin, double xmax, double epsilon, int maxiter,
  double (*func)(void *user_data, double x),
  double (*dfunc)(void *user_data, double x),
  void *user_data)
{
  double f, df;
  int iter;

  f = func(user_data, x);
  df = dfunc(user_data, x);

  for(iter = 0; iter < maxiter; iter++) {
    // Newton iteration
    assert(df != 0);
    x -= f / df;

    // Clamp iterate
    x = MAX(xmin, MIN(xmax, x));

    // Evaulate function
    f = func(user_data, x);
    df = dfunc(user_data, x);

    if(fabs(f) <= epsilon) {
      break;
    }
  }

  return x;
}

int sign(double x)
{
  if(x < 0) {
    return -1;
  }
  else if(x > 0) {
    return 1;
  }
  else {
    return 0;
  }
}

extern int fresnel_counter;

double newton_bisection(
  double xmin, double xmax, double epsilon,
  double (*fdf)(void *user_data, double x, double *df),
  void *user_data)
{
  double fxmin, fxmax, x, fx, dfx, junk;
  int iter;

  fresnel_counter = 0;
  iter = 0;

  // Initial evalution
  fx = fxmin = fdf(user_data, x = xmin, &dfx);
  fxmax = fdf(user_data, xmax, &junk);

  while(fabs(fx) > epsilon) {
    //printf("Iter %d [%g, %g]\n", iter, xmin, xmax);
    //printf("  f(%g) = %g, f(%g) = %g\n", xmin, fxmin, xmax, fxmax);

    // Calculate Newton iterate
    //printf("Newton result: %g = %g - %g/%g\n", x - fx/dfx, x, fx, dfx);
    x -= fx/dfx;
    
    // Reject Newton iterate and go with bisection in this case
    //if(isnan(x) || fabs(dfx/fx) < 1e-10 || x <= xmin || x >= xmax) {
      x = (xmax + xmin)/2;
      //}

    // Update bounds
    fx = fdf(user_data, x, &dfx);

    if(sign(fxmin) == sign(fx)) {
      xmin = x;
      fxmin = fx;
    }
    else {
      xmax = x;
      fxmax = fx;
    }

    if(iter > 128) {
      break;
    }

    iter++;
  }

  //printf("Newton result f(%g) = %g\n", x, fx);
  //printf("  %d Iterations, %d fresnels\n", iter, fresnel_counter);
  return x;
}
