#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <gmp.h>
#include <mpfr.h>

#include "fresnel.h"

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif


double fresnel_cosine_series_high_precision(double x, double epsilon)
{
  double log_term, log_max_term;
  int n, N;

  mpfr_prec_t precision;
  mpfr_t px, fact, n2, gamma_arg, term, C;

  double result;

  // Find number of terms necessary
  log_term = INFINITY;
  log_max_term = -INFINITY;
  for(N = 0; log_term > log(epsilon); N++) {
    log_term = 2 * N * log(x) - log(2*N + 0.5) - lgamma(2*N+1);
    log_max_term = MAX(log_term, log_max_term);
  }
  printf("N = %d\n", N);
  printf("max_logterm = %lg\n", log_max_term);

  // Calculate necessary precision
  precision = (mpfr_prec_t)(log_max_term / log(2)) + 64;
  //precision *= 4;
  printf("precision = %ld\n", (long int)precision);

  // Initialize variables
  mpfr_inits2(precision, px, fact, n2, gamma_arg, term, C, NULL);

  mpfr_set_d(C, 0, MPFR_RNDZ);
  mpfr_set_d(px, x, MPFR_RNDN);

  // Calculate sum
  for(n = 0; n < N; n++) {
    mpfr_set_d(n2, 2*n, MPFR_RNDN);

    // term = 2 * n * log(x)
    mpfr_log(term, px, MPFR_RNDN);
    mpfr_mul_d(term, term, 2 * n, MPFR_RNDN);

    // term -= log(gamma(2*n+1))
    mpfr_set_d(gamma_arg, 2*n+1, MPFR_RNDN);
    mpfr_lngamma(fact, gamma_arg, MPFR_RNDN);
    mpfr_sub(term, term, fact, MPFR_RNDN);

    // term = exp(term)
    mpfr_exp(term, term, MPFR_RNDN);

    // temp /= 2*n+0.5
    mpfr_div_d(term, term, 2*n + 0.5, MPFR_RNDN);

    log_term = 2 * n * log(x) - log(2*n + 0.5) - lgamma(2*n+1);

      //mpfr_set_d(term, log_term, MPFR_RNDN);
      //mpfr_exp(term, term, MPFR_RNDN);

    // term *= (-1)^n
    if(n % 2 == 1) {
      mpfr_neg(term, term, MPFR_RNDN);
    }

    mpfr_add(C, C, term, MPFR_RNDN);

    if(n != 0 && n % 64 == 0) {
      mpfr_printf("%4d: C = %Rg, term = %Rg, log(term) = %g\n", n, C, term, log_term);
    }
  }

  result = mpfr_get_d(C, MPFR_RNDN) * sqrt(x / (2*M_PI));
  
  // Cleanup
  mpfr_clears(px, fact, n2, gamma_arg, term, C, NULL);

  return result;
}

double fresnel_sine_series_exact(double x, double epsilon)
{
  int n;
  mpq_t fact, fact_term, x2, xpow, C, term, abs_term, max_term;
  double result;

  mpq_inits(fact, fact_term, x2, xpow, C, abs_term, term, max_term, NULL);
  mpq_set_ui(fact, 1, 1);
  mpq_set_ui(xpow, 1, 1);

  // x2 = x^2
  mpq_set_d(x2, x);
  mpq_mul(x2, x2, x2);
  mpq_canonicalize(x2);

  n = 0;
  do {
    // term = 1/(2n + 1.5)
    mpq_set_ui(term, 2, 4*n+3);

    // term /= fact
    mpq_div(term, term, fact);

    // term *= xpow
    mpq_mul(term, term, xpow);

    // max_term = max(abs(term), max_term)
    if(mpq_cmp(term, max_term) > 0) {
      mpq_set(max_term, term);
    }

    // term *= (-1)^n
    if(n % 2 == 1) {
      mpq_neg(term, term);
    }

    // C += term
    mpq_add(C, C, term);


    // xpow *= x2;
    mpq_mul(xpow, xpow, x2);

    // fact *= (2*n + 2) * (2*n + 1);
    mpq_set_ui(fact_term, (2*n + 3) * (2*n + 2), 1);
    mpq_mul(fact, fact, fact_term);

    n++;
  } while(fabs(mpq_get_d(term)) > epsilon);


  result = mpq_get_d(C) * x * sqrt(x / (2*M_PI));

  // Cleanup
  mpq_clears(fact, x2, xpow, C, term, NULL);

  return result;
}

double fresnel_cosine_series_exact(double x, double epsilon)
{
  int n;
  mpq_t fact, fact_term, x2, xpow, C, term, abs_term, max_term;
  double result;
  mpf_t Cf, Cf2, termf, mul;
  mp_bitcnt_t prec;

  prec = 2306 + 64;

  if(x == 0) {
    return 0;
  }

  // Intitialize
  // fact = xpow = 1, C = 0
  mpq_inits(fact, fact_term, x2, xpow, C, abs_term, term, max_term, NULL);
  mpf_init2(termf, prec);
  mpf_init2(Cf, prec);
  mpf_init2(mul, prec);
  mpq_set_ui(fact, 1, 1);
  mpq_set_ui(xpow, 1, 1);

  // x = x^2
  mpq_set_d(x2, x);
  mpq_mul(x2, x2, x2);
  mpq_canonicalize(x2);

  //printf("%lg\n", mpq_get_d(x2));

  n = 0;
  do {
    // term = 1/(2n + 0.5)
    mpq_set_ui(term, 2, 4*n+1);

    // term /= fact
    mpq_div(term, term, fact);

    // term *= xpow
    mpq_mul(term, term, xpow);

    // max_term = max(abs(term), max_term)
    if(mpq_cmp(term, max_term) > 0) {
      mpq_set(max_term, term);
    }

    // term *= (-1)^n
    if(n % 2 == 1) {
      mpq_neg(term, term);
    }

    // C += term
    mpq_add(C, C, term);

    // Cf += float(term)
    mpf_set_q(termf, term);
    mpf_add(Cf, Cf, termf);

    //printf("termf = ");
    //mpf_out_str(stdout, 10, 10, termf);
    //printf("\n");

    //printf("Cf = ");
    //mpf_out_str(stdout, 10, 10, Cf);
    //printf(" C = ");

    // xpow *= x2;
    mpq_mul(xpow, xpow, x2);

    // fact *= (2*n + 2) * (2*n + 1);
    mpq_set_ui(fact_term, (2*n + 2) * (2*n + 1), 1);
    mpq_mul(fact, fact, fact_term);

    if(n != 0 && n % 512 == 0) {
      //printf("%4d: C = %g, term = %g\n", n, mpq_get_d(C), mpq_get_d(term));
      //gmp_printf("%4d: Cf = %Fg, term = %Fg\n", n, Cf, termf);
    }

    n++;
  } while(fabs(mpq_get_d(term)) > epsilon);

  // Print out largest term...
  /*
  mpf_set_q(termf, max_term);
  printf("Largest term: ");
  mpf_out_str(stdout, 10, 10, termf);
  printf("\n");
  */

  //printf("mpf error: %lg\n", mpq_get_d(C) - mpf_get_d(Cf));

  /*
  mpf_set_d(mul, sqrt(x / (2*M_PI)));
  mpf_mul(Cf, Cf, mul);
  printf("Cf = %lg = ", mpf_get_d(Cf));
  mpf_out_str(stdout, 10, 10, Cf);
  printf("\n");

  printf("Done in %d iters\n", n);
  */
  result = mpq_get_d(C) * sqrt(x / (2*M_PI));

  // Cleanup
  mpf_clears(Cf, termf, NULL);
  mpq_clears(fact, x2, xpow, C, term, NULL);

  return result;
}
