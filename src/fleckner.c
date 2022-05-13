#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "fresnel.h"
#include "fresnel_table.h"

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

static double Kseries(int N, double rho, double gamma)
{
  double S, K0, K1, Kn, Kn1, Kn2;
  int n;

  double an, cn, dn;

  //printf("rho = %lg, gamma = %lg\n", rho, gamma);

  if(N < 0) {
    return 0;
  }

  // Base cases
  S = K0 = 1 - cos(rho);
  if(N == 0) {
    return S;
  }
  //printf("K0 = %lg, S = %lg\n", K0, S);


  K1 = sin(rho) / gamma - (rho/gamma) * cos(rho);
  S -= 0.5 * K1;
  if(N == 1) {
    return S;
  }
  //printf("K1 = %lg, S = %lg\n", K1, S);

  Kn2 = K0;
  Kn1 = K1;

  for(n = 2; n < N; n++) {
    // Calculate K_n
    cn = -n*(n-1)/(gamma * gamma);
    dn = -pow(rho/gamma, n) * cos(rho) + (n/gamma) * pow(rho/gamma, n-1) * sin(rho);
    Kn = cn * Kn2 + dn;

    // Calculate a_n
    an = exp(lgamma(n + 0.5) - lgamma(n+1));
    an /= sqrt(M_PI);
    if(n % 2 == 1) {
      an = -an;
    }

    // Accumulate sum
    S += an * Kn;
    //printf("K%d = %lg, an = %lg, S = %lg\n", n, Kn, an, S);

    // Shift K's
    Kn2 = Kn1;
    Kn1 = Kn;
  }

  return S;
}

double fresnel_max_fleckner()
{
  return M_PI * FRESNEL_TABLE_ENTRIES;
}

double fresnel_sine_fleckner(double x1, double epsilon)
{
  double phim, m, rho, gamma;
  double N;
  double Sphim, sgn;
  int im;

  if(x1 <= 6*M_PI + 1e-3) {
    return fresnel_sine_series(x1, 100, epsilon);
  }
  
  // Split input x1
  rho = modf(x1/M_PI, &m) * M_PI;
  gamma = phim = M_PI * m;
  im = (int)m;

  assert(fabs((phim + rho) - x1) < 1e-8);
  assert(rho - 1e-8 < M_PI);
  assert(fabs(phim - im * M_PI) < 1e-8);

  if(im >= FRESNEL_TABLE_ENTRIES) {
    Sphim = 0.5; // Assume approximately convergent...
    return 0.5; // TODO: Do bettter?
  }
  else {
    Sphim = fresnel_table_S[im];
  }

  // Determine # of series terms needed
  N = ceil(log(M_PI/epsilon) / log(m));
  N = MAX(N, 5);
  //printf("N = %lg\n", N);

  // Calculate
  sgn = (im % 2 == 0) ? 1.0 : -1.0;
  return Sphim + sgn / sqrt(2 * M_PI * phim) * Kseries((int)N, rho, gamma);
}

double fresnel_cosine_fleckner(double x0, double epsilon)
{
  double thetak, k, rho, gamma;
  double N;
  double Cthetak, sgn;
  int ik;

  if(x0 <= 6.5*M_PI + 1e-3) {
    return fresnel_cosine_series(x0, 100, epsilon);
  }
  
  // Split input x0
  k = floor(x0/M_PI - 0.5);
  ik = (int)k;

  gamma = thetak = (2 * k + 1) * M_PI/2;
  rho = x0 - thetak;

  assert(fabs((thetak + rho) - x0) < 1e-8);
  assert(rho - 1e-8 < M_PI);
  //assert(fmod(thetak - M_PI/2, M_PI) < 1e-8);
  assert(fabs((2*ik + 1) * M_PI/2 - thetak) < 1e-8);

  // TODO: use table
  if(ik >= FRESNEL_TABLE_ENTRIES) {
    Cthetak = 0.5; // Assume approximately convergent...
    return 0.5; // TODO: Do bettter?
  }
  else {
    Cthetak = fresnel_table_C[ik];
  }

  // Determine # of series terms needed
  N = ceil(log(M_PI/epsilon) / log(k + 0.5));
  N = MAX(N, 5);
  //printf("N = %lg\n", N);

  // Calculate
  sgn = ((ik + 1) % 2 == 0) ? 1.0 : -1.0;
  return Cthetak + sgn / sqrt(2 * M_PI * thetak) * Kseries((int)N, rho, gamma);
}
