#include <cstdlib>
#include <ctime>
#include <cmath>

#include <cairo.h>

#include "timing.h"
#include "cornu.h"
#include "../external/kimia/euler.hh"

typedef struct {
  double P[2];
  double theta;

} endpoint;

#define N (4096)
endpoint e1[N], e2[N];

static double frand(double a, double b)
{
  return (rand() / (double)RAND_MAX) * (b - a) + a;
}

// Sample a random cornu endpoint
static endpoint rand_endpoint()
{
  endpoint e;

  e.P[0] = frand(0, 1);
  e.P[1] = frand(0, 1);
  e.theta = frand(0, 2 * M_PI);

  return e;
}


int main(int argc, char **argv)
{
  cornu_spiral cornu;
  unsigned int i, seed;

  struct timeval timer;
  double duration;

  // seed = 0
  seed = time(NULL);
  //seed = 1328042666;
  printf("seed = %d\n", seed);
  srand(seed);

  for(i = 0; i < N; i++) {
    e1[i] = rand_endpoint();
    e2[i] = rand_endpoint();
  }

  timer_start(&timer);
  for(i = 0; i < N; i++) {
    cornu_spiral_angles(&cornu, e1[i].P, e1[i].theta, e2[i].P, e2[i].theta);
  }
  duration = timer_interval(&timer);

  printf("Walton: %d spirals in ", N);
  time_print(stdout, duration);
  printf(" (");
  time_print(stdout, duration/N);
  printf(" per spiral)\n");

  timer_start(&timer);
  for(i = 0; i < N; i++) {
    EulerSpiral spiral(Point2D<double>(e1[i].P[0], e1[i].P[1]), e1[i].theta,
		       Point2D<double>(e2[i].P[0], e2[i].P[1]), e2[i].theta);

    
  }
  duration = timer_interval(&timer);

  printf("Kimia: %d spirals in ", N);
  time_print(stdout, duration);
  printf(" (");
  time_print(stdout, duration/N);
  printf(" per spiral)\n");

  return 0;
}
