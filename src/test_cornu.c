#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <cairo.h>
#include <cairo-pdf.h>

#include "timing.h"
#include "cornu.h"

static const int width = 640;
static const int height = 320;

typedef struct {
  double P[2];
  double theta;

} endpoint;

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

static void plot(const char *path, const cornu_spiral *cornu);

int main(int argc, char **argv)
{
  endpoint e1, e2;
  cornu_spiral cornu;
  unsigned int seed;

  if(argc > 1) {
    if(strcmp(argv[1], "time") == 0) {
      seed = time(NULL);
    }
    else {
      seed = atoi(argv[1]);
    }
  }
  else {
    seed = 0;
    //seed = time(NULL);
    //seed = 1328042666;
  }

  printf("seed = %d\n", seed);
  srand(seed);

  e1 = rand_endpoint();
  e2 = rand_endpoint();

  //e1.P[0] = 0;  e1.P[1] = 0; e1.theta = -M_PI/6;
  //e2.P[0] = 3;  e2.P[1] = 0; e2.theta = -M_PI/2;

  // Failure case?
  //e1.P[0] = 215.1339038;  e1.P[1] = 119.9445352; e1.theta = -1.911657459;
  //e2.P[0] = 232.8759227;  e2.P[1] = 211.0513945; e2.theta = -1.626260606;

  //e1.P[0] = -1;  e1.P[1] = 0; e1.theta = M_PI/2;
  //e2.P[0] =  1;  e2.P[1] = 0; e2.theta = M_PI;

  printf("%g %g %g\n", e1.P[0], e1.P[1], e1.theta);
  printf("%g %g %g\n", e2.P[0], e2.P[1], e2.theta);

  // Calculate spiral
  cornu_spiral_angles(&cornu, e1.P, e1.theta, e2.P, e2.theta);

  // Print useful info
  printf("phi1 = %lg\n", cornu.phi1);
  printf("phi2 = %lg\n", cornu.phi2);

  if(cornu.line) {
    printf("special case: line\n");
  }
  else if(cornu.arc) {
    printf("special case: arc\n");
    printf("center: (%lg, %lg)\n", cornu.arc_center[0], cornu.arc_center[1]);
    printf("radius: %lg\n", cornu.arc_r);
    printf("theta: [%lg, %lg]\n", cornu.theta_min, cornu.theta_max);
  }
  else {
    printf("shape = %c\n", cornu.shape);
    printf("theta: [%lg, %lg]\n", cornu.theta_min, cornu.theta_max);
    printf("P0: (%lg, %lg)\n", cornu.P0[0], cornu.P0[1]);
    printf("T0: (%lg, %lg)\n", cornu.T0[0], cornu.T0[1]);
    printf("a = %lg\n", cornu.a);
    if(cornu.reverse) {
      printf("reverse\n");
    }
    if(cornu.reflect) {
      printf("reflect\n");
    }
  }

  // Plot
  plot("cornu.pdf", &cornu);
  
  return 0;
}

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

static void bounds(int N, double *pts, double *min, double *max)
{
  int i, k;

  min[0] = INFINITY;
  min[1] = INFINITY;

  max[0] = -INFINITY;
  max[1] = -INFINITY;

  for(i = 0; i < N; i++) {
    for(k = 0; k < 2; k++) {
      min[k] = MIN(min[k], pts[2*i+k]);
      max[k] = MAX(max[k], pts[2*i+k]);
    }
  }
}

static double plot_setup(
  cairo_surface_t *surf, cairo_t *cr, double width, double height,
  double *inmin, double *inmax)
{
  double min[2], max[2], size[2], inner_width, inner_height, start_ratio, end_ratio, scale;
  size_t i, j;

  cairo_translate(cr, 0, height);
  cairo_scale(cr, 1, -1);

  // Adjust bounding box
  for(i = 0; i < 2; i++) {
    min[i] = inmin[i];
    max[i] = inmax[i];

    size[i] = max[i] - min[i];
    min[i] -= size[i] * 0.1;
    max[i] += size[i] * 0.1;
    size[i] = max[i] - min[i];
  }
  start_ratio = (double)width/height;
  end_ratio = size[0]/size[1];
 
  // Find right ratio
  if(end_ratio < start_ratio) {
    // height is limiting dimension
    scale = height / size[1];

    inner_width = end_ratio * height;
    cairo_translate(cr, (width - inner_width) / 2, 0);
  }
  else {
    // width is limiting dimension
    scale = width / size[0];

    inner_height = width / end_ratio;
    cairo_translate(cr, 0, (height - inner_height) / 2);
  }

  cairo_scale(cr, scale, scale);
  cairo_set_line_width(cr, cairo_get_line_width(cr) * scale);

  cairo_translate(cr, -min[0], -min[1]);

  return scale;
}

static void draw_path(cairo_t *cr, int n, const double *p)
{
  int i;

  if(n < 2) {
    return;
  }

  cairo_move_to(cr, p[0], p[1]);
  for(i = 1; i < n; i++) {
    cairo_line_to(cr, p[2*i], p[2*i+1]);
  }
  cairo_stroke(cr);
}

static void plot(const char *path, const cornu_spiral *cornu)
{
  const int N = 32;
  int i;
  double *points, min[2], max[2], scale, D;
  cairo_surface_t *surface;
  cairo_t *cr;
  

  points = (double *)malloc(sizeof(double) * 2 * (N + 2));
  cornu_spiral_points(cornu, 1e-2, N, points);

  points[2*N] = cornu->P1[0];
  points[2*N+1] = cornu->P1[1];
  points[2*N+2] = cornu->P2[0];
  points[2*N+3] = cornu->P2[1];

  for(i = 0; i < N; i++) {
    printf("(%lg, %lg)\n", points[2*i], points[2*i+1]);
  }

  // Set up drawing
  bounds(N+2, points, min, max);
  surface = cairo_pdf_surface_create(path, width, height);
  //surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
  cr = cairo_create(surface);

  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_rectangle(cr, 0, 0, width, height);
  cairo_fill(cr);

  cairo_set_line_width(cr, 0.8);

  scale = plot_setup(surface, cr, width, height, min, max);
  //printf("scale = %g\n", scale);
  cairo_set_line_width(cr, 1.0 / scale);

  // Draw the spiral segment
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_move_to(cr, points[0], points[1]);
  for(i = 1; i < N; i++) {
    cairo_line_to(cr, points[2*i + 0], points[2*i + 1]);
  }
  cairo_stroke(cr);

  // Draw endpoints
  cairo_set_source_rgb(cr, 1, 0, 0);
  cairo_move_to(cr, cornu->P1[0], cornu->P1[1]);
  cairo_line_to(cr, cornu->P1[0] + cornu->T1[0] * 25/scale, cornu->P1[1] + cornu->T1[1] * 25/scale);
  cairo_stroke(cr);

  cairo_set_source_rgb(cr, 0, 0, 0.75);
  cairo_move_to(cr, cornu->P2[0], cornu->P2[1]);
  cairo_line_to(cr, cornu->P2[0] + cornu->T2[0] * 25/scale, cornu->P2[1] + cornu->T2[1] * 25/scale);
  cairo_stroke(cr);

  // Draw connecting line
  /*
  cairo_set_source_rgb(cr, 0, 0, 0);

  double dash = 5/scale;
  cairo_set_dash(cr, &dash, 1, 0);
  
  cairo_move_to(cr, cornu->P1[0], cornu->P1[1]);
  cairo_line_to(cr, cornu->P2[0], cornu->P2[1]);
  cairo_stroke(cr);
  cairo_set_dash(cr, NULL, 0, 0);
  */

  // Write to file
  //cairo_surface_write_to_png(surface, path);

  // Cleanup
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
  free(points);
}
