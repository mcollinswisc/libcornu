#ifndef __CORNU_H__
#define __CORNU_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double P1[2];
  double T1[2];
  double Tangle1;

  double P2[2];
  double T2[2];
  double Tangle2;

  double phi1, phi2;

  //// Unknown terms

  /**
   * Whether "C-shaped" or "S-shaped"
   */
  char shape;

  /**
   * Parameter for the start of the segment.
   */
  double theta_min;

  /**
   * Parameter for the end of the segment.
   */
  double theta_max;

  /**
   * Scale factor of the cornu spiral.
   */
  double a;

  /**
   * Location + tangent of the inflection point.
   */
  double P0[2], T0[2];

  // Modifications
  int reverse;
  int reflect;

  // Degenerate cases
  int line;
  int arc;
  double arc_center[2], arc_r;

} cornu_spiral;

/**
 * Calculate Cornu spiral interpolation from tangent coordinates
 */
void cornu_spiral_tangents(cornu_spiral *C,
  const double P1[2], const double T1[2],
  const double P2[2], const double T2[2]);

/**
 * Calculate Cornu spiral interpolation given angle of tangents from x-axis.
 */
void cornu_spiral_angles(cornu_spiral *C,
  const double P1[2], double Tangle1,
  const double P2[2], double Tangle2);

/**
 * Find a sequence of points for the cornu spiral
 */
  void cornu_spiral_points(const cornu_spiral *C, double epsilon, int N, double *p);

#ifdef __cplusplus
}
#endif

#endif // #ifndef __CORNU_H__
