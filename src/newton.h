#ifndef __NEWTON_H__
#define __NEWTON_H__

/**
 * Solve $f(x) = 0$ where the solution is known to be in [xmin, xmax]
 *
 * Finds a root of a one-dimensional function with derivative provided.
 * Clamps iterates to the range [xmin, xmax].  Halts when the function
 * magnitude is less than epsilon, or maxiter has been reached.
 * 
 */
double clamped_newton(
  double x0, double xmin, double xmax, double epsilon, int maxiter,
  double (*f)(void *user_data, double x),
  double (*df)(void *user_data, double x),
  void *user_data);

/**
 * Solve $f(x) = 0$ where the solution is known to be in [xmin, xmax]
 *
 * Finds a root of a one-dimensional function with derivative provided.
 * Clamps iterates to the range [xmin, xmax].  Works similar to bisection
 * method, except that it used Newton's method to generate the midpoint
 * where possible.  Returns a point such that |f(x)| < epsilon.
 *
 * @arg The function to solve.  Should return the function value f(x),
 *      and place the derivative in df
 */
double newton_bisection(
  double xmin, double xmax, double epsilon,
  double (*fdf)(void *user_data, double x,  double *df),
  void *user_data);


#endif // #ifndef __NEWTON_H__
