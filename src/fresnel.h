#ifndef __FRESNEL_H__
#define __FRESNEL_H__

double fresnel_sine_integral(double theta, int N);
double fresnel_cosine_integral(double theta, int N);

double fresnel_sine_series(double theta, int N, double epsilon);
double fresnel_cosine_series(double theta, int N, double epsilon);


/**
 * Compute the Fresnel sine integral using  exact arithmetic via power series
 */
double fresnel_sine_series_exact(double theta, double epsilon);

/**
 * Compute the Fresnel cosine integral using  exact arithmetic via power series
 */
double fresnel_cosine_series_exact(double theta, double epsilon);

double fresnel_cosine_series_high_precision(double x, double epsilon);


/**
 * The maximum range for fleckner integrals
 */
double fresnel_max_fleckner();

/**
 * Compute the Fresnel sine integral using the method of Fleckner 1968
 */
double fresnel_sine_fleckner(double theta, double epsilon);

/**
 * Compute the Fresnel cosine integral using the method of Fleckner 1968
 */
double fresnel_cosine_fleckner(double theta, double epsilon);

// TODO: incremental integrals

#endif // #ifndef __FRESNEL_H__
