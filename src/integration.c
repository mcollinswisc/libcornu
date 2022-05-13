#include "integration.h"

double integrate_composite_simpsons(
  double (*f)(void *user_data, double x),
  double a, double b,
  int N,
  void *user_data)
{
	double I, even, odd, h;
	int i;
	
	h = (b - a)/N;

  // Calculate integral by composite Simpson's rule
  // Even sum
	even = 0;
	for(i = 0; i < N/2 - 1; i++) {
		even += (*f)(user_data, a + (2*i) * h);
	}
	
	// Odd sum
	odd = 0;
	for(i = 1; i < N/2; i++) {
		odd += (*f)(user_data, a + (2*i - 1) * h);
	}
	
	// Endpoints
	I = (*f)(user_data, a) + 2*even + 4*odd + (*f)(user_data, b);
	I *= h/3;
	
	return I;
}
