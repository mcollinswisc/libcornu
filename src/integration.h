#ifndef __INTEGRATION_H__
#define __INTEGRATION_H__



double integrate_composite_simpsons(
  double (*f)(void *user_data, double x),
  double a, double b,
  int N,
  void *user_data);

#endif // #ifndef __INTEGRATION_H__
