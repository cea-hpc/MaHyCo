#include <array>                  // for array, array<>::value_type
#include "EucclhydRemap.h"        // for EucclhydRemap
#include "types/MathFunctions.h"  // for det

RealArray2D<2, 2> EucclhydRemap::inverse(RealArray2D<2, 2> a) {
  double alpha = 1.0 / MathFunctions::det(a);
  return {{{{a[1][1] * alpha, -a[0][1] * alpha}},
           {{-a[1][0] * alpha, a[0][0] * alpha}}}};
}

double EucclhydRemap::divideNoExcept(double a, double b) {
  if (MathFunctions::fabs(b) < 1.0E-12)
    return 0.0;
  else
    return a / b;
}

double EucclhydRemap::crossProduct2d(RealArray1D<2> a, RealArray1D<2> b) {
  return a[0] * b[1] - a[1] * b[0];
}
