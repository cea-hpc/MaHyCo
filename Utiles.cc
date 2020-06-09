

KOKKOS_INLINE_FUNCTION
RealArray2D<2, 2> inverse(RealArray2D<2, 2> a) {
  double alpha = 1.0 / MathFunctions::det(a);
  return {{{{a[1][1] * alpha, -a[0][1] * alpha}},
           {{-a[1][0] * alpha, a[0][0] * alpha}}}};
}

KOKKOS_INLINE_FUNCTION
double divideNoExcept(double a, double b) {
  if (MathFunctions::fabs(b) < 1.0E-12)
    return 0.0;
  else
    return a / b;
}

template <size_t N, size_t M>
KOKKOS_INLINE_FUNCTION RealArray2D<N, M> tensProduct(RealArray1D<N> a,
                                                     RealArray1D<M> b) {
  RealArray2D<N, M> res;
  for (size_t i = 0; i < N; i++)
    for (size_t j = 0; j < M; j++) res[i][j] = a[i] * b[j];
  return res;
}

KOKKOS_INLINE_FUNCTION
double crossProduct2d(RealArray1D<2> a, RealArray1D<2> b) {
  return a[0] * b[1] - a[1] * b[0];
}
