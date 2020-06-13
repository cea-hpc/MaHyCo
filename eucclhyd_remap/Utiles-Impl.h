
#ifndef UTILES_IMPL_H
#define UTILES_IMPL_H

template <size_t N, size_t K>
RealArray2D<N, K> EucclhydRemap::tensProduct(RealArray1D<N> a,
                                             RealArray1D<K> b) {
  RealArray2D<N, K> res;
  for (size_t i = 0; i < N; i++)
    for (size_t j = 0; j < K; j++) res[i][j] = a[i] * b[j];
  return res;
}

#endif  // UTILES_IMPL_H