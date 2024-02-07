#ifndef ACCENV_SORT_SMALL_H
#define ACCENV_SORT_SMALL_H

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;

namespace Accenv {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*!
 * \brief Positionne la valeur minimale au début du tableau à l'intérieur 
 *        d'un kernel, le reste des éléments pouvant être permutés
 *
 * @param[inout] arr : [in] tableau à modifier [out] tableau avec valeur minimale en premier
 *
 */
/*---------------------------------------------------------------------------*/
template<typename T>
ARCCORE_HOST_DEVICE inline void
min_at_begin(SmallSpan<T> arr)
{
  if (arr.size() <= 1)
    return;
  
  Integer imin = 0;
  T vmin = arr[0];

  for (Integer i = 1; i < arr.size(); ++i) {
    if (arr[i] < vmin) {
      vmin = arr[i];
      imin = i;
    } 
  }
  if (imin!=0) {
    T vtmp; // on swap arr[0] <==> arr[imin]
    vtmp = arr[0];
    arr[0] = arr[imin];
    arr[imin] = vtmp;
  }
}

/*---------------------------------------------------------------------------*/
/*!
 * \brief Tri d'un tableau de petite taille (<10) à l'intérieur d'un kernel
 *        (on évite d'appeler std::sort sur GPU)
 *
 * @param[inout] arr : [in] tableau à trier [out] tableau trié
 *
 */
/*---------------------------------------------------------------------------*/
template<typename ArrayType>
ARCCORE_HOST_DEVICE inline void
sort_small(ArrayType& arr)
{
  SmallSpan<Real> arr0(arr.data(), arr.size());
  const Integer sz=arr.size();
  const Integer sz_m1=sz-1;
  for(Integer i=0; i<sz_m1 ; ++i) {
    min_at_begin(arr0.subSpan(i,sz-i));
  }
}

} // namespace Accenv

#endif
