#ifndef CARTESIAN_CART_INTERVAL3_T
#define CARTESIAN_CART_INTERVAL3_T

#include <array>

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Encapsulation d'un ensemble d'entiers décrit par un intervalle dans Z^3
 */
/*---------------------------------------------------------------------------*/
template<typename IdType>
class Interval3T {
 public:
  using IdType3 = IdType[3];
  using AIdType = std::array<IdType, 3>;

 public:
  Interval3T(const IdType3 beg, const IdType3 end)
  {
    if (beg[0]<end[0] && beg[1]<end[1] && beg[2]<end[2]) {
      // L'ensemble n'est pas vide et l'on récupère les bornes
      for(Integer d(0) ; d < 3 ; ++d) {
        m_beg[d] = beg[d];
        m_end[d] = end[d];
      }
    } else {
      // L'ensemble est vide et l'on affecte tout à 0
      for(Integer d(0) ; d < 3 ; ++d) {
        m_beg[d] = 0;
        m_end[d] = 0;
      }
    }
    m_size = (m_end[0]-m_beg[0])*(m_end[1]-m_beg[1])*(m_end[2]-m_beg[2]);
  }

  Interval3T(const Interval3T<IdType>& rhs) {
    for(Integer d(0) ; d < 3 ; ++d) {
      m_beg[d] = rhs.m_beg[d];
      m_end[d] = rhs.m_end[d];
    }
    m_size = rhs.m_size;
  }

  //! Retourne les bornes inférieures incluses
  AIdType lowerBounds() const {
    return m_beg;
  }

  //! Retourne les bornes supérieures exclues
  AIdType upperBounds() const {
    return m_end;
  }

  //! Nombre d'éléments par direction
  AIdType size3() const {
    return {m_end[0]-m_beg[0], m_end[1]-m_beg[1], m_end[2]-m_beg[2]};
  }

  //! Nombre d'éléments dans l'ensemble
  IdType size() const {
    return m_size;
  }

 protected:
  // Bornes de l'ensemble [m_beg[0], m_end[0][ x [m_beg[1], m_end[1][ x [m_beg[2], m_end[2][
  AIdType m_beg;
  AIdType m_end;
  IdType m_size; //! Nb d'éléments dans l'ensmble = (m_end[0]-m_beg[0])*(m_end[1]-m_beg[1])*(m_end[2]-m_beg[2])
};

}
#endif
