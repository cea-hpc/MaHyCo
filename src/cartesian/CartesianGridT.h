#ifndef CARTESIAN_CARTESIAN_GRID_H
#define CARTESIAN_CARTESIAN_GRID_H

#include "arcane/utils/ArcaneGlobal.h"
#include "cartesian/CartesianNumberingT.h"

namespace Cartesian {
  

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Encapsulation d'une grille cartesienne avec les mailles, noeuds, faces
 * d'une dimension au plus 3
 */
/*---------------------------------------------------------------------------*/
template<typename IdType>
class CartesianGridT
{
 public:
  //! Type pour les triplets cartésiens (i,j,k) et les triplets des dimensions (ni,nj,nk)
  using IdType3 = IdType[3];

  //! Type de la numérotation cartésienne associé à IdType
  using CartesianNumbering = CartesianNumberingT<IdType>;

  //! Type tableau numérotations cartésiennes sur les 3 dimensions
  using CartesianNumbering3 = CartesianNumbering[3];

 public:
  //! param[in] ncells_dir Nombre de mailles dans chaque direction
  //  param[in] dimension  Dimension de la grille cartésienne (au plus 3)
  CartesianGridT(const IdType3 &ncells_dir, Integer dimension)
  : m_dimension(dimension)
  {

    for(Integer d(0) ; d < m_dimension ; ++d) {
      m_ncells_dir[d] = ncells_dir[d];
      m_nnodes_dir[d] = m_ncells_dir[d]+1;
    }
    m_cart_num_cell.initNumbering(m_ncells_dir, dimension);

    for(Integer d(m_dimension) ; d < 3 ; ++d) {
      m_ncells_dir[d] = 1;
      m_nnodes_dir[d] = 1;
    }
    m_cart_num_node.initNumbering(m_nnodes_dir, dimension);

    // Pour les faces, on va numéroter les faces selon X, puis selon Y et selon Z
    // On distingue les faces par leurs orientations (dnorm)
    IdType nfaces_norm_first = 0; // Premier numero de face selon la normale dnorm
    for(Integer dnorm(0) ; dnorm < m_dimension ; ++dnorm) {

      // Les directions orthogonales à la normale
      Integer d1 = (dnorm+1) % 3;
      Integer d2 = (dnorm+2) % 3;

      // Dans la direction de la normale, on a m_ncells_dir[dnorm]+1 faces
      // et dans les directions m_ncells_dir[d] faces
      m_nfaces_norm_dir[dnorm][dnorm] = m_ncells_dir[dnorm]+1;
      m_nfaces_norm_dir[dnorm][d1] = m_ncells_dir[d1];
      m_nfaces_norm_dir[dnorm][d2] = m_ncells_dir[d2];

      m_cart_num_face[dnorm].initNumbering(m_nfaces_norm_dir[dnorm], dimension, nfaces_norm_first);
      nfaces_norm_first += m_cart_num_face[dnorm].nbItem();
    }

    for(Integer dnorm(m_dimension) ; dnorm < 3 ; ++dnorm) {
      m_nfaces_norm_dir[dnorm][0] = 1;
      m_nfaces_norm_dir[dnorm][1] = 1;
      m_nfaces_norm_dir[dnorm][2] = 1;
    }
  }

  //! Référence en lecture sur la numérotation cartésienne aux mailles
  const CartesianNumbering& cartNumCell() const {
    return m_cart_num_cell;
  }

  //! Référence en lecture sur la numérotation cartésienne aux noeuds
  const CartesianNumbering& cartNumNode() const {
    return m_cart_num_node;
  }

  //! Référence en lecture sur la numérotation cartésienne aux faces dans la direction \a dir
  const CartesianNumbering& cartNumFace(Integer dir) const {
    ARCANE_ASSERT(dir < m_dimension, ("La direction doit être strictement inférieure à la dimension"));
    return m_cart_num_face[dir];
  }

  //! Référence en lecture sur les 3 numérotations cartésiennes aux faces
  const CartesianNumbering3& cartNumFace3() const {
    return m_cart_num_face;
  }

  //! Pointeur sur la numérotation cartésienne aux mailles
  CartesianNumbering* cartNumCellPtr() {
    return &m_cart_num_cell;
  }

  //! Pointeur sur la numérotation cartésienne aux noeuds
  CartesianNumbering* cartNumNodePtr() {
    return &m_cart_num_node;
  }

  //! Pointeur sur la numérotation cartésienne aux faces dans la direction \a dir
  CartesianNumbering* cartNumFacePtr(Integer dir) {
    ARCANE_ASSERT(dir < m_dimension, ("La direction doit être strictement inférieure à la dimension"));
    return &(m_cart_num_face[dir]);
  }

  //! Pointeur sur les 3 numérotations cartésiennes aux faces
  CartesianNumbering3* cartNumFace3Ptr()  {
    return &m_cart_num_face;
  }


  //! Dimension du maillage cartésien
  Integer dimension() const {
    return m_dimension;
  }


 protected:

  IdType3 m_ncells_dir = {1, 1, 1};  // Nb de mailles par direction
  IdType3 m_nnodes_dir = {1, 1, 1};  // Nb de noeuds par direction
  IdType3 m_nfaces_norm_dir[3]; //! m_nfaces_norm_dir[dnorm] = dimension de la grille de faces normales à dnorm

  Integer m_dimension = 0;

  CartesianNumbering m_cart_num_cell;
  CartesianNumbering m_cart_num_node;
  CartesianNumbering3 m_cart_num_face;
};


}

#endif

