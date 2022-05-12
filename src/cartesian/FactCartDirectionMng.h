#ifndef CARTESIAN_FACT_CART_DIRECTION_MNG_H
#define CARTESIAN_FACT_CART_DIRECTION_MNG_H

#include "arcane/IMesh.h"

#include "cartesian/CartTypes.h"
#include "cartesian/CartItemGroup.h"
#include "cartesian/CartesianGridT.h"
#include "cartesian/CartCellDirectionMng.h"
#include "cartesian/CartFaceDirectionMng.h"
#include "cartesian/CartNodeDirectionMng.h"
#include "cartesian/CartesianMeshProperties.h"

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Fabrique de Cart{Cell|Face}DirectionMng
 */
/*---------------------------------------------------------------------------*/
class FactCartDirectionMng {
 public:
  //! Type d'une grille cartésienne en numérotation locale
  using CartesianGrid = CartesianGridT<LocalIdType>;

 public:
  FactCartDirectionMng(IMesh *mesh) 
  : m_mesh (mesh),
  m_cart_mesh_prop (m_mesh),
  m_dimension (m_cart_mesh_prop.dimension()) {

    m_is_pure_cartesian_mesh = m_cart_mesh_prop.isPureCartesianMesh();

    if (m_is_pure_cartesian_mesh) {
      for(Integer d(0) ; d < m_dimension ; ++d) {
        m_ncells_dir[d] = m_cart_mesh_prop.nbCellDir(d);
      }
      m_cartesian_grid = new CartesianGrid(m_ncells_dir, m_dimension);
    } else {
      m_cartesian_grid = nullptr;
    }
  }

  virtual ~FactCartDirectionMng() {
    delete m_cartesian_grid;
  }

  //! Retourne vrai si l'on travaille sur un maillage réellement cartésien, faux sinon
  bool isPureCartesianMesh() const {
    return m_is_pure_cartesian_mesh;
  }

  //! Fabrique une instance de CartCellDirectionMng pour la direction dir
  CartCellDirectionMng cellDirection(Integer dir) const {
    ARCANE_ASSERT(m_cartesian_grid != nullptr, ("Le maillage doit être cartésien"));
    return CartCellDirectionMng(m_mesh->itemsInternal(IK_Cell).data(), dir, *m_cartesian_grid);
  }

  //! Fabrique une instance de CartFaceDirectionMng pour la direction dir
  CartFaceDirectionMng faceDirection(Integer dir) const {
    ARCANE_ASSERT(m_cartesian_grid != nullptr, ("Le maillage doit être cartésien"));
    return CartFaceDirectionMng(m_mesh->itemsInternal(IK_Face).data(), dir, *m_cartesian_grid);
  }

  //! Fabrique une instance de CartNodeDirectionMng pour la direction dir
  CartNodeDirectionMng nodeDirection(Integer dir) const {
    ARCANE_ASSERT(m_cartesian_grid != nullptr, ("Le maillage doit être cartésien"));
    return CartNodeDirectionMng(m_mesh->itemsInternal(IK_Node).data(), dir, *m_cartesian_grid);
  }

  //! Retourne l'instance de CartesianGrid si isPureCartesianMesh() == true, nullptr sinon
  CartesianGrid* cartesianGrid() {
    return m_cartesian_grid;
  }

 private:
  IMesh *m_mesh = nullptr; //! Le maillage 
  CartesianMeshProperties m_cart_mesh_prop;  //! Propriétés cartésiennes du maillage (si ce dernier est cartésien)
  Integer m_dimension = 0; //! Dimension du maillage m_mesh

  LocalIdType3 m_ncells_dir = {1, 1, 1}; //! Nb de mailles cartésiennes dans chaque direction, mailles fantomes comprises

  bool m_is_pure_cartesian_mesh = false; //! vrai si m_mesh est réellement cartésien, faux sinon
  CartesianGrid* m_cartesian_grid = nullptr; //! Grille cartésienne de mailles, faces, noeuds si m_is_pure_cartesian_mesh == true
};

}

#endif

