#ifndef CARTESIAN_NUMBERING_CONVERTER_T_H
#define CARTESIAN_NUMBERING_CONVERTER_T_H

#include "cartesian/CartesianNumberingT.h"
#include "cartesian/CartesianGridT.h"
#include "arcane/Item.h"
#include "arcane/utils/ArcaneGlobal.h"


namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Permet de passer d'une numérotation cartésienne d'un ItemType0 vers un ItemType1
 */
/*---------------------------------------------------------------------------*/
// Type général qui ne doit janais être utilisé et qui doit être surchargé
template<typename ItemType0, typename ItemType1>
class NumberingConverterT {
 public:
  //! Type de grille cartésienne sur des ids locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;

 public:
  NumberingConverterT(Integer dir, const CartesianGrid &cart_grid) 
  : m_dir(dir), m_cart_grid (cart_grid) {
  }

  void initDelta() {
    m_delta = 0;
  }

  //! Retourne le delta à ajouter au local id de la face courante pour obtenir le local id de la maille next()
  LocalIdType delta() const {
    ARCANE_ASSERT(false, ("Pas d'implémentation générique de NumberingConverterT"));
    return m_delta;
  }

  //! Calcul du delta de conversion à partir de (j,k)
  LocalIdType computeDelta(LocalIdType , LocalIdType) const {
    return 0;
  }

  //! Calcul et mise à jour du delta de conversion à partir de (j,k)
  void updateDelta(LocalIdType , LocalIdType) {
    return;
  }

 private:
  Integer m_dir; //! Direction dans laquelle on veut passer de Face à Cell
  const CartesianGrid &m_cart_grid;

  LocalIdType m_delta = 0;
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Spécialisation conversion Face => Cell
 * Le but est de calculer delta_d tq : cell_id(i,j,k) = face_id(d, i,j,k) + delta_d
 * Pour dir = x : delta_x = FC - FF_x - (j+k*Nj)
 * Pour dir = y : delta_y = FC - FF_y -    k*Ni
 * Pour dir = z : delta_z = FC - FF_z
 * avec FC   = first ID Cell
 *      FF_d = first ID Face dans la direction d
 *      Ni (resp. Nj) = Nb de mailles selon i (resp. j)
 */
/*---------------------------------------------------------------------------*/
template<>
class NumberingConverterT<Face, Cell> {
 public:
  //! Type de grille cartésienne sur des ids locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;

 public:
  NumberingConverterT(Integer dir, const CartesianGrid &cart_grid) 
  : m_dir(dir) {
    const auto& cart_numb_cell = cart_grid.cartNumCell();
    const auto& cart_numb_face = cart_grid.cartNumFace(m_dir);

    m_first_delta = cart_numb_cell.firstId() - cart_numb_face.firstId();  // FC-FF_d

    if (m_dir == MD_DirX) {
      m_nitems = cart_numb_cell.nbItemDir(MD_DirY); // Nj
    } else if (m_dir == MD_DirY) {
      m_nitems = cart_numb_cell.nbItemDir(MD_DirX); // Ni
    } else {
      m_nitems = 0; // ne servira pas
    }
  }

  ARCCORE_HOST_DEVICE NumberingConverterT(const NumberingConverterT<Face,Cell>& rhs)
  : m_dir (rhs.m_dir),
  m_first_delta (rhs.m_first_delta),
  m_delta (rhs.m_delta),
  m_nitems (rhs.m_nitems)
  {
  }

  void initDelta() {
    m_delta = m_first_delta;
  }

  //! Retourne le delta à ajouter au local id de la face courante pour obtenir le local id de la maille next()
  LocalIdType delta() const {
    return m_delta;
  }

  //! Calcul du delta de conversion à partir de (j,k)
  ARCCORE_HOST_DEVICE LocalIdType computeDelta(LocalIdType j, LocalIdType k) const {
    LocalIdType delta;
    if (m_dir == MD_DirX) {
      delta = m_first_delta - (j+k*m_nitems);  // Nj
    } else if (m_dir == MD_DirY) {
      delta = m_first_delta - k*m_nitems;  // Ni
    } else {
      delta = m_first_delta;
    }
    return delta;
  }

  //! Calcul et mise à jour du delta de conversion à partir de (j,k)
  void updateDelta(LocalIdType j, LocalIdType k) {
    m_delta=computeDelta(j,k);
  }

 private:
  Integer m_dir; //! Direction dans laquelle on veut passer de Face à Cell

  LocalIdType m_first_delta;
  LocalIdType m_delta = 0;
  LocalIdType m_nitems = 0; 
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Spécialisation conversion Cell => Face
 * Le but est de calculer delta_d tq : face_id(d, i,j,k) = cell_id(i,j,k) + delta_d
 * Pour dir = x : delta_x = FF_x - FC + (j+k*Nj)
 * Pour dir = y : delta_y = FF_y - FC +    k*Ni
 * Pour dir = z : delta_z = FF_z - FC
 * avec FC   = first ID Cell
 *      FF_d = first ID Face dans la direction d
 *      Ni (resp. Nj) = Nb de mailles selon i (resp. j)
 */
/*---------------------------------------------------------------------------*/
template<>
class NumberingConverterT<Cell, Face> {
 public:
  //! Type de grille cartésienne sur des ids locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;

 public:
  NumberingConverterT(Integer dir, const CartesianGrid &cart_grid) 
  : m_dir(dir) {
    const auto& cart_numb_cell = cart_grid.cartNumCell();
    const auto& cart_numb_face = cart_grid.cartNumFace(m_dir);

    m_first_delta = cart_numb_face.firstId() - cart_numb_cell.firstId();  // FF_d - FC

    if (m_dir == MD_DirX) {
      m_nitems = cart_numb_cell.nbItemDir(MD_DirY); // Nj
    } else if (m_dir == MD_DirY) {
      m_nitems = cart_numb_cell.nbItemDir(MD_DirX); // Ni
    } else {
      m_nitems = 0; // ne servira pas
    }
  }

  ARCCORE_HOST_DEVICE NumberingConverterT(const NumberingConverterT<Cell,Face>& rhs)
  : m_dir (rhs.m_dir),
  m_first_delta (rhs.m_first_delta),
  m_delta (rhs.m_delta),
  m_nitems (rhs.m_nitems)
  {
  }

  void initDelta() {
    m_delta = m_first_delta;
  }

  //! Retourne le delta à ajouter au local id de la face courante pour obtenir le local id de la maille next()
  LocalIdType delta() const {
    return m_delta;
  }

  //! Calcul du delta de conversion à partir de (j,k)
  ARCCORE_HOST_DEVICE LocalIdType computeDelta(LocalIdType j, LocalIdType k) const {
    LocalIdType delta;
    if (m_dir == MD_DirX) {
      delta = m_first_delta + (j+k*m_nitems);  // Nj
    } else if (m_dir == MD_DirY) {
      delta = m_first_delta + k*m_nitems;  // Ni
    } else {
      delta = m_first_delta;
    }
    return delta;
  }

  //! Calcul et mise à jour du delta de conversion à partir de (j,k)
  void updateDelta(LocalIdType j, LocalIdType k) {
    m_delta=computeDelta(j,k);
  }

 private:
  Integer m_dir; //! Direction dans laquelle on veut passer de Face à Cell

  LocalIdType m_first_delta;
  LocalIdType m_delta = 0;
  LocalIdType m_nitems = 0; 
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Spécialisation conversion Cell => Node
 * Le but est de calculer delta tq : node_id(i,j,k) = cell_id(i,j,k) + delta
 * delta = FN-FC + j+k*(Ni+Nj+1)
 * avec FC = first ID Cell
 *      FN = first ID Node
 *      Ni (resp. Nj) = Nb de mailles selon i (resp. j)
 */
/*---------------------------------------------------------------------------*/
template<>
class NumberingConverterT<Cell, Node> {
 public:
  //! Type de grille cartésienne sur des ids locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;

 public:
  NumberingConverterT(Integer dir, const CartesianGrid &cart_grid) 
  : m_dir(dir) {
    const auto& cart_numb_cell = cart_grid.cartNumCell();
    const auto& cart_numb_node = cart_grid.cartNumNode();

    m_first_delta = cart_numb_node.firstId() - cart_numb_cell.firstId(); // FN-FC
    m_nitems = cart_numb_cell.nbItemDir(MD_DirX)+cart_numb_cell.nbItemDir(MD_DirY)+1;  // Ni+Nj+1
  }

  void initDelta() {
    m_delta = m_first_delta;
  }

  //! Retourne le delta à ajouter au local id de la maille courante (i,j,k) pour obtenir le local id du premier noeud (i,j,k)
  LocalIdType delta() const {
    return m_delta;
  }

  //! Calcul du delta de conversion à partir de (j,k)
  LocalIdType computeDelta(LocalIdType j, LocalIdType k) const {
    return m_first_delta + j+k*m_nitems;
  }

  //! Calcul et mise à jour du delta de conversion à partir de (j,k)
  void updateDelta(LocalIdType j, LocalIdType k) {
    m_delta=computeDelta(j,k);
  }

 private:
  Integer m_dir; //! Direction dans laquelle on veut passer de Cell à Node

  LocalIdType m_first_delta;
  LocalIdType m_delta = 0;
  LocalIdType m_nitems = 0; 
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Spécialisation conversion Node => Cell
 * Le but est de calculer delta tq : cell_id(i,j,k) = node_id(i,j,k) + delta
 * delta = FC-FN -j-k*(Ni+Nj+1)
 * avec FC = first ID Cell
 *      FN = first ID Node
 *      Ni (resp. Nj) = Nb de mailles selon i (resp. j)
 */
/*---------------------------------------------------------------------------*/
template<>
class NumberingConverterT<Node, Cell> {
 public:
  //! Type de grille cartésienne sur des ids locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;

 public:
  NumberingConverterT(Integer dir, const CartesianGrid &cart_grid) 
  : m_dir(dir) {
    const auto& cart_numb_cell = cart_grid.cartNumCell();
    const auto& cart_numb_node = cart_grid.cartNumNode();

    m_first_delta = cart_numb_cell.firstId() - cart_numb_node.firstId(); // FC-FN
    m_nitems = cart_numb_cell.nbItemDir(MD_DirX)+cart_numb_cell.nbItemDir(MD_DirY)+1;  // Ni+Nj+1
  }

  void initDelta() {
    m_delta = m_first_delta;
  }

  //! Retourne le delta à ajouter au local id du noeud courant (i,j,k) pour obtenir le local id de la maille (i,j,k)
  LocalIdType delta() const {
    return m_delta;
  }

  //! Calcul du delta de conversion à partir de (j,k)
  LocalIdType computeDelta(LocalIdType j, LocalIdType k) const {
    return m_first_delta -j-k*m_nitems;
  }

  //! Calcul et mise à jour du delta de conversion à partir de (j,k)
  void updateDelta(LocalIdType j, LocalIdType k) {
    m_delta=computeDelta(j,k);
  }

 private:
  Integer m_dir; //! Direction dans laquelle on veut passer de Node à Cell

  LocalIdType m_first_delta;
  LocalIdType m_delta = 0;
  LocalIdType m_nitems = 0; 
};

}

#endif

