// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
/*---------------------------------------------------------------------------*/
/* CellDirectionMng.cc                                         (C) 2000-2014 */
/*                                                                           */
/* Infos sur les mailles d'une direction X Y ou Z d'un maillage structuré.   */
/*---------------------------------------------------------------------------*/
#ifndef CARTESIAN_CELLDIRECTIONMNG_H
#define CARTESIAN_CELLDIRECTIONMNG_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/ArcaneTypes.h"

#include "arcane/Item.h"
#include "arcane/ItemEnumerator.h"

#include "cartesian/CartesianGlobal.h"

#include "cartesian/CartTypes.h"
#include "cartesian/CartesianNumberingT.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class CellDirectionMng;
class ICartesianMesh;

/*!
 * \brief Position des noeuds de la maille par direction pour les maillages
 * cartésiens.
 */
enum eCellNodePosition
{
  CNP_NextLeft = 0,
  CNP_NextRight = 1,
  CNP_PreviousRight = 2,
  CNP_PreviousLeft = 3,

  CNP_TopNextLeft = 4,
  CNP_TopNextRight = 5,
  CNP_TopPreviousRight = 6,
  CNP_TopPreviousLeft = 7
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneCartesianMesh
 * \brief Maille avant et après une maille suivant une direction.
 *
 * Les instances de cette classe sont temporaires et construites via
 * CellDirectionMng::cell().
 */
class ARCANE_CEA_EXPORT DirCell
{
  //! Type tableau sur pointeurs d'ItemInternal (implémentation d'un Item)
  using ItemInternalPtr = Item::ItemInternalPtr;

 public:
  DirCell(LocalIdType cell_id, LocalIdType idx_dir, LocalIdType ncellsm1_dir, LocalIdType delta_dir, 
    const ItemInternalPtr* cell_internals, ItemInternalPtr null_cell_internal)
  : m_cell_id (cell_id), 
  m_idx_dir (idx_dir),
  m_ncellsm1_dir (ncellsm1_dir),
  m_delta_dir (delta_dir),
  m_cell_internals (cell_internals),
  m_null_cell_internal (null_cell_internal)
  {
  }

 public:
  //! Maille avant
  Cell previous() const { return Cell(m_idx_dir>0              ? m_cell_internals[m_cell_id-m_delta_dir] : m_null_cell_internal); }
  //! Maille après
  Cell next()     const { return Cell(m_idx_dir<m_ncellsm1_dir ? m_cell_internals[m_cell_id+m_delta_dir] : m_null_cell_internal); }
 private:
  LocalIdType m_cell_id; // le local id de la maille
  LocalIdType m_idx_dir; // indice cartesien de la maille dans la direction dir
  LocalIdType m_ncellsm1_dir; // Nb de mailles -1 dans la direction dir
  LocalIdType m_delta_dir; // +-delta a appliquer sur m_cell_id pour passer a la maille suivante/precedente
  const ItemInternalPtr* m_cell_internals;
  ItemInternalPtr m_null_cell_internal;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneCartesianMesh
 * \brief Maille avec info directionnelle des noeuds.
 *
 * Les instances de cette classe sont temporaires et construites via
 * CellDirectionMng::cellNode().
 */
class ARCANE_CEA_EXPORT DirCellNode
{
 public:
  DirCellNode(Cell c,Int32* nodes_indirection)
  : m_cell(c), m_nodes_indirection(nodes_indirection){}

 public:

  //! Maille associée
  Cell cell() const { return m_cell; }

  //! Noeud devant à gauche dans la direction
  Node nextLeft() const { return m_cell.node(m_nodes_indirection[CNP_NextLeft]); }
  //! Noeud devant à droite dans la direction
  Node nextRight() const { return m_cell.node(m_nodes_indirection[CNP_NextRight]); }
  //! Noeud derrière à droite dans la direction
  Node previousRight() const { return m_cell.node(m_nodes_indirection[CNP_PreviousRight]); }
  //! Noeud derrière à gauche dans la direction
  Node previousLeft() const { return m_cell.node(m_nodes_indirection[CNP_PreviousLeft]); }

  //! Noeud devant à gauche dans la direction
  Node topNextLeft() const { return m_cell.node(m_nodes_indirection[CNP_TopNextLeft]); }
  //! Noeud devant à droite dans la direction
  Node topNextRight() const { return m_cell.node(m_nodes_indirection[CNP_TopNextRight]); }
  //! Noeud derrière à droite dans la direction
  Node topPreviousRight() const { return m_cell.node(m_nodes_indirection[CNP_TopPreviousRight]); }
  //! Noeud derrière à gauche dans la direction
  Node topPreviousLeft() const { return m_cell.node(m_nodes_indirection[CNP_TopPreviousLeft]); }

 private:
  
  Cell m_cell;
  Int32* m_nodes_indirection;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneCartesianMesh
 * \brief Maille avec info directionnelle des faces.
 *
 * Les instances de cette classe sont temporaires et construites via
 * CellDirectionMng::cellFace().
 */
class ARCANE_CEA_EXPORT DirCellFace
{
 public:
  DirCellFace(Cell c,Int32 next_face_index,Int32 previous_face_index)
  : m_cell(c), m_next_face_index(next_face_index), m_previous_face_index(previous_face_index)
  {
  }

 public:

  //! Maille associée
  Cell cell() const { return m_cell; }

  //! Face connectée à la maille d'après la maille courante dans la direction
  Face next() const { return m_cell.face(m_next_face_index); }

  //! Face connectée à la maille d'avant la maille courante dans la direction
  Face previous() const { return m_cell.face(m_previous_face_index); }

  //! Indice locale dans la maille de la face next()
  Int32 nextLocalIndex() const { return m_next_face_index; }

  //! Indice locale dans la maille de la face previous()
  Int32 previousLocalIndex() const { return m_previous_face_index; }

 private:
  
  Cell m_cell;
  Int32 m_next_face_index;
  Int32 m_previous_face_index;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneCartesianMesh
 * \brief Infos sur les mailles d'une direction spécifique X,Y ou Z d'un maillage structuré.
 *
 * Cette classe contient les informations pour retourner la liste
 * des mailles dans une direction donnée et pour ces mailles
 * connaitre la maille avant et après dans cette direction.
 *
 * Les instances de cette classe sont gérées par un ICartesianMesh et sont
 * temporaires. Il ne faut pas les conserver d'une itération à l'autre
 * car elles sont invalidées si le maillage change.
 *
 * Cette classe à une sémantique par référence.
 *
 * Par exemple, pour itérer sur les mailles de la direction X:
 
 \code
 * ICartesianMesh* cartesian_mesh = ...;
 * CellDirectionMng cdm(cartesian_mesh->cellDirection(MD_DirX));
 * ENUMERATE_CELL(icell,cdm.allCells()){
 *   DirCell dir_cell(cdm[icell]);
 *   Cell next = dir_cell.next();
 *   Cell prev = dir_cell.previous();
 * }
 \endcode
 *
 */
class ARCANE_CEA_EXPORT CellDirectionMng
{
  friend class CartesianMesh;
  class Impl;
  static const int MAX_NB_NODE = 8;

 private:

  //! Type pour la numérotation cartésienne sur des identifiants locaux
  using CartesianNumbering = CartesianNumberingT<LocalIdType>;

 public:
  
  //! Créé une instance vide. L'instance n'est pas valide tant que init() n'a pas été appelé.
  CellDirectionMng();
  CellDirectionMng(const CellDirectionMng& rhs);
  ~CellDirectionMng();

  //! Maille direction correspondant à la maille \a c.
  DirCell cell(Cell c)
  {
    return _cell(c.localId());
  }

  //! Maille avec infos directionnelles aux noeuds correspondant à la maille \a c.
  DirCellNode cellNode(Cell c)
  {
    return DirCellNode(c,m_nodes_indirection);
  }

  //! Maille avec infos directionnelles aux faces correspondant à la maille \a c.
  DirCellFace cellFace(Cell c)
  {
    return DirCellFace(c,m_next_face_index,m_previous_face_index);
  }

  //! Groupe de toutes les mailles dans la direction.
  CellGroup allCells() const;

  /*!
   * \brief Groupe de toutes les mailles internes dans la direction.
   *
   * Une maille est considérée comme interne si sa maille
   * avant ou après n'est pas nulle.
   */

  CellGroup innerCells() const;
  /*!
   * \brief Groupe de toutes les mailles externes dans la direction.
   *
   * Une maille est considérée comme externe si sa maille
   * avant ou après est nulle.
   */
  CellGroup outerCells() const;

  //! Maille direction correspondant à la maille \a c.
  DirCell operator[](Cell c)
  {
    return _cell(c.localId());
  }

  //! Maille direction correspondant à l'itérateur de la maille \a icell.
  DirCell operator[](CellEnumerator icell)
  {
    return _cell(icell.itemLocalId());
  }

  /*!
   * \brief Nombre global de mailles dans cette direction.
   *
   * \note La valeur retournée n'est valide que si le
   * maillage a été créé avec un générateur spécifique, tel
   * le SodMeshGenerator ou le CartesianMeshGenerator. Si ce n'est
   * pas le cas, la valeur retournée vaut (-1)
   */
  Int64 globalNbCell() const { return m_global_nb_cell; }

  /*!
   * \brief Nombre de mailles propres dans cette direction.
   *
   * \note La valeur retournée n'est valide que si le
   * maillage a été créé avec un générateur spécifique, tel
   * le SodMeshGenerator ou le CartesianMeshGenerator. Si ce n'est
   * pas le cas, la valeur retournée vaut (-1)
   */
  Int32 ownNbCell() const { return m_own_nb_cell; }

  /*!
   * \brief Offset dans cette direction du sous-domaine.
   *
   * En supposant que le maillage cartésien global est découpé en
   * plusieurs sous-domaines rectangulaires qui forment une grille,
   * cette méthode retourne la position dans cette grille de ce sous-domaine
   * pour cette direction.
   *
   * \warning L'utilisation de cette méthode suppose que chaque
   * sous-domaine est parallélépipédique (en 3D) ou rectangulaire (en 2D)
   * ce qui n'est pas forcément le cas, notamment avec des mécanismes
   * d'équilibrage de charge par migration de maille.
   *
   * \note La valeur retournée n'est valide que si le
   * maillage a été créé avec un générateur spécifique, tel que
   * le CartesianMeshGenerator. Si ce n'est pas le cas,
   * la valeur retournée vaut (-1)
   */
  Int32 subDomainOffset() const { return m_sub_domain_offset; }

  /*!
   * \brief Offset dans cette direction de la première maille propre de ce sous-domaine.
   *
   * En supposant que le maillage cartésien global est découpé en
   * plusieurs sous-domaines rectangulaires qui forment une grille,
   * cette méthode retourne la position dans cette grille de la première
   * maille propre de ce sous-domaine pour cette direction.
   *
   * \warning L'utilisation de cette méthode suppose que chaque
   * sous-domaine est parallélépipédique (en 3D) ou rectangulaire (en 2D)
   * ce qui n'est pas forcément le cas, notamment avec des mécanismes
   * d'équilibrage de charge par migration de maille.
   *
   * \note La valeur retournée n'est valide que si le
   * maillage a été créé avec un générateur spécifique, tel que
   * le CartesianMeshGenerator. Si ce n'est pas le cas,
   * la valeur retournée vaut (-1)
   */
  Int64 ownCellOffset() const { return m_own_cell_offset; }

 private:
  
  //! Maille direction correspondant à la maille de numéro local \a local_id
  DirCell _cell(Int32 local_id)
  {
    LocalIdType idx_dir; // Trouver un moyen d'eviter les ifs sans utiliser des methodes virtuelles
    if (m_dir==0)
      idx_dir=m_p_cart_numb_cell->idxDir0(local_id); // couteux
    else if (m_dir==1)
      idx_dir=m_p_cart_numb_cell->idxDir1(local_id); // couteux
    else 
      idx_dir=m_p_cart_numb_cell->idxDir2(local_id); // couteux

    return DirCell(local_id, /*idx3[m_dir]*/idx_dir, m_ncellsm1_dir, m_delta_dir, m_cell_internals.data(), m_null_cell_internal);
  }

  void setNodesIndirection(Int32ConstArrayView nodes_indirection);

 public:

  /*!
   * \internal
   * \brief Usage interne à Arcane. Calcul les entités internes et externes.
   * Suppose que init() a été appelé.
   */
  void computeInnerAndOuterItems(const ItemGroup& items);

  /*!
   * \internal
   * Initialise l'instance.
   */
  void init(ICartesianMesh* cm,eMeshDirection dir,CartesianNumbering* p_cart_numb_cell);

  /*!
   * \internal
   * Détruit les ressources associées à l'instance.
   */
  void destroy();

  //! Valeur de la direction
  eMeshDirection direction() const
  {
    return m_direction;
  }

  /*!
   * \internal
   * Positionne les indices locaux de la face vers la maille d'après et d'avant.
   */
  void setLocalFaceIndex(Int32 next_index,Int32 previous_index)
  {
    m_next_face_index = next_index;
    m_previous_face_index = previous_index;
  }

 private:

  eMeshDirection m_direction;
  Integer m_dir;  //! m_direction casté en Integer
  Impl* m_p;
  Int32 m_nodes_indirection[MAX_NB_NODE];
  Int32 m_next_face_index;
  Int32 m_previous_face_index;
  Int32 m_sub_domain_offset;
  Int32 m_own_nb_cell;
  Int64 m_global_nb_cell;
  Int64 m_own_cell_offset;

  CartesianNumbering* m_p_cart_numb_cell=nullptr; //! Pointeur sur la numérotation cartésienne aux mailles
  ItemInternalList m_cell_internals;  //! Tableau dimensionne au nb total de Cell, chaque case pointe vers un ItemInternal
  LocalIdType m_ncellsm1_dir=0;
  LocalIdType m_delta_dir=0;

  Cell m_null_cell;  //! Une maille nulle (inexistante)
  ItemInternal* m_null_cell_internal=nullptr;  //! m_null_cell.internal()
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  

