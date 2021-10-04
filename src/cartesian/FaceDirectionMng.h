// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
/*---------------------------------------------------------------------------*/
/* FaceDirectionMng.cc                                         (C) 2000-2012 */
/*                                                                           */
/* Infos sur les faces d'une direction X Y ou Z d'un maillage structuré.     */
/*---------------------------------------------------------------------------*/
#ifndef CARTESIAN_FACEDIRECTIONMNG_H
#define CARTESIAN_FACEDIRECTIONMNG_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/ArcaneTypes.h"
#include "arcane/cea/CeaGlobal.h"

#include "arcane/Item.h"
#include "arcane/ItemEnumerator.h"

#include "cartesian/CartTypes.h"
#include "cartesian/CartesianGridT.h"
#include "cartesian/NumberingConverterT.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class FaceDirectionMng;
class ICartesianMesh;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneCartesianMesh
 * \brief Infos sur maille avant et après une face suivant une direction.
 *
 * Les instances de cette classe sont temporaires et construites via
 * FaceDirectionMng::face().
 */
class ARCANE_CEA_EXPORT DirFace
{
  //! Type tableau sur pointeurs d'ItemInternal (implémentation d'un Item)
  using ItemInternalPtr = Item::ItemInternalPtr;

 public:
  DirFace(LocalIdType face_idx_dir, LocalIdType nfacesm1_dir, LocalIdType next_cell_id, LocalIdType delta_prev_cell, 
    const ItemInternalPtr* cell_internals, ItemInternalPtr null_cell_internal) 
  : m_face_idx_dir (face_idx_dir),
  m_nfacesm1_dir (nfacesm1_dir),
  m_next_cell_id (next_cell_id),
  m_delta_prev_cell (delta_prev_cell),
  m_cell_internals (cell_internals),
  m_null_cell_internal (null_cell_internal)
  {
  }
 public:
  //! Maille avant
  Cell previousCell() const { 
    return Cell(m_face_idx_dir>0 ? m_cell_internals[m_next_cell_id-m_delta_prev_cell] : m_null_cell_internal); 
  }
  //! Maille après
  Cell nextCell() const { 
    return Cell(m_face_idx_dir<m_nfacesm1_dir ? m_cell_internals[m_next_cell_id] : m_null_cell_internal); 
  }
 private:
  LocalIdType m_face_idx_dir; //! m_face_ijk[m_dir] pré-calculé
  LocalIdType m_nfacesm1_dir; //! Nb de faces -1 dans la direction m_dir
  LocalIdType m_next_cell_id; //! L'identifiant de la maille qui suit la face
  LocalIdType m_delta_prev_cell; //! Delta à retirer pour passer de la maille suivante à précédente
  const ItemInternalPtr* m_cell_internals;
  ItemInternalPtr m_null_cell_internal;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneCartesianMesh
 * \brief Infos sur les face d'une direction spécifique X,Y ou Z
 * d'un maillage structuré.
 */
class ARCANE_CEA_EXPORT FaceDirectionMng
{
  friend class CartesianMesh;
  class Impl;

 private:

  //! Type pour la numérotation cartésienne sur des identifiants locaux
  using CartesianNumbering = CartesianNumberingT<LocalIdType>;
  //! Type pour une grille cartésienne sur des identifiants locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;
  //! Type pour convertir numérotation Face => Cell
  using NumbConvFaceCell = NumberingConverterT<Face,Cell>;

 public:
  
  //! Créé une instance vide. L'instance n'est pas valide tant que init() n'a pas été appelé.
  FaceDirectionMng();
  FaceDirectionMng(const FaceDirectionMng& rhs);
  ~FaceDirectionMng();

  //! Face direction correspondant à la face \a f.
  DirFace face(Face f)
  {
    return _face(f.localId());
  }

  //! Groupe de toutes les faces dans la direction.
  FaceGroup allFaces() const;

  /*!
   * \brief Groupe de toutes les faces internes dans la direction.
   *
   * Une face est considérée comme interne si sa maille
   * devant et derrière n'est pas nulle.
   */

  FaceGroup innerFaces() const;

  /*!
   * \brief Groupe de toutes les faces externes dans la direction.
   *
   * Une face est considérée comme externe si sa maille
   * devant ou derrière est nulle.
   */
  FaceGroup outerFaces() const;

  //! Face direction correspondant à la face \a f.
  DirFace operator[](Face f)
  {
    return _face(f.localId());
  }

  //! Face direction correspondant à l'itérateur de la face \a iface
  DirFace operator[](FaceEnumerator iface)
  {
    return _face(iface.itemLocalId());
  }

 private:
  
  //! Face direction correspondant à la face de numéro local \a local_id
  DirFace _face(Int32 local_id)
  {
    LocalIdType3 face_idx3; 
    m_p_cart_numb_face_dir->ijk(local_id, face_idx3); // couteux

    // Passage numérotation Face => Cell
    LocalIdType next_cell_id = local_id+m_p_numb_conv->computeDelta(face_idx3[1], face_idx3[2]);
    return DirFace(face_idx3[m_dir], m_nfacesm1_dir, next_cell_id, m_delta_prev_cell, m_cell_internals.data(), m_null_cell_internal);
  }

 public:

  /*!
   * \internal
   * \brief Usage interne à Arcane. Calcul les entités internes et externes.
   * Suppose que init() a été appelé.
   */
  void computeInnerAndOuterItems(const ItemGroup& items);

  /*!
   * \internal
   * \brief Usage interne à Arcane. Calcul les entités internes et externes, implémentation cartésienne.
   * Suppose que init() a été appelé.
   */
  void computeInnerAndOuterItemsFromFamily(IItemFamily* face_family);

  /*!
   * \internal
   * Initialise l'instance.
   */
  void init(ICartesianMesh* cm,eMeshDirection dir,
    CartesianGrid* p_cart_grid, NumbConvFaceCell* p_numb_conv);

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

 private:

  eMeshDirection m_direction;
  Integer m_dir;  //! m_direction casté en Integer
  Impl* m_p;

  CartesianNumbering* m_p_cart_numb_face_dir=nullptr; //! Pointeur sur la numérotation cartésienne aux faces pour m_dir
  NumbConvFaceCell* m_p_numb_conv=nullptr;  //! Pointeur sur le convertisseur ids Face => Cell
  ItemInternalList m_cell_internals;  //! Tableau dimensionne au nb total de Cell, chaque case pointe vers un ItemInternal
  LocalIdType m_nfacesm1_dir=0;
  LocalIdType m_delta_prev_cell=0;

  Cell m_null_cell;  //! Une maille nulle (inexistante)
  ItemInternal* m_null_cell_internal=nullptr;  //! m_null_cell.internal()
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  

