// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
/*---------------------------------------------------------------------------*/
/* ICartesianMesh.h                                            (C) 2000-2012 */
/*                                                                           */
/* Interface d'un maillage cartésien.                                        */
/*---------------------------------------------------------------------------*/
#ifndef CARTESIAN_ICARTESIANMESH_H
#define CARTESIAN_ICARTESIANMESH_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/ArcaneTypes.h"
#include "arcane/cea/CeaGlobal.h"
#include "arcane/IMesh.h"
#include "arcane/ArcaneTypes.h"
#include "arcane/utils/ITraceMng.h"

//using IMesh = Arcane::IMesh;
//using eMeshDirection = Arcane::eMeshDirection;
//using ITraceMng = Arcane::ITraceMng;
using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class CellDirectionMng;
class FaceDirectionMng;
class NodeDirectionMng;
class CartesianConnectivity;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneCartesianMesh
 * \brief Interface d'un maillage cartésien.
 */
class ARCANE_CEA_EXPORT ICartesianMesh
{
 public:

  virtual ~ICartesianMesh() {} //<! Libère les ressources

  /*!
   * \brief Récupère ou créé la référence associée à \a mesh.
   *
   * Si aucun gestionnaire de matériau n'est associé à \a mesh, il
   * sera créé lors de l'appel à cette méthode si \a create vaut \a true.
   * Si \a create vaut \a false est qu'aucune gestionnaire n'est associé
   * au maillage, un pointeur nul est retourné.
   * L'instance retournée reste valide tant que le maillage \a mesh existe.
   */
  static ICartesianMesh* getReference(IMesh* mesh,bool create=true);

 public:

  virtual void build() =0;

 public:

  //! Maillage associé à ce maillage cartésien
  virtual IMesh* mesh() const =0;

  //! Gestionnaire de trace associé.
  virtual ITraceMng* traceMng() const =0;

  //! Liste des mailles dans la direction \a dir
  virtual CellDirectionMng cellDirection(eMeshDirection dir) =0;

  //! Liste des mailles dans la direction \a dir (0, 1 ou 2)
  virtual CellDirectionMng cellDirection(Integer idir) =0;

  //! Liste des faces dans la direction \a dir
  virtual FaceDirectionMng faceDirection(eMeshDirection dir) =0;

  //! Liste des faces dans la direction \a dir (0, 1 ou 2)
  virtual FaceDirectionMng faceDirection(Integer idir) =0;

  //! Liste des noeuds dans la direction \a dir
  virtual NodeDirectionMng nodeDirection(eMeshDirection dir) =0;

  //! Liste des noeuds dans la direction \a dir (0, 1 ou 2)
  virtual NodeDirectionMng nodeDirection(Integer idir) =0;

  /*!
   * \brief Calcule les infos pour les accès par direction.
   *
   * Actuellement, les restrictions suivantes existent:
   * - calcule uniquement les infos sur les entités mailles.
   * - suppose que la maille 0 est dans un coin (ne fonctionne que
   * pour le meshgenerator).
   * - cela ne marche que en sequentiel.
   */
  virtual void computeDirections() =0;

  //! Informations sur la connectivité
  virtual CartesianConnectivity connectivity() =0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

extern "C++" ARCANE_CEA_EXPORT ICartesianMesh*
arcaneCreateCartesianMesh(IMesh* mesh);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  

