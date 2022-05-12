#ifndef CARTESIAN_CARTESIAN_MESH_PROPERTIES_H
#define CARTESIAN_CARTESIAN_MESH_PROPERTIES_H

#include "cartesian/CartTypes.h"
#include "arcane/IMesh.h"
#include "arcane/IGhostLayerMng.h"
#include "arcane/Properties.h"
#include <algorithm>

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Recuperation des propriétés cartésiennes d'un maillage comme son nb de mailles Nx x Ny x Nz
 */
/*---------------------------------------------------------------------------*/
class CartesianMeshProperties {
 public:
  CartesianMeshProperties(IMesh* mesh)
  : m_dimension (mesh->dimension()) {
    // Pré-calcul du nb de mailles dans chaque direction
    LocalIdType nghost_layer = mesh->ghostLayerMng()->nbGhostLayer(); // Nb de couches de mailles fantomes

    // Le maillage sera cartésien ssi own_offset_dir != -1 pour toutes les directions
    m_is_pure_cartesian_mesh = true;

    String str_dir[3];
    str_dir[0] = String("X");
    str_dir[1] = String("Y");
    str_dir[2] = String("Z");

    Properties* mesh_properties = mesh->properties();
    for(Integer d(0) ; d < m_dimension ; ++d) {

      // Hors mailles fantomes, les indices globaux des mailles du sous-domaine
      // sont dans [own_offset_dir[d], own_offset_dir[d]+own_ncells_dir[d][
      LocalIdType own_ncells_dir = mesh_properties->getInt32WithDefault(String("OwnNbCell")+str_dir[d],-1);
      UniqueIdType own_offset_dir = mesh_properties->getInt64WithDefault(String("OwnCellOffset")+str_dir[d],-1);
      m_glob_ncells_dir[d] = mesh_properties->getInt64WithDefault(String("GlobalNbCell")+str_dir[d],-1);

      _computeDir(own_ncells_dir, own_offset_dir, m_glob_ncells_dir[d], nghost_layer, m_ncells_dir[d], m_is_pure_cartesian_mesh);
    }
  }

  //! Dimension du maillage
  Integer dimension() const { return m_dimension; }

  //! Retourne vrai si l'on travaille sur un maillage réellement cartésien, faux sinon
  bool isPureCartesianMesh() const {
    return m_is_pure_cartesian_mesh;
  }

  //! Retourne le nombre de mailles dans la direction dir (mailles fantomes comprises).
  // A appeler que si isPureCartesianMesh() est vrai, valeur quelconque sinon
  LocalIdType nbCellDir(Integer dir) {
    ARCANE_ASSERT(dir<m_dimension, ("Direction depasse dimension du maillage"));
    return m_ncells_dir[dir];
  }

  //! Retourne le nombre de mailles pour toutes les directions (mailles fantomes comprises), et 1 pour les directions >= dimension().
  // A appeler que si isPureCartesianMesh() est vrai
  const LocalIdType3& nbCell3() {
    return m_ncells_dir;
  }

  //! Retourne le nombre de mailles pour toutes les directions du domaine global, et 1 pour les directions >= dimension().
  // A appeler que si isPureCartesianMesh() est vrai
  const UniqueIdType3& globalNbCell3() {
    return m_glob_ncells_dir;
  }

 private:
  void _computeDir(LocalIdType own_ncells_dir, UniqueIdType own_offset_dir, UniqueIdType glob_ncells_dir, LocalIdType nghost_layer,
      LocalIdType& ncells_dir, bool& is_pure_cartesian_mesh) {
    if (own_offset_dir != -1) {
      // Nb de mailles du sous-domaine dans la direction d = nb de mailles interieures + des couches fantomes
      UniqueIdType u_nghost_layer = UniqueIdType(nghost_layer); // pour opérations entre UniqueIdType
      ncells_dir = own_ncells_dir;
      // Bord gauche du domaine global ssi own_offset_dir == 0
      // il faut ajouter des mailles fantomes a gauche sans faire déborder vers des indices négatifs
      ncells_dir += LocalIdType(std::min(u_nghost_layer, own_offset_dir));

      // Bord droit du domaine global ssi own_offset_dir+own_ncells_dir == glob_ncells_dir, 
      // il faut ajouter des mailles fantomes a droite sans faire déborder vers des indices sup. à glob_ncells_dir
      UniqueIdType own_ncells_after = glob_ncells_dir - (own_offset_dir+LocalIdType(own_ncells_dir));
      ncells_dir += LocalIdType(std::min(u_nghost_layer, own_ncells_after));
    } else {
      is_pure_cartesian_mesh = false;
    }
  }
 private:
  Integer m_dimension = 0; //! Dimension du maillage cartésien
  bool m_is_pure_cartesian_mesh = false; //! vrai s'il s'agit réellement d'un maillage cartésien, faux sinon
  LocalIdType3 m_ncells_dir = {1, 1, 1}; //! Nb de mailles cartésiennes dans chaque direction, mailles fantomes comprises
  UniqueIdType3 m_glob_ncells_dir = {1, 1, 1}; //! Nb de mailles cartésiennes dans chaque direction pour le domaine global
};
  
}

#endif

