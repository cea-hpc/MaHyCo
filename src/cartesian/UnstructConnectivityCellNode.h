#ifndef CARTESIAN_UNSTRUCT_CONNECTIVITY_CELL_NODE_H
#define CARTESIAN_UNSTRUCT_CONNECTIVITY_CELL_NODE_H

#include "arcane/IMesh.h"
#include "arcane/ItemTypes.h"

namespace Cartesian {
  

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Connectivité maille => noeuds (hypothèse certésienne : toutes les mailles ont 2**dim noeuds)
 */
/*---------------------------------------------------------------------------*/
class UnstructConnectivityCellNode {
 public:
  UnstructConnectivityCellNode(IMesh* mesh) {

    Integer dim = mesh->dimension();
    m_nb_node = 1 << dim; //! Nombre de noeuds d'une maille quelconque pour la dimension de la grille cartésienne
  }

  // Nombre de noeuds d'une maille quelconque pour la dimension de la grille cartésienne
  inline Integer nbNode() const {
    return m_nb_node;
  }

  //! Retourne l'objet permettant de récupérer les noeuds de la maille identifiée par son itérateur
  inline Cell cellConnectivity(const CellEnumerator &ci) {
    return *ci;
  }

  //! Retourne l'objet permettant de récupérer les noeuds de la maille identifiée par elle-même !
  inline Cell cellConnectivity(const Cell &cell) {
    return cell;
  }

 private:
  Integer m_nb_node;  //! Nb de noeuds d'une maille quelconque pour la dimension de la grille cartésienne
};

}

#endif

