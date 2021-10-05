#ifndef CARTESIAN_CARTESIAN_FACE_ID_H
#define CARTESIAN_CARTESIAN_FACE_ID_H

#include "cartesian/CartTypes.h"
#include "cartesian/CartesianGridT.h"
#include "arcane/Item.h"
#include "arcane/IItemFamily.h"
#include "arcane/ItemUniqueId.h"
#include "arcane/utils/Array.h"
#include "arcane/Properties.h"
#include "arcane/IItemInternalSortFunction.h"
#include "arcane/IMesh.h"
#include "arcane/utils/ITraceMng.h"
#include "arcane/utils/String.h"
#include <algorithm>

namespace Cartesian {
  
/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Permet la construction d'unique ids cartésiens sur les faces
 */
/*---------------------------------------------------------------------------*/
class CartesianFaceUniqueIdBuilder {
 public:
  //!
  struct UniqueIdItem {
    UniqueIdType unique_id = -1;
    ItemInternal* item = nullptr;
  };

 public:
  CartesianFaceUniqueIdBuilder(IItemFamily* face_family) 
  : m_face_family(face_family)
  {
  }

  //! Pour affecter des uniques Id cartesiens sur les faces
  static void getFaceUniqueIdItems(IItemFamily* face_family, UniqueArray<UniqueIdItem> &id_items) {
    ItemInternalArrayView items(face_family->itemsInternal());

    Integer dim = face_family->mesh()->dimension();
    Integer nb_nodes_face = 1 << (dim-1); // Nb de noeuds sur une face
    
    auto trace_mng = face_family->traceMng();

    // Recuperation de la taille de la grille du domaine global
    Properties* mesh_properties = face_family->mesh()->properties();
    UniqueIdType3 glob_ncells_dir = {1, 1, 1};
    const char *str_prop[3] = {
      "GlobalNbCellX",
      "GlobalNbCellY",
      "GlobalNbCellZ"
    };
    for(Integer d(0) ; d < dim ; ++d) {
      glob_ncells_dir[d] = mesh_properties->getInt64WithDefault(String(str_prop[d]),-1);
    }

    // Grille cartesienne sur la globalite du domaine de calcul
    CartesianGridT<UniqueIdType> u_cart_grid(glob_ncells_dir, dim);

    // Numerotations cartesiennes globales pour les noeuds et les faces
    const auto &u_cart_num_node = u_cart_grid.cartNumNode();
    const auto &u_cart_num_face3 = u_cart_grid.cartNumFace3(); // u_cart_num_face3[dir] = numérotation dans la direction dir

    Integer nb_faces_to_do = items.size();
    id_items.resize(nb_faces_to_do); // dimensionnement du tableau en entrée id_items

    for( Integer i=0; i<nb_faces_to_do; ++i ) {
      ItemInternal* ii_cur = items[i];

      Face face(ii_cur);

      // On ne connait pas la direction de la face
      // En recuperant les noeuds et en identifiant la direction dans laquelle les coordonnees
      // sont identiques, on aura la direction normale de la face
      bool same_coord[3] = {true, true, true};
      UniqueIdType3 uniq_node_ijk[4];
      UniqueIdType min_node_unique_id = std::numeric_limits<UniqueIdType>::max();
      Integer min_inode = -1;
      Integer face_dir = -1;

      for(Integer inode = 0 ; inode < nb_nodes_face ; ++inode) {
        Node node(face.node(inode));
        UniqueIdType node_unique_id = node.uniqueId().asInt64();

        if (node_unique_id < min_node_unique_id) {
          min_node_unique_id = node_unique_id;
          min_inode = inode;
        }

        // Passage unique_id => (I,J,K) pour le noeud
        u_cart_num_node.ijk(node_unique_id, uniq_node_ijk[inode]);

        if (inode > 0) {
          for(Integer d(0) ; d < dim ; ++d) {
            if (uniq_node_ijk[/*inode=*/0][d] != uniq_node_ijk[inode][d]) {
              same_coord[d] = false;
            }
          }
        }
      }
      Integer nb_equal_dir = 0;
      for(Integer d(0) ; d < dim ; ++d) {
        if (same_coord[d]) {
          face_dir = d;
          nb_equal_dir++;
        }
      }
      if (nb_equal_dir != 1) {
        trace_mng->fatal() << "Impossible, la face n'a pas exactement une direction : " << nb_equal_dir;
      }

      // La face est de direction face_dir et ses coordonnes cartesiennes vont etre celles du noeud min_inode
      UniqueIdType3 uniq_face_ijk;
      uniq_face_ijk[0] = uniq_node_ijk[min_inode][0];
      uniq_face_ijk[1] = uniq_node_ijk[min_inode][1];
      uniq_face_ijk[2] = uniq_node_ijk[min_inode][2];

      // Passage (I,J,K) => unique_face_id
      auto unique_face_id = u_cart_num_face3[face_dir].id(uniq_face_ijk);

      // On cherche à récupérer ce couple
      id_items[i].unique_id = unique_face_id;
      id_items[i].item = ii_cur;
    }
  }

  //! Pour modifier les uniques Id cartesiens sur les faces
  void computeFaceUniqueId() {
    // Récupération des couples (unique_id, ptr_item)
    UniqueArray<UniqueIdItem> id_items;
    getFaceUniqueIdItems(m_face_family, id_items);

    // On ecrase avec les items tries
    Integer nb_faces_to_do = id_items.size();
    for( Integer i=0; i<nb_faces_to_do; ++i ) {
      id_items[i].item->setUniqueId(id_items[i].unique_id);
    }
  }

 private:

  IItemFamily* m_face_family;  // les faces
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Pour trier les faces locales selon une numerotation cartesienne globale (unique)
 */
/*---------------------------------------------------------------------------*/
class CartesianFaceSorter
: public IItemInternalSortFunction
{
  // Pour faire le tri selon unique_id
  using UniqueIdItem = CartesianFaceUniqueIdBuilder::UniqueIdItem;

 public:
  CartesianFaceSorter(IItemFamily* face_family) 
  : m_face_family(face_family), m_name("CartesianFaceSorter")
  {
  }

  virtual ~CartesianFaceSorter() {}

  const String& name() const override
  {
    return m_name;
  }

  void sortItems(ItemInternalMutableArrayView items) override
  {
    // Récupération des couples (unique_id, ptr_item)
    UniqueArray<UniqueIdItem> id_items;
    CartesianFaceUniqueIdBuilder::getFaceUniqueIdItems(m_face_family, id_items);

    // Le tri lui-meme
    std::sort(id_items.begin(), id_items.end(), 
        [](const UniqueIdItem &ii1, const UniqueIdItem &ii2) {
        return ii1.unique_id < ii2.unique_id;
        });

    // On ecrase avec les items par les items tries
    Integer nb_faces_to_do = items.size();
    for( Integer i=0; i<nb_faces_to_do; ++i ) {
      items[i] = id_items[i].item;
    }
  }

 private:

  IItemFamily* m_face_family;  // les faces
  String m_name;
};

}

#endif

