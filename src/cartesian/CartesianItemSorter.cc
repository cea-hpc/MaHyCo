#include "cartesian/CartesianItemSorter.h"
#include "cartesian/CartesianFaceId.h"
#include "arcane/IItemFamily.h"
#include "arcane/IMesh.h"

namespace Cartesian {

CartesianItemSorter::CartesianItemSorter(IMesh* mesh) : m_mesh(mesh) {}
CartesianItemSorter::~CartesianItemSorter() {}

/*! \brief Trie les faces selon un ordre cartésien
*/
void CartesianItemSorter::sortFaces() {
  // MODIFICATION DES UNIQUE IDS DES FACES PUIS TRI PAR DEFAUT SELON LES UNIQUE IDS
  // CONSEQUENCE : LES IDS LOCAUX ET GLOBAUX SONT MODIFIES ET CARTESIENS
  Arcane::IItemFamily *face_family = m_mesh->faceFamily();

  // On construit de nouvels id globaux pour les faces
  Cartesian::CartesianFaceUniqueIdBuilder cart_face_uniq_id_builder(face_family);
  cart_face_uniq_id_builder.computeFaceUniqueId();

  // On retrie (par defaut, le tri se fait selon l'unique Id)
  face_family->compactItems(/*do_sort=*/true);    
}

}

