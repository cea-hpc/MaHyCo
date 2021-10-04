#ifndef CARTESIAN_CART_ITEM_ENUMERATOR_T_H
#define CARTESIAN_CART_ITEM_ENUMERATOR_T_H

#include "cartesian/CartesianNumberingT.h"
#include "cartesian/CartesianGridT.h"
#include "cartesian/NumberingConverterT.h"
#include "arcane/Item.h"


/*!
 * \brief Devrait être dans static const bool Arcane::Item::null(LocalIdType local_id)
 */
class ItemId {
 public:
  ARCCORE_HOST_DEVICE static bool null(LocalIdType local_id) {
    return local_id == -1; // NULL_ITEM_ID
  }
};

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * A partir d'un ITEM_TYPE, définit les types complémentaires ItemType1 et ItemType2
 * Il faut spécialier pour chaque ITEM_TYPE
 */
/*---------------------------------------------------------------------------*/
// TODO : cas général à ne jamais utiliser, comment interdire son usage ?
template<typename ITEM_TYPE>
class TypesNumbConv {
 public:
  using ItemType1 = ITEM_TYPE;
  using ItemType2 = ITEM_TYPE;
};

// Spécialisations, à utiliser
template<>
class TypesNumbConv<Cell> {
 public:
  using ItemType1 = Face;
  using ItemType2 = Node;
};

template<>
class TypesNumbConv<Face> {
 public:
  using ItemType1 = Cell;
  using ItemType2 = Node;
};

template<>
class TypesNumbConv<Node> {
 public:
  using ItemType1 = Cell;
  using ItemType2 = Face;
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Implementation d'un iterateur sur item cartesien
 */
/*---------------------------------------------------------------------------*/
template<typename ITEM_TYPE>
class CartItemEnumeratorT {
 public:
  typedef typename ITEM_TYPE::LocalIdType ItemLocalIdType;

  //! Type de grille cartésienne sur des ids locaux
  using CartesianGrid = CartesianGridT<LocalIdType>;

  //! Type de numérotation cartésienne cohérente avec CartesianGrid
  using CartesianNumbering = typename CartesianGrid::CartesianNumbering;

  //! Type tableau sur pointeurs d'ItemInternal (implémentation d'un Item)
  using ItemInternalPtr = Item::ItemInternalPtr;

 private:
  using ItemType1 = typename TypesNumbConv<ITEM_TYPE>::ItemType1;
  using ItemType2 = typename TypesNumbConv<ITEM_TYPE>::ItemType2;

  using NumbConv1 = NumberingConverterT<ITEM_TYPE, ItemType1>;
  using NumbConv2 = NumberingConverterT<ITEM_TYPE, ItemType2>;

 public:
  // Affectation de l'itérateur
  CartItemEnumeratorT(const ItemInternalPtr* internals, Integer dir, 
    const CartesianGrid &cart_grid, const CartesianNumbering &cart_numb, 
    const LocalIdType3 &beg, const LocalIdType3 &end)
  : m_internals (internals), 
  m_dir (dir),
  m_cart_numbering (cart_numb), 
  m_beg(beg), m_end(end),
  m_numb_conv1 (dir, cart_grid),
  m_numb_conv2 (dir, cart_grid) {
    
    // On initialise les indices de boucles proprement dit
    m_item_ijk[0] = m_beg[0];
    m_item_ijk[1] = m_beg[1];
    m_item_ijk[2] = m_beg[2];

    // On calcule directement le premier id, ensuite on incrementera
    m_item_id = m_cart_numbering.id(m_item_ijk);

    // Initialisation des convertisseurs pour les types d'items complémenentaires
    m_numb_conv1.initDelta();
    m_numb_conv2.initDelta();
  }

  bool hasNext() const {
    bool is_last = (m_item_ijk[0] == m_beg[0] && m_item_ijk[1] == m_beg[1] && m_item_ijk[2] == m_end[2]);
    return !is_last;
  }

  void operator++() {
    m_item_ijk[0]++;  // i++
    m_item_id++;

    if (m_item_ijk[0] == m_end[0]) {  // Au bout de la ligne
      m_item_ijk[0] = m_beg[0];  // i = m_beg[0], on revient au debut de la ligne
      m_item_ijk[1]++;  // j++ , passage a la ligne suivante
      m_item_id = m_cart_numbering.id(m_item_ijk); // on recalcule l'id, plus simple
      m_numb_conv1.updateDelta(m_item_ijk[1], m_item_ijk[2]);
      m_numb_conv2.updateDelta(m_item_ijk[1], m_item_ijk[2]);

      if (m_item_ijk[1] == m_end[1]) {  // Au bout de la ligne et de la colonne
        m_item_ijk[1] = m_beg[1];  // j = m_beg[1], on revient au debut de la colonne
        m_item_ijk[2]++;  // k++ , passage au plan suivant
        m_item_id = m_cart_numbering.id(m_item_ijk); // on recalcule l'id, plus simple
        m_numb_conv1.updateDelta(m_item_ijk[1], m_item_ijk[2]);
        m_numb_conv2.updateDelta(m_item_ijk[1], m_item_ijk[2]);
      }
    }
  }

  //! Construit l'instance de ITEM_TYPE pour l'item sur lequel pointe l'itérateur
  ITEM_TYPE operator*() const {
    return ITEM_TYPE(m_internals, m_item_id); 
  }

  //! Le local id encapsulé selon le type de l'item
  inline ItemLocalIdType itemLocalId() const {
    return ItemLocalIdType(m_item_id);
  }

  //! Le local id (entier) de l'item courant pointé par l'énumérateur (=itemLocalId().localId())
  inline LocalIdType localId() const {
    return m_item_id;
  }

  //! Retourne le triplet (i,j,k) repérant l'item cartésien en cours
  inline const LocalIdType3 &itemIdx() const {
    return m_item_ijk;
  }

  //! Retourne itemIdx()[dir]
  inline LocalIdType itemIdxDir(Integer dir) const {
    return m_item_ijk[dir];
  }

  //! Récupération du local id correspondant à ItemType1
  inline LocalIdType localIdConv(ItemType1*) const {
    return m_item_id + m_numb_conv1.delta();
  }

  //! Récupération du local id correspondant à ItemType2
  inline LocalIdType localIdConv(ItemType2*) const {
    return m_item_id + m_numb_conv2.delta();
  }

 protected:
  const ItemInternalPtr* m_internals;  //! Tableau dimensionne au nb total d'items ITEM_TYPE, chaque case pointe vers un ItemInternal
  Integer m_dir; //! Direction privilégiée dans laquelle on peut demander des items voisins
  const CartesianNumbering &m_cart_numbering; //! Passage id <=> (i,j,k)

  //! Défintion de l'intervalle sur lequel itérer
  const LocalIdType3 &m_beg;
  const LocalIdType3 &m_end;

  NumbConv1 m_numb_conv1; //! Permet de passer de la numérotation de ITEM_TYPE vers ItemType1
  NumbConv2 m_numb_conv2; //! Permet de passer de la numérotation de ITEM_TYPE vers ItemType2

  LocalIdType m_item_id; //! L'identifiant local de l'item courant
  LocalIdType3 m_item_ijk; //! Les indices cartésiens (triplet) de l'item courant
};

typedef CartItemEnumeratorT<Cell> CartCellEnumerator;
typedef CartItemEnumeratorT<Face> CartFaceEnumerator;
typedef CartItemEnumeratorT<Node> CartNodeEnumerator;

}

#endif

