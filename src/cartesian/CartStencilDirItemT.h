#ifndef CARTESIAN_CART_STENCIL_DIR_ITEM_T_H
#define CARTESIAN_CART_STENCIL_DIR_ITEM_T_H

#include "cartesian/CartTypes.h"

namespace Cartesian {
 
/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Permet d'accéder aux items [MinLayer,MaxLayer] autour de l'item courant 
 * dans une direction dir
 * Ex pour [MinLayer=-3,MaxLayer=+4] : item courant => id0
 *              ----- ----- ----- ----- ----- ----- ----- -----
 *  localId    |  -1 |  -1 | im1 | id0 | ip1 | ip2 | ip3 | ip4 |  ---> dir
 *              ----- ----- ----- ----- ----- ----- ----- -----
 *  ilayer        -3    -2    -1     0    +1    +2    +3    +4
 *  valid                  [validMin() ,             validMax()]
 *
 * ilayer = 0 réfère à l'item courant. 
 * [MinLayer,MaxLayer] peut ne pas contenir 0 
 *   => l'intervalle est décalé par rapport à l'item
 */
/*---------------------------------------------------------------------------*/
template<typename ItemIdType, Integer MinLayer, Integer MaxLayer>
class AsymStencilDirItemT {
 public:
  static constexpr Integer min_layer = MinLayer;
  static constexpr Integer max_layer = MaxLayer;
  static constexpr Integer extent = MaxLayer-MinLayer+1;

  /*!
   * item_id : le local id de l'item de base
   * idx_dir : indice cartesien de l'item dans la direction dir
   * m_nitemsm1_dir : Nb d'items -1 dans la direction dir
   * delta_dir : +-delta a appliquer sur item_id pour passer a l'item suivant/precedent
   */
  ARCCORE_HOST_DEVICE AsymStencilDirItemT(LocalIdType item_id, LocalIdType idx_dir, LocalIdType nitemsm1_dir, LocalIdType delta_dir)
  {
    // Domaine de validité : [idx_min, idx_max] (en numérotation "idx_dir")
    //  [idx_min, idx_max] = [idx_dir+min_layer, idx_dir+max_layer] \inter [0, nitemsm1_dir]
    LocalIdType idx_min = std::max(0,idx_dir+min_layer);
    LocalIdType idx_max = std::min(nitemsm1_dir,idx_dir+max_layer);
    // Passage en numérotation relative au stencil dans [min_layer,max_layer]
    m_sten_min = idx_min - idx_dir; // \in [min_layer,0] si min_layer<=0
    m_sten_max = idx_max - idx_dir; // \in [0,max_layer] si max_layer>=0

    // On fait le calcul quitte à invalider après
    for(Integer ilayer=min_layer ; ilayer<=max_layer ; ilayer++) {
      LocalIdType value=item_id+ilayer*delta_dir;
      bool is_valid = (m_sten_min<=ilayer && ilayer<=m_sten_max);
      m_sten_id[ilayer-min_layer]=(is_valid ? value : -1);
    }
  }

  //! [MinLayer
  ARCCORE_HOST_DEVICE static constexpr Integer minLayer() {
    return min_layer;
  }

  //! MaxLayer]
  ARCCORE_HOST_DEVICE static constexpr Integer maxLayer() {
    return max_layer;
  }

  //! Nombre total d'items possiblement impliqués dans le stencil (item de base compris)
  // = maxLayer()-minLayer()+1
  ARCCORE_HOST_DEVICE Integer eXtent() const {
    return extent;
  }

  //! Nombre d'items valides dans le stencil ( <=eXtent() )
  // = validMax() - validMin() +1
  ARCCORE_HOST_DEVICE Integer nbValidItem() const {
    return m_sten_max-m_sten_min+1;
  }

  //! Local id de l'item à la position ilayer \in [minLayer(),maxLayer()]
  // par rapport à l'item de base
  // -1 si l'item n'existe
  ARCCORE_HOST_DEVICE ItemIdType operator()(Integer ilayer) const {
    /*
     *          ---- ---- ---- ---- ---- ---- ---- ----
     *  ilayer | -3 | -2 | -1 |  0 | +1 | +2 | +3 | +4 |  ---> dir
     *          ---- ---- ---- ---- ---- ---- ---- ----
     *  Ind:      0    1    2    3    4    5    6    7 
     */
#ifndef ARCCORE_DEVICE_CODE
    ARCANE_ASSERT(min_layer<=ilayer && ilayer<=max_layer, ("ilayer out of bounds"));
#endif
    return ItemIdType{m_sten_id[ilayer-min_layer]};
  }

  //! this->operator()(ilayer) == -1 ssi ilayer < validMin()
  ARCCORE_HOST_DEVICE Integer validMin() const {
    return m_sten_min;
  }

  //! this->operator()(ilayer) == -1 ssi ilayer > validMax()
  ARCCORE_HOST_DEVICE Integer validMax() const {
    return m_sten_max;
  }

 protected:

  //! Intervalle de validité [m_sten_min,m_sten_max] \in [min_layer,max_layer]
  Integer m_sten_min;
  Integer m_sten_max;

  LocalIdType m_sten_id[extent];  //! les ids locaux des items dans le stencil
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Permet d'accéder aux NLayer items autour de l'item courant dans une direction
 * Ex pour NLayer=3 : item courant => id0
 *              ----- ----- ----- ----- ----- ----- -----
 *  localId    |  -1 |  -1 | im1 | id0 | ip1 | ip2 | ip3 |  ---> dir
 *              ----- ----- ----- ----- ----- ----- -----
 *  ilayer        -3    -2    -1     0    +1    +2    +3
 *  valid                  [validMin() ,       validMax()]
 */
/*---------------------------------------------------------------------------*/
template<typename ItemIdType, Integer NLayer>
class CartStencilDirItemT 
: public AsymStencilDirItemT<ItemIdType,-NLayer,+NLayer> {

  using SuperType = AsymStencilDirItemT<ItemIdType,-NLayer,+NLayer>;
 public:
  /*!
   * item_id : le local id de l'item central
   * idx_dir : indice cartesien de l'item dans la direction dir
   * m_nitemsm1_dir : Nb d'items -1 dans la direction dir
   * delta_dir : +-delta a appliquer sur item_id pour passer a l'item suivant/precedent
   */
  ARCCORE_HOST_DEVICE CartStencilDirItemT(LocalIdType item_id, LocalIdType idx_dir, LocalIdType nitemsm1_dir, LocalIdType delta_dir)
  : SuperType(item_id, idx_dir, nitemsm1_dir, delta_dir)
  {
  }

  //! Nombre de couches d'items avant (ou après) l'item central
  ARCCORE_HOST_DEVICE static constexpr Integer nLayer() {
    return NLayer;
  }

  //! Local Id de l'item central
  ARCCORE_HOST_DEVICE ItemIdType centralId() const {
    return this->operator()(0);
  }

  //! Local id de l'item précédent et adjacent à itemId()
  ARCCORE_HOST_DEVICE ItemIdType previousId() const {
    return this->operator()(-1);
  }

  //! Local id de l'item suivant et adjacent à itemId()
  ARCCORE_HOST_DEVICE ItemIdType nextId() const {
    return this->operator()(+1);
  }
  
  /*
   *  C'est pas beau, je n'arrive pas à utiliser l'opérateur (ilayer) ...
   */ 
  
  //! Local id de l'item précédent et adjacent à itemId()
  ARCCORE_HOST_DEVICE ItemIdType prev_previousId() const {
    return this->operator()(-2);
  }
  
  //! Local id de l'item suivant et adjacent à itemId()
  ARCCORE_HOST_DEVICE ItemIdType next_nextId() const {
    return this->operator()(+2);
  }
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Pour une diretion dir donnée,
 * permet d'accéder aux NLayer items avant l'item de référence (base)
 * et aux NLayer-1 items après l'item courant
 * Ex pour NLayer=3 : item courant => id0
 *              ----- ----- -----|----- ----- ----- 
 *  localId    |  -1 |  -1 | im1 | id0 | ip1 | ip2 |  ---> dir
 *              ----- ----- -----|----- ----- ----- 
 *  ilayer        -3    -2    -1 |  +1    +2    +3
 *            previousId(ilayer) |   nextId(ilayer)
 *
 * L'item courant correspond à nextId(+1) d'où Pos (Positive)
 */
/*---------------------------------------------------------------------------*/
template<typename ItemIdType, Integer NLayer>
class PosAsymStencilDirItemT 
: public AsymStencilDirItemT<ItemIdType,-NLayer,+NLayer-1> {

  using SuperType = AsymStencilDirItemT<ItemIdType,-NLayer,+NLayer-1>;
  using SuperType::m_sten_min;
  using SuperType::m_sten_max;
 public:
  /*!
   * item_id : le local id de l'item de référence (base)
   * idx_dir : indice cartesien de l'item dans la direction dir
   * m_nitemsm1_dir : Nb d'items -1 dans la direction dir
   * delta_dir : +-delta a appliquer sur item_id pour passer a l'item suivant/precedent
   */
  ARCCORE_HOST_DEVICE PosAsymStencilDirItemT(LocalIdType item_id, LocalIdType idx_dir, LocalIdType nitemsm1_dir, LocalIdType delta_dir)
  : SuperType(item_id, idx_dir, nitemsm1_dir, delta_dir)
  {
  }

  //! Nombre de couches d'items avant (ou après)
  ARCCORE_HOST_DEVICE static constexpr Integer nLayer() {
    return NLayer;
  }

  //! Item de référence qui a servi de base
  ARCCORE_HOST_DEVICE ItemIdType baseId(Integer ilayer) const {
    return this->operator()(0);
  }

  //! Local id d'un item précédent : ilayer \in [-NLayer,-1]
  ARCCORE_HOST_DEVICE ItemIdType previousId(Integer ilayer) const {
#ifndef ARCCORE_DEVICE_CODE
    ARCANE_ASSERT(-NLayer<=ilayer && ilayer<=-1, ("ilayer not in [-NLayer,-1]"));
#endif
    return this->operator()(ilayer);
  }

  //! Local id d'un item suivant : ilayer \in [+1,+NLayer]
  ARCCORE_HOST_DEVICE ItemIdType nextId(Integer ilayer) const {
#ifndef ARCCORE_DEVICE_CODE
    ARCANE_ASSERT(1<=ilayer && ilayer<=NLayer, ("ilayer not in [+1,+NLayer]"));
#endif
    return this->operator()(ilayer-1);
  }

  //! this->previousId()(ilayer) == -1 ssi ilayer < validMin()
  ARCCORE_HOST_DEVICE Integer validMin() const {
    return m_sten_min;
  }

  //! this->nextId()(ilayer) == -1 ssi ilayer > validMax()
  ARCCORE_HOST_DEVICE Integer validMax() const {
    return m_sten_max+1;
  }
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Pour une direction dir donnée,
 * permet d'accéder aux NLayer-1 items avant l'item de référence (base)
 * et aux NLayer items après l'item courant
 * Ex pour NLayer=3 : item courant => id0
 *              ----- ----- -----|----- ----- ----- 
 *  localId    |  -1 | im1 | id0 | ip1 | ip2 | ip3 |  ---> dir
 *              ----- ----- -----|----- ----- ----- 
 *  ilayer        -3    -2    -1 |  +1    +2    +3
 *            previousId(ilayer) |   nextId(ilayer)
 *
 * L'item courant correspond à previousId(-1) d'où Neg (Negative)
 */
/*---------------------------------------------------------------------------*/
template<typename ItemIdType, Integer NLayer>
class NegAsymStencilDirItemT 
: public AsymStencilDirItemT<ItemIdType,-NLayer+1,+NLayer> {

  using SuperType = AsymStencilDirItemT<ItemIdType,-NLayer+1,+NLayer>;
  using SuperType::m_sten_min;
  using SuperType::m_sten_max;
 public:
  /*!
   * item_id : le local id de l'item de référence (base)
   * idx_dir : indice cartesien de l'item dans la direction dir
   * m_nitemsm1_dir : Nb d'items -1 dans la direction dir
   * delta_dir : +-delta a appliquer sur item_id pour passer a l'item suivant/precedent
   */
  ARCCORE_HOST_DEVICE NegAsymStencilDirItemT(LocalIdType item_id, LocalIdType idx_dir, LocalIdType nitemsm1_dir, LocalIdType delta_dir)
  : SuperType(item_id, idx_dir, nitemsm1_dir, delta_dir)
  {
  }

  //! Nombre de couches d'items avant (ou après)
  ARCCORE_HOST_DEVICE static constexpr Integer nLayer() {
    return NLayer;
  }

  //! Item de référence qui a servi de base
  ARCCORE_HOST_DEVICE ItemIdType baseId(Integer ilayer) const {
    return this->operator()(0);
  }

  //! Local id d'un item précédent : ilayer \in [-NLayer,-1]
  ARCCORE_HOST_DEVICE ItemIdType previousId(Integer ilayer) const {
#ifndef ARCCORE_DEVICE_CODE
    ARCANE_ASSERT(-NLayer<=ilayer && ilayer<=-1, ("ilayer not in [-NLayer,-1]"));
#endif
    return this->operator()(ilayer+1);
  }

  //! Local id d'un item suivant : ilayer \in [+1,+NLayer]
  ARCCORE_HOST_DEVICE ItemIdType nextId(Integer ilayer) const {
#ifndef ARCCORE_DEVICE_CODE
    ARCANE_ASSERT(1<=ilayer && ilayer<=NLayer, ("ilayer not in [+1,+NLayer]"));
#endif
    return this->operator()(ilayer);
  }

  //! this->previousId()(ilayer) == -1 ssi ilayer < validMin()
  ARCCORE_HOST_DEVICE Integer validMin() const {
    return m_sten_min-1;
  }

  //! this->nextId()(ilayer) == -1 ssi ilayer > validMax()
  ARCCORE_HOST_DEVICE Integer validMax() const {
    return m_sten_max;
  }
};

}

#endif
