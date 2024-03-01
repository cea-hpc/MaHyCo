#ifndef ACCENV_INDEX_SELECTER_H
#define ACCENV_INDEX_SELECTER_H

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/accelerator/Filter.h"

#include "accenv/SingletonIAccEnv.h"

namespace Accenv {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*!
 * \brief Construction d'un sous-ensemble d'indexes à partir d'un critère
 * 
 */
class IndexSelecter {
 public:
  IndexSelecter(ISubDomain* sd) :
    m_sd (sd),
    // --------------------------------------------------------
    m_is_acc_avl (AcceleratorUtils::isAvailable(m_sd)),
    m_mem_h (m_is_acc_avl ? eMemoryRessource::HostPinned : eMemoryRessource::Host),
    m_mem_d (m_is_acc_avl ? eMemoryRessource::Device : eMemoryRessource::Host),
    m_lid_select_d (platform::getDataMemoryRessourceMng()->getAllocator(m_mem_d)),
    m_lid_select_h (platform::getDataMemoryRessourceMng()->getAllocator(m_mem_h))
  {
  }

  ~IndexSelecter() {}

  /*!
   * \brief Définit l'intervalle [0,nb_idx[ sur lequel va s'opérer la sélection
   */
  void resize(Int32 nb_idx) {
    m_nb_idx = nb_idx;
    m_lid_select_d.resize(m_nb_idx);
    m_lid_select_h.resize(m_nb_idx);
  }

  /*!
   * \brief Termine et synchronise la sélection pour l'environnment env_id
   * \return Si host_view, retourne une vue HOST sur les éléments sélectionnés, sinon vue DEVICE
   */
  template<typename PredicateType>
  ConstArrayView<Int32> syncSelectIf(Ref<RunQueue> rqueue_async, PredicateType pred, bool host_view=false) 
  {
    // On sélectionne dans [0,m_nb_idx[ les indices i pour lesquels pred(i) est vrai
    //  et on les copie dans out_lid_select.
    //  Le nb d'indices sélectionnés est donné par nbOutputElement()
    SmallSpan<Int32> out_lid_select(m_lid_select_d.data(), m_nb_idx);

    ax::GenericFilterer gen_filterer(rqueue_async.get());
    gen_filterer.applyWithIndex(m_nb_idx, pred,
	[=] ARCCORE_HOST_DEVICE (Int32 input_index, Int32 output_index) -> void
	{
	  out_lid_select[output_index] = input_index;
	});
    Int32 nb_idx_selected = gen_filterer.nbOutputElement();

    if (nb_idx_selected && host_view) 
    {
      prof_acc_begin("cpyD2H_lid_sel");
      // Copie asynchrone Device to Host (m_lid_select_d ==> m_lid_select_h)
      rqueue_async->copyMemory(ax::MemoryCopyArgs(m_lid_select_h.subView(0, nb_idx_selected).data(),
                                                  m_lid_select_d.subView(0, nb_idx_selected).data(),
                                                  nb_idx_selected * sizeof(Int32))
                                   .addAsync());

      rqueue_async->barrier();
      prof_acc_end("cpyD2H_lid_sel");
    }
    else
    {
      rqueue_async->barrier();
    }

    ConstArrayView<Int32> lid_select_view = (
      host_view ? 
      m_lid_select_h.subConstView(0, nb_idx_selected) : 
      m_lid_select_d.subConstView(0, nb_idx_selected));

    return lid_select_view;
  }

 private:
  ISubDomain* m_sd = nullptr;

  bool m_is_acc_avl = false;  // indique si l'accélérateur est disponible ou non
  eMemoryRessource m_mem_h; // identification de l'allocateur HOST
  eMemoryRessource m_mem_d; // identification de l'allocateur DEVICE
  UniqueArray<Int32> m_lid_select_d; // liste des identifiants sélectionnés avec un Filterer (alloué sur DEVICE)
  UniqueArray<Int32> m_lid_select_h; // liste des identifiants sélectionnés avec un Filterer (alloué sur HOST)

  Int32 m_nb_idx=0;  //!< Intervalle [0, m_nb_idx[ sur lequel on va opérer la sélection
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}  // namespace Accenv

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#endif
