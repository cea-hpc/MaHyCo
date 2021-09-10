#ifndef ACCELERATOR_UTILS_H
#define ACCELERATOR_UTILS_H

#include "arcane/IApplication.h"
#include "arcane/accelerator/Reduce.h"
#include "arcane/accelerator/Runner.h"
#include "arcane/accelerator/Views.h"
#include "arcane/accelerator/Accelerator.h"
#include "arcane/accelerator/RunCommandLoop.h"
#include "arcane/accelerator/RunCommandEnumerate.h"

/*---------------------------------------------------------------------------*/
/* Pour les accélérateurs                                                    */
/*---------------------------------------------------------------------------*/

namespace ax = Arcane::Accelerator;
#if 0
template<typename ItemType>
class ItemRunCommand
{
 public:
  ItemRunCommand(ax::RunCommand& command,const ItemVectorViewT<ItemType>& items)
  : m_command(command), m_items(items)
  {
  }
  ax::RunCommand& m_command;
  ItemVectorViewT<ItemType> m_items;
};

template<typename ItemType> ItemRunCommand<ItemType>
operator<<(ax::RunCommand& command,const ItemGroupT<ItemType>& items)
{
  return ItemRunCommand<ItemType>(command,items.view());
}

template<typename ItemType> ItemRunCommand<ItemType>
operator<<(ax::RunCommand& command,const ItemVectorViewT<ItemType>& items)
{
  return ItemRunCommand<ItemType>(command,items);
}

template<typename ItemType,typename Lambda>
void operator<<(ItemRunCommand<ItemType>&& nr,Lambda f)
{
  run(nr.m_command,nr.m_items,std::forward<Lambda>(f));
}
template<typename ItemType,typename Lambda>
void operator<<(ItemRunCommand<ItemType>& nr,Lambda f)
{
  run(nr.m_command,nr.m_items,std::forward<Lambda>(f));
}

#define RUNCOMMAND_ENUMERATE(ItemNameType,iter_name,item_group)  \
  item_group << [=] ARCCORE_HOST_DEVICE (ItemNameType##LocalId iter_name)
#endif
/*---------------------------------------------------------------------------*/
/* Pour le profiling sur accélérateur                                        */
/*---------------------------------------------------------------------------*/

#if defined(ARCANE_COMPILING_CUDA) && defined(PROF_ACC)
#define USE_PROF_ACC
#endif

#if defined(USE_PROF_ACC)

#warning "PROF_ACC : instrumentation avec nvtx"
#include <nvtx3/nvToolsExt.h>

#ifndef PROF_ACC_BEGIN
#define PROF_ACC_BEGIN(__name__) nvtxRangePushA(__name__)
#endif

#ifndef PROF_ACC_END
#define PROF_ACC_END nvtxRangePop()
#endif

#else

//#warning "Pas d'instrumentation"
#ifndef PROF_ACC_BEGIN
#define PROF_ACC_BEGIN(__name__)
#endif

#ifndef PROF_ACC_END
#define PROF_ACC_END 
#endif

#endif

/*---------------------------------------------------------------------------*/
/* Pour créer une vue sur les valeurs d'un environnement                     */
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::Materials;

template<typename value_type>
ArrayView<value_type> envView(CellMaterialVariableScalarRef<value_type>& var_menv, IMeshEnvironment* env) {
  // Pour rappel, [0] fait référence à .globalVariable() d'où le +1
  return var_menv._internalValue()[env->id()+1];
}

/*---------------------------------------------------------------------------*/
/* Pour gérer un nombre dynamique de RunQueue asynchrones                    */
/*---------------------------------------------------------------------------*/
class MultiAsyncRunQueue {
 public:

  MultiAsyncRunQueue(ax::Runner& runner, Integer asked_nb_queue) {
    // au plus 32 queues (32 = nb de kernels max exécutables simultanément)
    m_nb_queue = std::min(asked_nb_queue, 32);
    m_queues.resize(m_nb_queue);
    for(Integer iq=0 ; iq<m_nb_queue ; ++iq) {
      m_queues[iq] = new ax::RunQueue(runner);
      m_queues[iq]->setAsync(true);
    }
  }

  virtual ~MultiAsyncRunQueue() {
    for(auto q : m_queues) {
      delete q;
    }
    m_nb_queue=0;
    m_queues.clear();
  }

  // Pour récupérer la iq%nbQueue()-ième queue d'exécution
  inline ax::RunQueue& queue(Integer iq) {
    return *(m_queues[iq%m_nb_queue]);
  }

  // Force l'attente de toutes les RunQueue
  void waitAllQueues() {
    for(auto q : m_queues) {
      q->barrier();
    }
  }

  // Nombre de RunQueue asynchrones
  inline Integer nbQueue() const {
    ARCANE_ASSERT(m_nb_queue==m_queues.size(), ("m_nb_queue est different de m_queues.size()"));
    return m_nb_queue;
  }

 protected:
  UniqueArray<ax::RunQueue*> m_queues; //!< toutes les RunQueue
  Integer m_nb_queue=0; //!< m_queues.size()
};

#ifdef ARCANE_HAS_CUDA
/*---------------------------------------------------------------------------*/
/* Pour indiquer que le contenu du tableau est accéder fréquemment           */
/*---------------------------------------------------------------------------*/
template<typename ViewType>
void mem_adv_set_read_mostly(ViewType view, int device) {
  cudaMemAdvise (view.data(), view.size(), cudaMemAdviseSetReadMostly,device);
}
#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


#endif
