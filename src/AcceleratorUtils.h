#ifndef ACCELERATOR_UTILS_H
#define ACCELERATOR_UTILS_H

#include "arcane/IApplication.h"
#include "arcane/accelerator/Reduce.h"
#include "arcane/accelerator/Runner.h"
#include "arcane/accelerator/Views.h"
#include "arcane/accelerator/Accelerator.h"

/*---------------------------------------------------------------------------*/
/* Pour les accélérateurs                                                    */
/*---------------------------------------------------------------------------*/

namespace ax = Arcane::Accelerator;

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
/*---------------------------------------------------------------------------*/


#endif
