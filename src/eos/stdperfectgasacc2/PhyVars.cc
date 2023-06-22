#include "eos/stdperfectgasacc2/PhyVarType.h"
#include "eos/stdperfectgasacc2/PhyVars.h"

#include "arcane/VariableView.h"
#include "arcane/ServiceBuilder.h"
#include "arcane/IItemFamily.h"
#include "arcane/ItemVector.h"
#include "arcane/materials/IMeshMaterialSynchronizeBuffer.h"

#include <accenv/IAccEnv.h>

#include "accenv/AcceleratorUtils.h"
#include "accenv/SingletonIAccEnv.h"
#include "accenv/ProfAcc.h"

#include <arcane/accelerator/RunCommandMaterialEnumerate.h>
#include <arcane/accelerator/MaterialVariableViews.h>
#include <arcane/accelerator/AsyncRunQueuePool.h>
#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/core/MeshVariableScalarRef.h>
#include <arcane/core/ISubDomain.h>

namespace Stdperfectgasacc2 {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class ImplPhyVars {
 public:
  ImplPhyVars(ax::Runner& runner) :
    m_runner (runner),
    m_async_queue_pool (runner)
  {
    m_buffer = Arcane::Materials::impl::makeOneBufferMeshMaterialSynchronizeBufferRef(eMemoryRessource::UnifiedMemory);
    m_main_event = ax::makeEventRef(m_runner);
  }

  ~ImplPhyVars() {}

  void addPhyVar(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector)
  {
    prof_acc_begin(var.name().localstr());
    auto phy_var_data = Arccore::makeRef(new PhyMatVarData(var, mat_cell_vector));
    m_phy_vars.push_back(Arccore::makeRef(new PhyVarType(phy_var_data)));
    prof_acc_end(var.name().localstr());
  }
  
  void asyncCopyFromAllArcVars(ax::RunQueue* async_queue)
  {
    prof_acc_begin("asyncCopyFromAllArcVars");
    // On décompte toutes les variables pour pré-allouer
    m_buffer->setNbRank(m_phy_vars.size());
    for(Arcane::Int32 i = 0 ; i<m_phy_vars.size() ; ++i) {
      m_buffer->setSendBufferSize(i, m_phy_vars[i]->phy_var_data->sizeInBytes());
    }
    m_buffer->allocate();

    // On récupère les zones mémoires déjà allouées
    for(Arcane::Int32 i = 0 ; i<m_phy_vars.size() ; ++i) {
      m_phy_vars[i]->phy_var_data->setBuffer(m_buffer->sendBuffer(i));
    }

    // On remplit les rawData()

    // On pré-construit les évenements manquants si nécessaire
    for(Arcane::Int32 i = m_event.size() ; i<m_phy_vars.size() ; ++i) {
      m_event.add(ax::makeEventRef(m_runner));
    }

    for(Arcane::Int32 i = 0 ; i<m_phy_vars.size() ; ++i) 
    {
      m_phy_vars[i]->phy_var_data->asyncCopyVarToRawData(m_async_queue_pool[i]);
    }
    for(Arcane::Int32 i = 0 ; i<m_phy_vars.size() ; ++i) 
    {
      m_async_queue_pool[i]->recordEvent(m_event[i]); // m_event[i] = événement de fin de copie sur m_async_queue_pool[i]
      async_queue->waitEvent(m_event[i]); // la queue principale async_queue doit se synchroniser avec m_event[i]
    }
    prof_acc_end("asyncCopyFromAllArcVars");
  }

  Arcane::Span<Arcane::Real> rawData(Arcane::Int32 i)
  {
    return m_phy_vars[i]->phy_var_data->rawData();
  }

  void asyncCopyIntoArcVar(ax::RunQueue* async_queue, Arcane::Int32 i) const
  {
    auto phy_var_data = m_phy_vars[i]->phy_var_data;
    prof_acc_begin(phy_var_data->name().localstr());
    phy_var_data->asyncCopyRawDataToVar(async_queue);
    prof_acc_end(phy_var_data->name().localstr());
  }

  void asyncCopyIntoArcVarList(ax::RunQueue* async_queue, std::initializer_list<Arcane::Int32> li) 
  {
    prof_acc_begin("asyncCopyIntoArcVarList");

    async_queue->recordEvent(m_main_event); // événement de toutes les calculs précédents sur async_queue
    for(Arcane::Int32 i : li) 
    {
      m_async_queue_pool[i]->waitEvent(m_main_event); // permet à m_async_queue_pool[i] de se synchroniser avec async_queue

      m_phy_vars[i]->phy_var_data->asyncCopyRawDataToVar(m_async_queue_pool[i]);
    }
    for(Arcane::Int32 i : li) 
    {
      m_async_queue_pool[i]->recordEvent(m_event[i]); // m_event[i] = événement de fin de copie sur m_async_queue_pool[i]
      async_queue->waitEvent(m_event[i]); // la queue principale async_queue doit se synchroniser avec m_event[i]
    }

    prof_acc_end("asyncCopyIntoArcVarList");
  }

  void clear()
  {
    m_phy_vars.clear();
  }

 protected:
  ax::Runner& m_runner;
  std::vector< Arcane::Ref<PhyVarType> > m_phy_vars;
  Arcane::Ref<Arcane::Materials::IMeshMaterialSynchronizeBuffer> m_buffer;
  ax::AsyncRunQueuePool m_async_queue_pool;
  Arcane::UniqueArray<Arcane::Ref<ax::RunQueueEvent> > m_event;
  Arcane::Ref<ax::RunQueueEvent> m_main_event;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

PhyVars::PhyVars(ax::Runner& runner) 
{
  ImplPhyVars* impl_phy_vars = new ImplPhyVars(runner);
  m_impl = static_cast<void*>(impl_phy_vars);
}

ImplPhyVars* get_ImplPhyVars_ptr(void* impl) {
  return reinterpret_cast<ImplPhyVars*>(impl);
}

PhyVars::~PhyVars() 
{
  ImplPhyVars* impl_phy_vars = get_ImplPhyVars_ptr(m_impl);
  delete impl_phy_vars;
  m_impl = nullptr;
}

void PhyVars::addPhyVar(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector) 
{
  get_ImplPhyVars_ptr(m_impl)->addPhyVar(var, mat_cell_vector);
}

void PhyVars::asyncCopyFromAllArcVars(ax::RunQueue* async_queue) 
{
  get_ImplPhyVars_ptr(m_impl)->asyncCopyFromAllArcVars(async_queue);
}

Arcane::Span<Arcane::Real> PhyVars::rawData(Arcane::Int32 i) 
{
  return get_ImplPhyVars_ptr(m_impl)->rawData(i);
}

void PhyVars::asyncCopyIntoArcVar(ax::RunQueue* async_queue, Arcane::Int32 i) const 
{
  get_ImplPhyVars_ptr(m_impl)->asyncCopyIntoArcVar(async_queue, i);
}

void PhyVars::asyncCopyIntoArcVarList(ax::RunQueue* async_queue, std::initializer_list<Arcane::Int32> li)  
{
  get_ImplPhyVars_ptr(m_impl)->asyncCopyIntoArcVarList(async_queue, li);
}

void PhyVars::clear() 
{
  get_ImplPhyVars_ptr(m_impl)->clear();
}

} // Stdperfectgasacc2

