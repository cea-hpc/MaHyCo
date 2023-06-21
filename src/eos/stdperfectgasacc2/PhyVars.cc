#include "eos/stdperfectgasacc2/PhyVarType.h"
#include "eos/stdperfectgasacc2/PhyVars.h"

#include "arcane/VariableView.h"
#include "arcane/ServiceBuilder.h"
#include "arcane/IItemFamily.h"
#include "arcane/ItemVector.h"
#include "arcane/materials/IMeshMaterialSynchronizeBuffer.h"

#include <accenv/IAccEnv.h>

#include "accenv/AcceleratorUtils.h"
#include "accenv/ProfAcc.h"

#include <arcane/accelerator/RunCommandMaterialEnumerate.h>
#include <arcane/accelerator/MaterialVariableViews.h>
#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/core/MeshVariableScalarRef.h>
#include <arcane/core/ISubDomain.h>

namespace Stdperfectgasacc2 {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class ImplPhyVars {
 public:
  ImplPhyVars() 
  {
    m_buffer = Arcane::Materials::impl::makeOneBufferMeshMaterialSynchronizeBufferRef(eMemoryRessource::UnifiedMemory);
  }

  ~ImplPhyVars() {}

  void addPhyVar(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector)
  {
    prof_acc_begin(var.name().localstr());
    auto phy_var_data = Arccore::makeRef(new PhyMatVarData(var, mat_cell_vector));
    m_phy_vars.push_back(Arccore::makeRef(new PhyVarType(phy_var_data)));
    prof_acc_end(var.name().localstr());
  }
  
  void copyFromAllArcVars()
  {
    prof_acc_begin("copyFromAllArcVars");
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
    for(Arcane::Int32 i = 0 ; i<m_phy_vars.size() ; ++i) {
      m_phy_vars[i]->phy_var_data->copyVarToRawData();
    }
    prof_acc_end("copyFromAllArcVars");
  }

  Arcane::Span<Arcane::Real> rawData(Arcane::Int32 i)
  {
    return m_phy_vars[i]->phy_var_data->rawData();
  }

  void copyIntoArcVar(Arcane::Int32 i) const
  {
    auto phy_var_data = m_phy_vars[i]->phy_var_data;
    prof_acc_begin(phy_var_data->name().localstr());
    phy_var_data->copyRawDataToVar();
    prof_acc_end(phy_var_data->name().localstr());
  }

  void clear()
  {
    m_phy_vars.clear();
  }

 protected:
  std::vector< Arcane::Ref<PhyVarType> > m_phy_vars;
  Arcane::Ref<Arcane::Materials::IMeshMaterialSynchronizeBuffer> m_buffer;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

PhyVars::PhyVars() 
{
  ImplPhyVars* impl_phy_vars = new ImplPhyVars();
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

void PhyVars::copyFromAllArcVars() 
{
  get_ImplPhyVars_ptr(m_impl)->copyFromAllArcVars();
}

Arcane::Span<Arcane::Real> PhyVars::rawData(Arcane::Int32 i) 
{
  return get_ImplPhyVars_ptr(m_impl)->rawData(i);
}

void PhyVars::copyIntoArcVar(Arcane::Int32 i) const 
{
  get_ImplPhyVars_ptr(m_impl)->copyIntoArcVar(i);
}

void PhyVars::clear() 
{
  get_ImplPhyVars_ptr(m_impl)->clear();
}

} // Stdperfectgasacc2

