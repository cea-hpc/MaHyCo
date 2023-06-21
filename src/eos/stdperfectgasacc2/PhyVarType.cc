#include "eos/stdperfectgasacc2/PhyVarType.h"
#include "accenv/AcceleratorUtils.h"
#include "accenv/ProfAcc.h"

#include <arcane/accelerator/RunCommandMaterialEnumerate.h>
#include <arcane/accelerator/MaterialVariableViews.h>
#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/core/MeshVariableScalarRef.h>
#include <arcane/core/ISubDomain.h>

namespace Stdperfectgasacc2 {

/*---------------------------------------------------------------------------*/
/* To access to MatCellVectorView on accelerator */
/*---------------------------------------------------------------------------*/
class MatCellContainer : 
  public Accelerator::impl::MatCommandContainerBase
{
 public:
  using MatCellVectorView = Arcane::Materials::MatCellVectorView;
  using ComponentItemVectorView = Arcane::Materials::ComponentItemVectorView;
  using IMeshEnvironment = Arcane::Materials::IMeshEnvironment;
  using ComponentItemLocalId = Arcane::Materials::ComponentItemLocalId;
  using MatVarIndex = Arcane::Materials::MatVarIndex;

 public:
  explicit MatCellContainer(Arcane::Materials::MatCellVectorView view) :
    Accelerator::impl::MatCommandContainerBase(view)
  {
  }

  //! Accesseur pour le i-ème élément de la liste
  constexpr ARCCORE_HOST_DEVICE ComponentItemLocalId operator[](Arcane::Int32 i) const
  {   
    return { ComponentItemLocalId(m_matvar_indexes[i]) };
  }   
};

/*---------------------------------------------------------------------------*/
/* PhyMatVarData */
/*---------------------------------------------------------------------------*/
PhyMatVarData::
PhyMatVarData(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector) : 
  IPhyVarData(),
  m_var (var),
  m_mat_cell_vector (mat_cell_vector)
{
}

PhyMatVarData::
~PhyMatVarData() {
}

Arcane::String PhyMatVarData::name() const {
  return m_var.name();
}

Arcane::Int32 PhyMatVarData::sizeInBytes() const {
  return m_mat_cell_vector->nbItem() * sizeof(Arcane::Real);
}

void PhyMatVarData::
setBuffer(Arcane::Span<Arcane::Byte> span_byte) {
  if (span_byte.size()%sizeof(Arcane::Real) != 0)
    throw FatalErrorException(A_FUNCINFO, "PhyMatVarData::setBuffer(...) : span_byte.size() n'est pas un multiple de sizeof(Arcane::Real)");
  m_raw_data = Arcane::Span<Arcane::Real>(reinterpret_cast<Arcane::Real*>(span_byte.data()), span_byte.size()/sizeof(Arcane::Real));
}

// Arcane var => raw_data, implem GPU API
void PhyMatVarData::
copyVarToRawData() 
{
  PROF_ACC_BEGIN(__FUNCTION__);

  auto command = ax::makeCommand(m_var.globalVariable().subDomain()->acceleratorMng()->defaultQueue());
  auto in_var = ax::viewIn(command, m_var);
  Arcane::Span<Arcane::Real> out_raw_data(m_raw_data);
  MatCellContainer mat_cell_cont(*m_mat_cell_vector);

  auto nelt = mat_cell_cont.size();

  command << RUNCOMMAND_LOOP1(iter, nelt) 
  {
    auto i = iter()[0];

    auto mc = mat_cell_cont[i];
    out_raw_data[i] = in_var[mc];
  };

  PROF_ACC_END;
}

// raw_data => Arcane var, implem GPU API
void PhyMatVarData::
copyRawDataToVar() 
{
  PROF_ACC_BEGIN(__FUNCTION__);

  auto command = ax::makeCommand(m_var.globalVariable().subDomain()->acceleratorMng()->defaultQueue());
  auto out_var = ax::viewOut(command, m_var);
  Arcane::Span<const Arcane::Real> in_raw_data(m_raw_data);
  MatCellContainer mat_cell_cont(*m_mat_cell_vector);

  auto nelt = mat_cell_cont.size();

  command << RUNCOMMAND_LOOP1(iter, nelt) 
  {
    auto i = iter()[0];

    auto mc = mat_cell_cont[i];
    out_var[mc] = in_raw_data[i];
  };

  PROF_ACC_END;
}

PhyVarRawData PhyMatVarData::
rawData() {
  return m_raw_data;
}

/*---------------------------------------------------------------------------*/
/*  */
/*---------------------------------------------------------------------------*/

PhyVarType::
PhyVarType()
{
  prof_acc_mark("DfltCtor");
}

/* Accelerator variant */
PhyVarType::
PhyVarType(Arcane::Ref<IPhyVarData> a_phy_var_data) :
  phy_var_data (a_phy_var_data)
{
  prof_acc_mark("Ctor2");
}

PhyVarType::
PhyVarType(const PhyVarType& other) :
  phy_var_data (other.phy_var_data)
{
  prof_acc_mark("CpyCtor");
}

PhyVarType::
~PhyVarType()
{
  prof_acc_mark("Dtor");
}

PhyVarType& PhyVarType::
operator=(const PhyVarType& other)
{
  prof_acc_mark("CpyOp");
  phy_var_data = other.phy_var_data;

  return *this;
}

} // Stdperfectgasacc2

