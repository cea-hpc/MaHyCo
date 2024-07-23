// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#include "eos/stdperfectgasacc1/PhyVarType.h"
#include "accenv/AcceleratorUtils.h"
#include "accenv/ProfAcc.h"

#include <arcane/accelerator/RunCommandMaterialEnumerate.h>
#include <arcane/accelerator/MaterialVariableViews.h>
#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/core/MeshVariableScalarRef.h>
#include <arcane/core/ISubDomain.h>

namespace Stdperfectgasacc1 {

#ifdef PHY_VAR_RAW_DATA_IS_STD

//#warning "__raw_data__ is a std::vector<Real>"
#define RAW_DATA_CTOR(__raw_data__, __size__) __raw_data__(__size__)

#else
Arccore::MemoryAllocationOptions alloc_opts() {
  Arcane::IMemoryAllocator* allocator = Arcane::platform::getDefaultDataAllocator();
  //Arcane::IMemoryAllocator* allocator = Arcane::platform::getDataMemoryRessourceMng()->getAllocator(Arcane::eMemoryRessource::Device);
  return Arccore::MemoryAllocationOptions(
      allocator,
      Arccore::eMemoryLocationHint::MainlyDevice);
}
//#warning "__raw_data__ is an Arcane::UniqueArray<Real>"
#define RAW_DATA_CTOR(__raw_data__, __size__) __raw_data__(alloc_opts(), __size__)

#endif

PhyVarType::
PhyVarType()
{
  prof_acc_mark("DfltCtor");
}

PhyVarType::
PhyVarType(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector) :
  var_name (var.name()),
  RAW_DATA_CTOR (raw_data, mat_cell_vector->size())
{
  prof_acc_mark("Ctor1");
  // Arcane var => raw_data
  copyVarToRawData(var, mat_cell_vector);
}

/* Accelerator variant */
PhyVarType::
PhyVarType(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector) :
  var_name (var.name()),
  RAW_DATA_CTOR (raw_data, mat_cell_vector->nbItem())
{
  prof_acc_mark("Ctor2");
  // Arcane var => raw_data
  copyVarToRawData(var, mat_cell_vector);
}

PhyVarType::
PhyVarType(const PhyVarType& other) :
  var_name (other.var_name),
  raw_data (other.raw_data)
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
  var_name = other.var_name;
  raw_data = other.raw_data;

  return *this;
}

// Arcane var => raw_data
void PhyVarType::
copyVarToRawData(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector)
{
  Arcane::Integer i(0);

  for(auto& mat_cell_i : *mat_cell_vector) {
    raw_data[i] = var[mat_cell_i];
    ++i;
  }
}

/* To access to MatCellVectorView on accelerator */
class MatCellContainer : 
  public Arcane::Accelerator::impl::MatCommandContainerBase
{
 public:
  using MatCellVectorView = Arcane::Materials::MatCellVectorView;
  using ComponentItemVectorView = Arcane::Materials::ComponentItemVectorView;
  using IMeshEnvironment = Arcane::Materials::IMeshEnvironment;
  using ComponentItemLocalId = Arcane::Materials::ComponentItemLocalId;
  using MatVarIndex = Arcane::Materials::MatVarIndex;

 public:
  explicit MatCellContainer(MatCellVectorView view) :
    Arcane::Accelerator::impl::MatCommandContainerBase(view)
  {
  }

  //! Accesseur pour le i-ème élément de la liste
  constexpr ARCCORE_HOST_DEVICE ComponentItemLocalId operator[](Arcane::Int32 i) const
  {   
    return { ComponentItemLocalId(m_matvar_indexes[i]) };
  }   
};

// Arcane var => raw_data, implem GPU API
void PhyVarType::
copyVarToRawData(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector)
{
  PROF_ACC_BEGIN(__FUNCTION__);

  auto command = ax::makeCommand(var.globalVariable().subDomain()->acceleratorMng()->defaultQueue());
  auto in_var = ax::viewIn(command, var);
  Arcane::Span<Arcane::Real> out_raw_data(raw_data.data(), raw_data.size());
  MatCellContainer mat_cell_cont(*mat_cell_vector);

  auto nelt = mat_cell_vector->nbItem();

  command << RUNCOMMAND_LOOP1(iter, nelt) 
  {
    auto i = iter()[0];

    auto mc = mat_cell_cont[i];
    out_raw_data[i] = in_var[mc];
  };

  PROF_ACC_END;
}

// raw_data => Arcane var
void PhyVarType::
copyRawDataToVar(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector) const
{
  for (size_t i(0) ; i < mat_cell_vector->size() ; ++i) {
    var[mat_cell_vector->operator[](i)] = raw_data[i];
  }
}

// raw_data => Arcane var, implem GPU API
void PhyVarType::
copyRawDataToVar(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector) const
{
  PROF_ACC_BEGIN(__FUNCTION__);

  auto command = ax::makeCommand(var.globalVariable().subDomain()->acceleratorMng()->defaultQueue());
  auto out_var = ax::viewOut(command, var);
  Arcane::Span<const Arcane::Real> in_raw_data(raw_data.data(), raw_data.size());
  MatCellContainer mat_cell_cont(*mat_cell_vector);

  auto nelt = mat_cell_vector->nbItem();

  command << RUNCOMMAND_LOOP1(iter, nelt) 
  {
    auto i = iter()[0];

    auto mc = mat_cell_cont[i];
    out_var[mc] = in_raw_data[i];
  };

  PROF_ACC_END;
}

} // Stdperfectgasacc1

