#ifndef _EOS_STDPERFECTGASACC1_PHY_VAR_TYPE_H
#define _EOS_STDPERFECTGASACC1_PHY_VAR_TYPE_H

#include "eos/stdperfectgasacc1/EosTypes.h"
#include "arcane/materials/IMeshMaterial.h"
#include "arcane/materials/MeshMaterialVariableRef.h"

namespace Stdperfectgasacc1 {

#define PHY_VAR_RAW_DATA_IS_STD

#ifdef PHY_VAR_RAW_DATA_IS_STD
typedef StdRealVector PhyVarRawData;
#else
typedef Arcane::UniqueArray<Arcane::Real> PhyVarRawData;
#endif

struct PhyVarType
{
  Arcane::String var_name;
  PhyVarRawData raw_data;

  PhyVarType();
  PhyVarType(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector);
  PhyVarType(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector);
  explicit PhyVarType(const PhyVarType&);

  ~PhyVarType();

  PhyVarType& operator=(const PhyVarType&);

  // Arcane var => raw_data
  void copyVarToRawData(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector);

  // Arcane var => raw_data, implem GPU API
  void copyVarToRawData(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector);


  // raw_data => Arcane var
  void copyRawDataToVar(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector) const;

  // raw_data => Arcane var, implem GPU API
  void copyRawDataToVar(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector) const;
};

} // Stdperfectgasacc1

#endif

