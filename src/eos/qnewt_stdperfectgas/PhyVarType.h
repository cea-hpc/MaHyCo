// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#ifndef _EOS_QNEWT_STDPERFECTGAS_PHY_VAR_TYPE_H
#define _EOS_QNEWT_STDPERFECTGAS_PHY_VAR_TYPE_H

#include "eos/qnewt_stdperfectgas/EosTypes.h"
#include "arcane/materials/IMeshMaterial.h"
#include "arcane/materials/MeshMaterialVariableRef.h"

namespace Qnewt_stdperfectgas {

typedef StdRealVector PhyVarRawData;

struct PhyVarType
{
  Arcane::String var_name;
  PhyVarRawData raw_data;

  PhyVarType();
  PhyVarType(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector);
  explicit PhyVarType(const PhyVarType&);

  ~PhyVarType();

  PhyVarType& operator=(const PhyVarType&);

  // Arcane var => raw_data
  void copyVarToRawData(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector);

  // raw_data => Arcane var
  void copyRawDataToVar(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector) const;
};

} // Qnewt_stdperfectgas

#endif

