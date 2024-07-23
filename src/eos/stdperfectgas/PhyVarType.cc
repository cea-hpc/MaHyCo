// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#include "eos/stdperfectgas/PhyVarType.h"

namespace Stdperfectgas {

PhyVarType::
PhyVarType()
{
}

PhyVarType::
PhyVarType(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector) :
  var_name (var.name()),
  raw_data (mat_cell_vector->size())
{
  // Arcane var => raw_data
  copyVarToRawData(var, mat_cell_vector);
}

PhyVarType::
PhyVarType(const PhyVarType& other) :
  var_name (other.var_name),
  raw_data (other.raw_data)
{
}

PhyVarType::
~PhyVarType()
{
}

PhyVarType& PhyVarType::
operator=(const PhyVarType& other)
{
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

// raw_data => Arcane var
void PhyVarType::
copyRawDataToVar(Arcane::Materials::MaterialVariableCellReal var, const StdMatCellVector* mat_cell_vector) const
{
  for (size_t i(0) ; i < mat_cell_vector->size() ; ++i) {
    var[mat_cell_vector->operator[](i)] = raw_data[i];
  }
}

} // Stdperfectgas

