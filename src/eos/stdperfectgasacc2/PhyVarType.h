#ifndef _EOS_STDPERFECTGASACC2_PHY_VAR_TYPE_H
#define _EOS_STDPERFECTGASACC2_PHY_VAR_TYPE_H

#include "eos/stdperfectgasacc2/EosTypes.h"
#include "arcane/materials/IMeshMaterial.h"
#include "arcane/materials/MeshMaterialVariableRef.h"
#include "arcane/core/materials/MaterialsCoreGlobal.h"

namespace Stdperfectgasacc2 {

typedef Arcane::Span<Arcane::Real> PhyVarRawData;

// Interface for Arcane variable to pack/unpack into raw data
class IPhyVarData
{
 public:
  IPhyVarData() {}
  virtual ~IPhyVarData() {}

  virtual Arcane::String name() const = 0;
  virtual Arcane::Int32 sizeInBytes() const = 0;
  virtual void setBuffer(Arcane::Span<Arcane::Byte> span_byte) = 0;
  virtual void copyVarToRawData() = 0;
  virtual void copyRawDataToVar() = 0;
  virtual PhyVarRawData rawData() = 0;
};

// Implementation for multi-mat Arcane variable to pack/unpack into raw data
class PhyMatVarData : public IPhyVarData
{
 public:
  PhyMatVarData(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector);
  virtual ~PhyMatVarData();

  Arcane::String name() const override;
  Arcane::Int32 sizeInBytes() const override;
  void setBuffer(Arcane::Span<Arcane::Byte> span_byte) override;
  void copyVarToRawData() override;
  void copyRawDataToVar() override;
  PhyVarRawData rawData() override;

 protected:
  Arcane::Materials::MaterialVariableCellReal m_var;
  const Arcane::Materials::MatCellVectorView* m_mat_cell_vector;
  PhyVarRawData m_raw_data;
};

struct PhyVarType
{
  Arcane::Ref<IPhyVarData> phy_var_data;

  PhyVarType();
  PhyVarType(Arcane::Ref<IPhyVarData> phy_var_data);
  explicit PhyVarType(const PhyVarType&);

  ~PhyVarType();

  PhyVarType& operator=(const PhyVarType&);
};

} // Stdperfectgasacc2

#endif

