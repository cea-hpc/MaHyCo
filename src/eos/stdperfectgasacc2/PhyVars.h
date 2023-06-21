#ifndef _EOS_STDPERFECTGASACC2_PHY_VARS_H
#define _EOS_STDPERFECTGASACC2_PHY_VARS_H

namespace Stdperfectgasacc2 {

class PhyVars
{
 public:
  PhyVars();
  ~PhyVars();

  void addPhyVar(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector);
  
  void copyFromAllArcVars();

  Arcane::Span<Arcane::Real> rawData(Arcane::Int32 i);

  void copyIntoArcVar(Arcane::Int32 i) const;

  void clear();

 protected:
  void* m_impl=nullptr;
};

} // Stdperfectgasacc2

#endif

