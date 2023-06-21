#ifndef _EOS_STDPERFECTGASACC2_PHY_VARS_H
#define _EOS_STDPERFECTGASACC2_PHY_VARS_H

#include "accenv/AcceleratorUtils.h"

namespace Stdperfectgasacc2 {

class PhyVars
{
 public:
  PhyVars();
  ~PhyVars();

  void addPhyVar(Arcane::Materials::MaterialVariableCellReal var, const Arcane::Materials::MatCellVectorView* mat_cell_vector);
  
  void asyncCopyFromAllArcVars(Arcane::Ref<ax::RunQueue> async_queue);

  Arcane::Span<Arcane::Real> rawData(Arcane::Int32 i);

  void asyncCopyIntoArcVar(Arcane::Ref<ax::RunQueue> async_queue, Arcane::Int32 i) const;

  void clear();

 protected:
  void* m_impl=nullptr;
};

} // Stdperfectgasacc2

#endif

