#ifndef _EOS_STDPERFECTGASACC1_EOS_TYPES_H
#define _EOS_STDPERFECTGASACC1_EOS_TYPES_H

#include <vector>
#include "arcane/materials/MatItem.h"
#include "eos/stdperfectgasacc1/Mallocator.h"

namespace Stdperfectgasacc1 {

typedef std::vector<Arcane::Materials::MatCell> StdMatCellVector;
typedef std::vector<Arcane::Real, Mallocator<Arcane::Real> > StdRealVector;

} // Stdperfectgasacc1

#endif

