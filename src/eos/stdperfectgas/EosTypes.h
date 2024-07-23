// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#ifndef _EOS_STDPERFECTGAS_EOS_TYPES_H
#define _EOS_STDPERFECTGAS_EOS_TYPES_H

#include <vector>
#include "arcane/materials/MatItem.h"

namespace Stdperfectgas {

typedef std::vector<Arcane::Materials::MatCell> StdMatCellVector;
typedef std::vector<Arcane::Real> StdRealVector;

} // Stdperfectgas

#endif

