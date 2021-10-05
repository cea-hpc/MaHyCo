#ifndef CARTESIAN_INTERFACE_CARTESIAN_NAMESPACE_H
#define CARTESIAN_INTERFACE_CARTESIAN_NAMESPACE_H

#include "cartesian/interface/CartesianDefines.h"

#ifdef USE_CARTESIAN_IMPL
using namespace Cartesian;
namespace CartesianInterface = Cartesian;
#else
using namespace Arcane;
namespace CartesianInterface = Arcane;
#endif

#endif

