#ifndef IACCENV_H
#define IACCENV_H

#include <arcane/ItemTypes.h>
#include "accenv/AcceleratorUtils.h"

using namespace Arcane;

class IAccEnv {
 public:
  IAccEnv() {}
  virtual ~IAccEnv() {}

 public:
  virtual void initAcc() = 0;
  virtual ax::Runner& runner() = 0;
};

#endif
