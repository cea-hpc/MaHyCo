#ifndef __ACCENV_SINGLETON_I_ACC_ENV_H
#define __ACCENV_SINGLETON_I_ACC_ENV_H

#include "arcane/ServiceBuilder.h"
#include "accenv/IAccEnv.h"

class SingletonIAccEnv
{
 public:
  static IAccEnv* accEnv(ISubDomain* sd)
  {
    return ServiceBuilder<IAccEnv>(sd).getSingleton();
  }
};

#endif

