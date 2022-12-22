#ifndef __ACCENV_SINGLETON_I_ACC_ENV_H
#define __ACCENV_SINGLETON_I_ACC_ENV_H

#include "arcane/ServiceBuilder.h"
#include "accenv/IAccEnv.h"

class SingletonIAccEnv
{
 public:

  //! Retourne l'instance du singleton IAccEnv (qui doit forcement exister)
  static IAccEnv* accEnv(ISubDomain* sd)
  {
    return ServiceBuilder<IAccEnv>(sd).getSingleton();
  }

  //! Retourne l'instance du singleton IAccEnv si celui-ci existe, nullptr sinon
  static IAccEnv* accEnvIfExists(ISubDomain* sd)
  {
    return ServiceBuilder<IAccEnv>(sd).getSingleton(SB_AllowNull);
  }
};

#endif

