#ifndef __ACCENV_ACC_ENV_MODULE_H
#define __ACCENV_ACC_ENV_MODULE_H

#include "AccEnv_axl.h"
#include <arcane/ModuleBuildInfo.h>
class IAccEnv;

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Acc_env {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class AccEnvModule : public ArcaneAccEnvObject
{
 public:
  AccEnvModule(const ModuleBuildInfo& mbi);
  ~AccEnvModule();

  //! Initialise les environnements pour les accélérateurs
  void accBuild() override;

  //! Initialise des informations liées au maillage pour les accélérateurs
  void initMesh() override;

  //! Initialise des informations liées au multi-environnement pour les accélérateurs
  void initMultiEnv() override;

  //! Amorce l'instrumentation pour le profiling
  void startInstrument() override;

  //! Termine l'instrumentation pour le profiling
  void stopInstrument() override;

 protected:
  //! Pour l'utilisation des accélérateurs
  IAccEnv* m_acc_env=nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}  // namespace Acc_env

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif

