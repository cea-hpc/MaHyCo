#include "AccEnvDefaultService.h"
#include <arcane/AcceleratorRuntimeInitialisationInfo.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void AccEnvDefaultService::
initAcc()
{
  info() << "Using accelerator";
  IApplication* app = subDomain()->application();
  initializeRunner(m_runner,traceMng(),app->acceleratorRuntimeInitialisationInfo());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ACCENVDEFAULT(AccEnvDefault, AccEnvDefaultService);

