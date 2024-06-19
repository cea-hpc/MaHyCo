// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "TimeHistoryModule.h"
#include "ITimeHistoryCategory.h"


using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/**
 * \brief Constructeur du module
 */
/*---------------------------------------------------------------------------*/
TimeHistoryModule:: TimeHistoryModule(const ModuleBuildInfo& mbi) : 
ArcaneTimeHistoryObject(mbi)
{
}


/*---------------------------------------------------------------------------*/
/**
 * \brief Initialisation du module
 */
/*---------------------------------------------------------------------------*/
void TimeHistoryModule::
init()
{
  for (ITimeHistoryCategory* bilan : options()->bilan()) {
    bilan->init();
  }
}

/*---------------------------------------------------------------------------*/
/**
 * \brief Ecriture des sorties bilan
 */
/*---------------------------------------------------------------------------*/
void TimeHistoryModule::
write()
{
  for (ITimeHistoryCategory* bilan : options()->bilan()) {
    bilan->write();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_TIMEHISTORY(TimeHistoryModule);  
