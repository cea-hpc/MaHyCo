// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

/*
Copyright 2000-2024 CEA

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

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
