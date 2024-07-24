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

#ifndef TIMEHISTORYMODULE_H
#define TIMEHISTORYMODULE_H

#include "ITimeHistoryCategory.h"
#include "time_history/TimeHistory_axl.h"

using namespace Arcane;

/**
 * Représente un module pour l'écriture des sorties TimeHistory
 * c'est à dire un bilan temporel des quantités thermo
 */


class TimeHistoryModule
: public ArcaneTimeHistoryObject
{
 public:
  /** Constructeur de la classe */
  TimeHistoryModule(const ModuleBuildInfo& mbi);
  
  /** Destructeur de la classe */
  ~TimeHistoryModule() {}

  /** Initialisation du module */
  virtual void init();
  
  /** Ecriture des sorties time history */
  virtual void write();


};

#endif  // TIMEHISTORYMODULE_H