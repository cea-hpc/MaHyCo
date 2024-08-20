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

#ifndef CELLWATCHINGSERVICE_H
#define CELLWATCHINGSERVICE_H

#include "arcane/ArcaneTypes.h"
#include "arcane/utils/Real3.h"
#include "arcane/core/materials/MeshMaterialVariableRef.h"
#include "arcane/materials/IMeshMaterialMng.h"
#include "time_history/ITimeHistoryCategory.h"
#include "time_history/CellWatching_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Classe implémentant l'écriture des grandeurs sur un noeud au cours du temps
 */
class CellWatchingService 
: public ArcaneCellWatchingObject
{
public:
  /** Constructeur de la classe */
  CellWatchingService(const ServiceBuildInfo & sbi)
    : ArcaneCellWatchingObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~CellWatchingService() {};

public:
  /** 
   *  Initialise l'écriture des grandeurs sur le noeud choisi
   */
  void init();
  /** 
   *  Ecriture des sorties
   */
  void write(); 

 private:
  IMeshMaterialMng* mm;
};

#endif  // CELLWATCHINGSERVICE_H
