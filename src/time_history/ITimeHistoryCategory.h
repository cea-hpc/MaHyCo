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

#ifndef ITIMEHISTORYCATEGORY_H
#define ITIMEHISTORYCATEGORY_H

#include "arcane/ArcaneTypes.h"

using namespace Arcane;

/**
 * Interface de l'écriture des sorties bilan temporel
 */
class ITimeHistoryCategory
{
public:
  /** Constructeur de la classe */
  ITimeHistoryCategory() {};
  /** Destructeur de la classe */
  virtual ~ITimeHistoryCategory() {};
  
public:
  /** 
   *  Initialise le service d'écriture des sorties bilan
   */
  virtual void init() = 0;
  /** 
   *  Ecrit les sorties bilan
   */
  virtual void write() = 0;

  
};

#endif  // ITIMEHISTORYCATEGORY_H
