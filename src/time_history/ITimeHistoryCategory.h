// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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
