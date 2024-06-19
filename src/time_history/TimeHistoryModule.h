// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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