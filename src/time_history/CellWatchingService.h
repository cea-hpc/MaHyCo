// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
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
