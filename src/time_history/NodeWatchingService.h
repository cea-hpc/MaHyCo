// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef NODEWATCHINGSERVICE_H
#define NODEWATCHINGSERVICE_H

#include "time_history/ITimeHistoryCategory.h"
#include "time_history/NodeWatching_axl.h"

using namespace Arcane;

/**
 * Classe implémentant l'écriture des grandeurs sur un noeud au cours du temps
 */
class NodeWatchingService 
: public ArcaneNodeWatchingObject
{
public:
  /** Constructeur de la classe */
  NodeWatchingService(const ServiceBuildInfo & sbi)
    : ArcaneNodeWatchingObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~NodeWatchingService() {};

public:
  /** 
   *  Initialise l'écriture des grandeurs sur le noeud choisi
   */
  virtual void init();
  /** 
   *  Ecriture des sorties
   */
  virtual void write();  
};

#endif  // NODEWATCHINGSERVICE_H
