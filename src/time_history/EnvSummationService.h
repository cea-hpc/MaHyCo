// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef ENVSUMMATION_H
#define ENVSUMMATION_H

#include "arcane/ArcaneTypes.h"
#include "arcane/utils/Real3.h"
#include "arcane/core/materials/MeshMaterialVariableRef.h"
#include "time_history/ITimeHistoryCategory.h"
#include "time_history/EnvSummation_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Classe implémentant l'écriture des grandeurs sur un noeud au cours du temps
 */
class EnvSummationService 
: public ArcaneEnvSummationObject
{
public:
  /** Constructeur de la classe */
  EnvSummationService(const ServiceBuildInfo & sbi)
    : ArcaneEnvSummationObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~EnvSummationService() {};

public:
  /** 
   *  Initialise l'écriture des grandeurs sur les milieux choisis
   */
  virtual void init();
  /** 
   *  Ecriture des sorties
   */
  virtual void write(); 
 
 private:
  UniqueArray<IMeshEnvironment*> m_environments;
};

#endif  // ENVSUMMATION_H
