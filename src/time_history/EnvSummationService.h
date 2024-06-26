// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef ENVSUMMATION_H
#define ENVSUMMATION_H

#include <memory>
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
  Real _computeAvgVarForEnv(
    IParallelMng* pm, IMeshEnvironment* env, const Integer nb_envcells, MaterialVariableCellReal& variable);

  Real _computeExtensiveVarForEnv(
    IParallelMng* pm, IMeshEnvironment* env, MaterialVariableCellReal& variable, MaterialVariableCellReal& globalization_var);
 
 private:
  UniqueArray<IMeshEnvironment*> m_environments;  /// liste des milieux à sortir

  UniqueArray<MaterialVariableCellReal*> m_avg_var_list; /// += var / nb_cells
  UniqueArray<MaterialVariableCellReal*> m_massic_var_list;  /// += mass * var / mass_tot
  UniqueArray<MaterialVariableCellReal*> m_volumic_var_list;  /// += vol * var / vol_tot
};

#endif  // ENVSUMMATION_H
