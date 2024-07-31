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

#include "EnvSummationService.h"
#include "arcane/IParallelMng.h"
#include "arcane/core/ITimeHistoryMng.h"
#include "arcane/core/GlobalTimeHistoryAdder.h"
#include "arcane/core/materials/EnvItemVector.h"


using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Initialisation de la surveillance temporelle du milieu                    */
/*---------------------------------------------------------------------------*/
void EnvSummationService::init()
{
  auto* mm = IMeshMaterialMng::getReference ( mesh() );

  // Récupération de la liste des milieux
  UniqueArray<String> env_list(options()->environment.view());
  // Si aucun milieu renseigné, on écrit tous les milieux
  if (env_list.size() == 0) {
    ENUMERATE_ENV(env, mm) {
      m_environments.add(*env);
    }
  } else {
  // Sinon, on écrit les milieux renseignés dans le jdd uniquement
    ENUMERATE_ENV(env, mm) {
      if (std::find_if(env_list.begin(), env_list.end(), 
          [&env](const String wanted_name){return wanted_name == (*env)->name();}) != env_list.end())
      m_environments.add(*env);
    }
  }

  // Récupération de la liste des variables
  UniqueArray<String> jdd_var_list(options()->variable.view());

  // Par défaut, toutes les variables sont moyennées sans pondération
  // (pour aller plus vite et brancher la non-régression)
  // Il faudra à terme trouver un mécanisme pour choisir la bonne moyenne
  // pour écrire les bilans
  m_avg_var_list.add(&m_cell_mass); // sortie systématiquement
  for (const String& var_name : jdd_var_list) {
    if (var_name == "Density") {
      m_volumic_var_list.add(&m_density);
    } else if (var_name == "Pressure") {
      m_avg_var_list.add(&m_pressure);
    } else if (var_name == "InternalEnergy") {
      m_massic_var_list.add(&m_internal_energy);
    } else {
      info() << "La variable " << var_name << " n'est pour l'instant pas prévue dans les sorties bilan par milieu";
      // Quand la classe aura un peu évolué, ce sera une erreur. Pour l'instant, juste un print dans le listing
    }
  }

}

/*---------------------------------------------------------------------------*/
/* Ecriture des sorties temporelles pour la surveillance du milieu           */
/*---------------------------------------------------------------------------*/
void EnvSummationService::write()
{
  Integer period = options()->periode;
  if (m_global_iteration()%period !=0) return;

  IParallelMng* pm = subDomain()->parallelMng();
  Integer my_rank = pm->commRank();
  String filename;
  Real global_var;

  // Création de l'objet permettant d'ajouter des valeurs dans un historique des valeurs.
  GlobalTimeHistoryAdder global_adder(subDomain()->timeHistoryMng());

  for (IMeshEnvironment* env : m_environments) {

    // Calcul du nombre de mailles dans le milieu
    EnvCellVector own_envcells(allCells().own(), env);
    Integer nb_envcells = own_envcells.view().nbItem();
    nb_envcells = pm->reduce(IParallelMng::eReduceType::ReduceSum, nb_envcells);

    if (nb_envcells == 0) continue; // milieu vide dans tous les sous-domaines

    for (auto var : m_avg_var_list) {
      global_var = _computeAvgVarForEnv(pm, env, nb_envcells, *var);
      filename = var->globalVariable().name() + String("_") + env->name();
      global_adder.addValue(TimeHistoryAddValueArg(filename), global_var);
    }

    for (auto var : m_massic_var_list) {
      global_var = _computeExtensiveVarForEnv(pm, env, *var, m_cell_mass);
      filename = var->globalVariable().name() + String("_") + env->name();
      global_adder.addValue(TimeHistoryAddValueArg(filename), global_var);
    }

    for (auto var : m_volumic_var_list) {
      global_var = _computeExtensiveVarForEnv(pm, env, *var, m_cell_volume);
      filename = var->globalVariable().name() + String("_") + env->name();
      global_adder.addValue(TimeHistoryAddValueArg(filename), global_var);
    }

  }
}

/*---------------------------------------------------------------------------*/
/* Ecriture des sorties temporelles pour la surveillance du milieu           */
/*---------------------------------------------------------------------------*/
Real EnvSummationService::_computeAvgVarForEnv(
  IParallelMng* pm, IMeshEnvironment* env, const Integer nb_envcells, MaterialVariableCellReal& variable)
{
  // Calcul de la moyenne des pressions du sous-domaine.
  Real accumulated_var = 0;
  ENUMERATE_ENVCELL (envcell, env) {
    if ((*envcell).globalCell().isOwn()) {
      accumulated_var += variable[envcell];
      // todo : trier le tableau pour avoir une somme répétable
    }
  }

  // Calcul de la pression moyenne globale.
  accumulated_var = pm->reduce(IParallelMng::eReduceType::ReduceSum, accumulated_var); 
  accumulated_var /= nb_envcells;

  return accumulated_var;
}

/*---------------------------------------------------------------------------*/
/** Ecriture des sorties temporelles pour la surveillance du milieu  
  * @param pm : gestionnaire du parallélisme Arcane
  * @param variable : grandeur à écrire
  * @param globalization_var : grandeur servant à la globalisation (masse ou volume)
  */
/*---------------------------------------------------------------------------*/
Real EnvSummationService::_computeExtensiveVarForEnv(
  IParallelMng* pm, IMeshEnvironment* env, MaterialVariableCellReal& variable, MaterialVariableCellReal& globalization_var)
{
  // Calcul de la moyenne des pressions du sous-domaine.
  Real accumulated_var = 0;
  Real accumulated_globalization_var = 0;
  ENUMERATE_ENVCELL (envcell, env) {
    if ((*envcell).globalCell().isOwn()) {
      accumulated_var += variable[envcell] * globalization_var[envcell];
      accumulated_globalization_var += globalization_var[envcell];
      // todo : trier le tableau pour avoir une somme répétable
    }
  }

  // Calcul de la pression moyenne globale.
  accumulated_var = pm->reduce(IParallelMng::eReduceType::ReduceSum, accumulated_var); 
  accumulated_globalization_var = pm->reduce(IParallelMng::eReduceType::ReduceSum, accumulated_globalization_var); 
  accumulated_var /= accumulated_globalization_var;

  return accumulated_var;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ENVSUMMATION(EnvSummation, EnvSummationService);
