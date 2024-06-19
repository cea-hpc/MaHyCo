// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "EnvSummationService.h"
#include "arcane/IParallelMng.h"
#include "arcane/core/ITimeHistoryMng.h"
#include "arcane/core/GlobalTimeHistoryAdder.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Initialisation de la surveillance temporelle du milieu                    */
/*---------------------------------------------------------------------------*/
void EnvSummationService::init()
{
  auto* mm = IMeshMaterialMng::getReference ( mesh() );

  UniqueArray<String> env_list = options()->environment.view();

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

  // Création de l'objet permettant d'ajouter des valeurs dans un historique des valeurs.
  GlobalTimeHistoryAdder global_adder(subDomain()->timeHistoryMng());

  for (IMeshEnvironment* env : m_environments) {
    // Calcul de la moyenne des pressions du sous-domaine.
    Real var_subdomain = 0;
    Integer nb_envcells = 0;
    ENUMERATE_ENVCELL (envcell, env) {
      if ((*envcell).globalCell().isOwn()) {
        var_subdomain += m_pressure[envcell];
        // todo : trier le tableau pour avoir une somme répétable
        nb_envcells += 1;
      }
    }
 
    // Calcul de la pression moyenne globale.
    var_subdomain = pm->reduce(IParallelMng::eReduceType::ReduceSum, var_subdomain);
    nb_envcells = pm->reduce(IParallelMng::eReduceType::ReduceSum, nb_envcells);
    if (nb_envcells == 0) continue; // milieu vide dans tous les sous-domaines
    var_subdomain /= nb_envcells;
 
    // Ajout de la valeur avg_pressure_global dans l'historique "avg_pressure". Historique global.
    String filename = String("pressure_") + env->name();
    global_adder.addValue(TimeHistoryAddValueArg(filename), var_subdomain);
  }
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ENVSUMMATION(EnvSummation, EnvSummationService);
