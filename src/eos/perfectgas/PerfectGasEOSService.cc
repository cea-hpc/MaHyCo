// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "PerfectGasEOSService.h"
#include "arcane/VariableView.h"
#include "arcane/ServiceBuilder.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/* Constructeur de la classe                                                 */
/*---------------------------------------------------------------------------*/
PerfectGasEOSService::PerfectGasEOSService(const ServiceBuildInfo & sbi)
  : ArcanePerfectGasEOSObject(sbi) {
  m_acc_env = ServiceBuilder<IAccEnv>(subDomain()).getSingleton();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::initEOS(IMeshEnvironment* env)
{
  Real adiabatic_cst = getAdiabaticCst(env);
  // Initialise l'énergie et la vitesse du son
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Real pressure = m_pressure[ev];
    Real density = m_density[ev];
    m_internal_energy[ev] = pressure / ((adiabatic_cst - 1.) * density);
    m_sound_speed[ev] = sqrt(adiabatic_cst * pressure / density);
  }
}

/*---------------------------------------------------------------------------*/
/* Formule PerfectGas pour calculer unitairement pression, vitesse du son et dp/de */
/* Cette formule est appelée dans applyEOS(...) et applyOneCellEOS(...)      */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE inline void compute_pressure_sndspd_PG(Real adiabatic_cst,
    Real density, Real internal_energy,
    Real& pressure, Real& sound_speed, Real& dpde) 
{
  pressure = (adiabatic_cst - 1.) * density * internal_energy;
  sound_speed = sqrt(adiabatic_cst * pressure / density);
  dpde = (adiabatic_cst - 1.) * density;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::applyEOS(IMeshEnvironment* env)
{
  Real adiabatic_cst = getAdiabaticCst(env);
  // Calcul de la pression et de la vitesse du son
#if 0
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;
    if (m_density[ev] == 0.) info() << ev.globalCell().localId() << " densité nulle";
    compute_pressure_sndspd_PG(adiabatic_cst,
        m_density[ev], m_internal_energy[ev],
        m_pressure[ev], m_sound_speed[ev], m_dpde[ev]);
  }
#elif 0
  Parallel::Foreach(env->envView(),[&](EnvItemVectorView view)
  {
    ENUMERATE_ENVCELL(ienvcell,view){
      compute_pressure_sndspd_PG(adiabatic_cst,
          m_density[ienvcell], m_internal_energy[ienvcell],
          m_pressure[ienvcell], m_sound_speed[ienvcell], m_dpde[ienvcell]);
    }
  });
#else
  // Mailles pures
  auto queue_pur = m_acc_env->newQueue();
  queue_pur.setAsync(true);
  {
    auto command = makeCommand(queue_pur);

    // Nombre de mailles pures de l'environnement
    Integer nb_pur = env->pureEnvItems().nbItem();

    // Pour les mailles pures, valueIndexes() est la liste des ids locaux des mailles
    Span<const Int32> in_cell_id(env->pureEnvItems().valueIndexes());

    auto in_density         = ax::viewIn(command, m_density.globalVariable());
    auto in_internal_energy = ax::viewIn(command, m_internal_energy.globalVariable());

    auto out_pressure    = ax::viewOut(command, m_pressure.globalVariable());
    auto out_sound_speed = ax::viewOut(command, m_sound_speed.globalVariable());
    auto out_dpde        = ax::viewOut(command, m_dpde.globalVariable());

    command << RUNCOMMAND_LOOP1(iter, nb_pur) {
      auto [ipur] = iter(); // ipur \in [0,nb_pur[
      CellLocalId cid(in_cell_id[ipur]); // accés indirect à la valeur de la maille

      Real pressure, sound_speed, dpde;

      compute_pressure_sndspd_PG(adiabatic_cst,
          in_density[cid], in_internal_energy[cid],
          pressure, sound_speed, dpde);

      out_pressure[cid] = pressure;
      out_sound_speed[cid] = sound_speed;
      out_dpde[cid] = dpde;

    }; // non-bloquant et asynchrone par rapport au CPU et autres queues
  }

  // Mailles mixtes
  auto queue_mix = m_acc_env->newQueue();
  queue_mix.setAsync(true);
  {
    auto command = makeCommand(queue_mix);

    // Nombre de mailles impures (mixtes) de l'environnement
    Integer nb_imp = env->impureEnvItems().nbItem();

    Span<const Real> in_density         (envView(m_density, env));
    Span<const Real> in_internal_energy (envView(m_internal_energy, env));

    Span<Real> out_pressure    (envView(m_pressure, env));
    Span<Real> out_sound_speed (envView(m_sound_speed, env));
    Span<Real> out_dpde        (envView(m_dpde, env));

    command << RUNCOMMAND_LOOP1(iter, nb_imp) {
      auto [imix] = iter(); // imix \in [0,nb_imp[

      compute_pressure_sndspd_PG(adiabatic_cst,
          in_density[imix], in_internal_energy[imix],
          out_pressure[imix], out_sound_speed[imix], out_dpde[imix]);

    }; // non-bloquant et asynchrone par rapport au CPU et autres queues
  }
  queue_pur.barrier();
  queue_mix.barrier();
#endif
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PerfectGasEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
  Real adiabatic_cst = getAdiabaticCst(env);
  // Calcul de la pression et de la vitesse du son
#if 0
    Real internal_energy = m_internal_energy[ev];
    Real density = m_density[ev];
    if (density == 0.) info() << ev.globalCell().localId() << " densité " << density;
    Real pressure = (adiabatic_cst - 1.) * density * internal_energy;
    m_pressure[ev] = pressure;
    m_sound_speed[ev] = sqrt(adiabatic_cst * pressure / density);
    m_dpde[ev] = (adiabatic_cst - 1.) * density;
#else
    if (m_density[ev] == 0.) info() << ev.globalCell().localId() << " densité nulle";
    compute_pressure_sndspd_PG(adiabatic_cst,
        m_density[ev], m_internal_energy[ev],
        m_pressure[ev], m_sound_speed[ev], m_dpde[ev]);
#endif
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real PerfectGasEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real PerfectGasEOSService::getTensionLimitCst(IMeshEnvironment* env) { return options()->limitTension();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_PERFECTGASEOS(PerfectGas, PerfectGasEOSService);
