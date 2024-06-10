// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "IpgModule.h"

#include "arcane/ITimeLoopMng.h"
#include "arcane/IMesh.h"
#include "arcane/IItemFamily.h"
#include "arcane/IParticleFamily.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/**
 * \brief Constructeur du module
 */
/*---------------------------------------------------------------------------*/
IpgModule:: IpgModule(const ModuleBuildInfo& mbi) : 
  ArcaneIpgObject(mbi)
{
  auto* item_family = mesh()->createItemFamily (eItemKind::IK_Particle, "AllParticles");

  // creation de groupes pour les particules :
  // celles qui sont actives suivront les equations des particules
  // on aura besoin ailleurs de créer d'autres groupes qui ne suivront pas necessairement les equations.
  activeParticlesGroup = item_family->createGroup("activeItem");

  m_particles_family = item_family->toParticleFamily();

}

/*---------------------------------------------------------------------------*/
/**
 * \brief Création des particules
 */
/*---------------------------------------------------------------------------*/
void IpgModule::
injectParticles()
{
  options()->getParticleInjectorType()->createParticles();
}

/*---------------------------------------------------------------------------*/
/**
 * \brief Initialisation éventuellement nécessaires dans start-init pour créer 
 * les particules.
 */
/*---------------------------------------------------------------------------*/
void IpgModule::
initInjectParticles()
{
  options()->getParticleInjectorType()->initParticles();
}

/*---------------------------------------------------------------------------*/
/**
 * \brief Initialisation des sorties
 */
/*---------------------------------------------------------------------------*/
void IpgModule::
initParticleOutput()
{
  options()->getIpgOutput()->initOutput();
}

/*---------------------------------------------------------------------------*/
/**
 * \brief Initialisation du solver pour les particules dans start-init
 */
/*---------------------------------------------------------------------------*/
void IpgModule::
initSolverParticles()
{
  options()->getSprayType()->initSolverParticles();
}

/*---------------------------------------------------------------------------*/
/**
 * \brief Mise à jour de la position des particules
 */
/*---------------------------------------------------------------------------*/
void IpgModule::
updateParticlePosition()
{
  options()->getSprayType()->updateParticlePosition();
}

/*---------------------------------------------------------------------------*/
/**
 * \brief Mise à jour de la vitesse des particules
 */
/*---------------------------------------------------------------------------*/
void IpgModule::
updateParticleVelocity()
{
  options()->getSprayType()->updateParticleVelocity();
}

/*---------------------------------------------------------------------------*/
/**
 * \brief Correction de la vitesse du fluide en présence de particules (traînée)
 */
/*---------------------------------------------------------------------------*/
void IpgModule::
correctFluidVelocity()
{
  options()->getSprayType()->correctFluidVelocity();
}

/*---------------------------------------------------------------------------*/
/**
 * \brief Ecriture des sorties
 */
/*---------------------------------------------------------------------------*/
void IpgModule::
writeParticleOutput()
{
  options()->getIpgOutput()->writeOutput(m_particles_family->allItems());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_IPG(IpgModule);
