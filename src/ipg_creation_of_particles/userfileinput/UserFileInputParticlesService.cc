// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "UserFileInputParticlesService.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* particule creation in compute-loop */
/*---------------------------------------------------------------------------*/

void UserFileInputParticlesService::createParticles(IParticleFamily* m_particles_family)
{
  // we fill in the list of particule injected in the simulation



  // we fill out the list of particules to be injected

}

/*---------------------------------------------------------------------------*/
/* initialisation in start-init: read the user file and store particule data */
/*---------------------------------------------------------------------------*/

void UserFileInputParticlesService::initParticles()
{
  // we fill in the familly of particules which will later be introduce in the simulation

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_USERFILEINPUTPARTICLES(UserFileInputParticles, UserFileInputParticlesService);
