// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "SprayFinService.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/** 
    Initialisation solver particules
*/
/*---------------------------------------------------------------------------*/
void SprayFinService::initSolverParticles() {

  IItemFamily* item_family = mesh()->findItemFamily (eItemKind::IK_Particle, "AllParticles");
  m_particles_family = item_family->toParticleFamily();
  activeParticlesGroup = item_family->findGroup("activeItem");

}

/*---------------------------------------------------------------------------*/
/** 
    correction de la vitesse du fluide (prise en compte de la traînée due aux particules).
*/
/*---------------------------------------------------------------------------*/
void SprayFinService::correctFluidVelocity() {

  const Real dt (m_global_deltat());

  // 1 - trouver les particules appartenant à chaque cellule *duale*

  // 2 - calcul de la force de trainée aux noeuds

  // 3 - MaJ de la vitesse du fluide aux noeuds
  
  // ENUMERATE_NODE ( inode, allNodes() ) {
  //   Node node = *inode;
  //   m_velocity[inode] += (dt / m_node_mass[inode]) * m_force[inode];
  // }
  // q. est-ce que je peux erase et utiliser m_force ?
  // ou est-ce que je dois creer une nouvelle variable ?

  
}

/*---------------------------------------------------------------------------*/
/** 
    mise à jour de la vitesse des particules
*/
/*---------------------------------------------------------------------------*/
void SprayFinService::updateParticleVelocity() {

}

/*---------------------------------------------------------------------------*/
/** 
      mise à jour de la position des particules
*/
/*---------------------------------------------------------------------------*/
void SprayFinService::updateParticlePosition() {

  // we need dt^{n+1/2}
  const Real dt ( 0.5 * ( m_global_old_deltat() + m_global_deltat() ) );
  
  ENUMERATE_PARTICLE (part_i, activeParticlesGroup) {
    m_particle_coord[part_i] += m_particle_velocity[part_i] * dt;
    info() << "Particle " << part_i.localId() << " vel = " << m_particle_velocity[part_i];
    info() << "Particle " << part_i.localId() << " coord = " << m_particle_coord[part_i];
  }
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_SPRAYFIN(SprayFin, SprayFinService);
