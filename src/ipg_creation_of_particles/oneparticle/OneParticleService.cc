// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "OneParticleService.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


void OneParticleService::createParticles()
{
  // pour l'instant on ne crée qu'une seule particule.
  // La particule est injectée à la position spécifié à l'instant t^n tel que  t^n <= t_init_particule < t^n+dt
  if ((options()->initTime() >= m_global_time()-m_global_deltat()) && (options()->initTime() < m_global_time())){
    // note: m_global_time() contient le temps t^{n+1} !!!

    // il faut ajouter la particule au groupe activeItem ET à la famille AllParticles
    IItemFamily* item_family = mesh()->findItemFamily (eItemKind::IK_Particle, "AllParticles");
    IParticleFamily* m_particles_family = item_family->toParticleFamily();
    ParticleGroup activeParticlesGroup = item_family->findGroup("activeItem");

    UniqueArray<Integer> lids({0}); //local Id
    UniqueArray<Int64> uids({0}); //unique Id
    info() << "Création des particules de localId " << lids.view();
    m_particles_family->addParticles(uids.view(), lids.view());
    m_particles_family->endUpdate();
    activeParticlesGroup.addItems(lids.view());

    ENUMERATE_PARTICLE (part_i, activeParticlesGroup) {
    // ENUMERATE_PARTICLE (part_i, m_particles_family->allItems()) {
      m_particle_coord[part_i] = options()->getInitCoord();
      m_particle_velocity[part_i] = options()->getInitVelocity();
      m_particle_temperature[part_i] = options()->getInitTemperature();
      m_particle_weight[part_i] = options()->getInitWeight();
      m_particle_radius[part_i] = options()->getInitRadius();
    }
  }

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ONEPARTICLE(OneParticle, OneParticleService);
