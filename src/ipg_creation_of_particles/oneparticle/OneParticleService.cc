// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "OneParticleService.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void OneParticleService::createParticles(IParticleFamily* m_particles_family)
{
  
  // ENUMERATE_PARTICLE (part_i, m_particles_family->allItems()) {
  //   m_particle_coord[part_i] = Real3{0.5, 0.05, 0};
  //   Real init_velocity = 500.;
  //   Real angle = 2 * part_i.localId() * 3.141592 / 3;
  //   m_particle_velocity[part_i] = Real3(init_velocity * std::cos(angle), init_velocity * std::sin(angle), 0);
  // }

  // pour l'instant on ne crée qu'une seule particule.
  if ((options()->initTime() >= m_global_time()) && (options()->initTime() < m_global_time()+m_global_deltat())){
    UniqueArray<Integer> lids({0}); //local Id
    UniqueArray<Int64> uids({0}); //unique Id
    info() << "Création des particules de localId " << lids.view();
    m_particles_family->addParticles(uids.view(), lids.view());
    m_particles_family->endUpdate();

    ENUMERATE_PARTICLE (part_i, m_particles_family->allItems()) { // fixme: une seule particule, enumerate non necessaire ?
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
