// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "OneParticleService.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void OneParticleService::createParticles(IParticleFamily* m_particles_family)
{
  
  UniqueArray<Integer> lids({0,1,2}); //local Id
  UniqueArray<Int64> uids({0,1,2}); //unique Id
  info() << "Création des particules de localId " << lids.view();
  m_particles_family->addParticles(uids.view(), lids.view());
  m_particles_family->endUpdate();

  ENUMERATE_PARTICLE (part_i, m_particles_family->allItems()) {
    m_particle_coord[part_i] = Real3{0.5, 0.05, 0};
    Real init_velocity = 500.;
    Real angle = 2 * part_i.localId() * 3.141592 / 3;
    m_particle_velocity[part_i] = Real3(init_velocity * std::cos(angle), init_velocity * std::sin(angle), 0);
  }

  // // pour l'instant on ne crée qu'une seule particule.
  // // plus tard il serait bon que l'on puisse en créer plusieurs comme ça
  // // peut-etre au autre service ManyParticlesService
  // UniqueArray<Integer> lids({0}); //local Id
  // UniqueArray<Int64> uids({0}); //unique Id
  // info() << "Création des particules de localId " << lids.view();
  // m_particles_family->addParticles(uids.view(), lids.view());       // m_particle_family created in IpgModule. Do I need it in class definition OneParticuleService ?
  // m_particles_family->endUpdate();

  // ENUMERATE_PARTICLE (part_i, m_particles_family->allItems()) {
  //   m_particle_coord[part_i] = Real3{0.5, 0.05, 0};
  //   m_particle_velocity[part_i] = options()->getInitVelocity();
  //   // m_particle_velocity[part_i] = options()->initVelocity(); // both work, I don't know which is best as programming practice
  // }
  
  // Real T = options()->initTemperature();
  // Real3 u = options()->initVelocity();
  // Real3 x = options()->initCoord();
  // Real w = options()->initWeight();
  // Real r = options()->initRadius();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ONEPARTICLE(OneParticle, OneParticleService);
