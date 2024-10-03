// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#include "SprayTresFinService.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/** 
    Initialisation solver particules
*/
/*---------------------------------------------------------------------------*/
void SprayTresFinService::initSolverParticles() {

  IItemFamily* item_family = mesh()->findItemFamily (eItemKind::IK_Particle, "AllParticles");
  m_particles_family = item_family->toParticleFamily();
  activeParticlesGroup = item_family->findGroup("activeItem");
  
}

/*---------------------------------------------------------------------------*/
/** 
    correction de la vitesse du fluide 
    (prise en compte de la traînée due aux particules -> Rien à faire pour 
    ce service).
*/
/*---------------------------------------------------------------------------*/
void SprayTresFinService::correctFluidVelocity() {
}

/*---------------------------------------------------------------------------*/
/** 
    mise à jour de la vitesse des particules
*/
/*---------------------------------------------------------------------------*/
void SprayTresFinService::updateParticleVelocity() {

  // dt^n
  const Real dt ( 0.5 * ( m_global_old_deltat() + m_global_deltat() ) );

  ENUMERATE_PARTICLE (ipart, activeParticlesGroup) {
    Real3 up_chapo = m_particle_velocity[ipart] + dt*options()->getGravity();
    Particle particle = *ipart;
    Node node_proche = findNodeOfParticle(particle); // cherche à quel noeud appartient la particule (noeud le + proche)
    Real3 Dp = computeDp(particle, node_proche);
    Real3 one={1., 1., 1.};

    m_particle_velocity[ipart] = ( up_chapo + dt*Dp*m_velocity[node_proche] )/( one + dt*Dp );
  }
}

/*---------------------------------------------------------------------------*/
/** 
      mise à jour de la position des particules
*/
/*---------------------------------------------------------------------------*/
void SprayTresFinService::updateParticlePosition() {

  // dt^{n+1/2}
  const Real dt (m_global_deltat());

  ENUMERATE_PARTICLE (part_i, activeParticlesGroup) {
    m_particle_coord[part_i] += m_particle_velocity[part_i] * dt;
    // info() << "Particle " << part_i.localId() << " vel = " << m_particle_velocity[part_i];
    // info() << "Particle " << part_i.localId() << " coord = " << m_particle_coord[part_i];
  }
 
}

/*---------------------------------------------------------------------------*/
/** 
      Calcul du coefficient de traînée.

      Pour le moment, on se limite au cas Re_p -> infinity
*/
/*---------------------------------------------------------------------------*/

Real SprayTresFinService::computeCd(){

  // TODO if necessary: add fluid viscosity for case Re_p < 1000
  
  // Re_p = 2.*m_density*m_particle_radius/fluid_viscosity*abs(m_velocity_n[]-m_particle_velocity);
  // if (Re_p < 1000.)
  //   return 24./Re_p * (1. + 1./6.* pow(Re_p, 2./3.));
  // else
  return 0.424;
}

/*---------------------------------------------------------------------------*/
/** 
      Calcul du préfacteur $D_p$ de la force de traînée $F = - D_p (u-u_p)$ due à une particule.
*/
/*---------------------------------------------------------------------------*/

Real3 SprayTresFinService::computeDp(Particle particle, Node node_proche){

  Real3 Dp;
  TypesIpg::eDrag type = options()->drag.type();

  if (type == TypesIpg::LinearDrag){   // linear drag
    Real coef=options()->drag.getCoef();
    Dp = Real3(coef, coef, coef);
  }
  else {  // quadratic drag

    Real fluid_node_density = m_node_mass[node_proche]/m_node_volume[node_proche];
    Real Cd=computeCd();

    Dp.x = 3./8.*fluid_node_density/m_particle_density[particle]*Cd/m_particle_radius[particle]*abs(m_velocity_n[node_proche].x-m_particle_velocity[particle].x);
    Dp.y = 3./8.*fluid_node_density/m_particle_density[particle]*Cd/m_particle_radius[particle]*abs(m_velocity_n[node_proche].y-m_particle_velocity[particle].y);
    Dp.z = 3./8.*fluid_node_density/m_particle_density[particle]*Cd/m_particle_radius[particle]*abs(m_velocity_n[node_proche].z-m_particle_velocity[particle].z);
  }
  return Dp;
}

/*---------------------------------------------------------------------------*/
/*
  Associe une particule à un noeud (== à une cellule duale)
 */
/*---------------------------------------------------------------------------*/

Node SprayTresFinService::findNodeOfParticle(Particle particule){

  // TODO: déplacer cette fonction dans ipg_creation_of_particles/utilsIpg.h
  
  Cell cell = particule.cell();

  // chercher à quel noeud appartient la particule
  // pour le moment, on considère que la particule appartient au noeud le plus proche
  Real distance_min=1e30;
  Node node_proche;
  for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
    /* Real distance = (node_coord[inode] - particle_coord[particule]).normL2(); */
    Real distance = (m_node_coord[inode] - m_particle_coord[particule]).normL2();
    if (distance < distance_min){
      distance_min = distance;
      node_proche=*inode;
    }
  }
  return node_proche;
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_SPRAYTRESFIN(SprayTresFin, SprayTresFinService);
