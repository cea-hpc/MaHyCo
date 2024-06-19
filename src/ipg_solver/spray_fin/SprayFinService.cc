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

  const Real dt (m_global_deltat());  //  dt^n

  // ### Ce qu'il faudrait faire idéalement :
  // 1 - trouver les particules appartenant à chaque *cellule duale*
  // 2 - calcul de la force de trainée aux noeuds
  // 3 - MaJ de la vitesse du fluide aux noeuds

  // ### Ce que l'on fait pour le moment :
  // 1 - calcul du noeud le plus proche de la particule (ce n'est pas forcément juste suivant la définition de la cellule duale)
  // 2 - calcul de la force de trainée aux noeuds
  // 3 - MaJ de la vitesse du fluide aux noeuds

  // 0 - Réinitialisation des vecteurs des forces et du denominateur.
  m_force.fill ( Real3::zero() );
  m_denom.fill ( Real3::zero() );
  ENUMERATE_NODE(inode, allNodes()){
    m_denom[inode] = 1.;
  }

  ENUMERATE_PARTICLE (ipart, activeParticlesGroup) {
    Particle particle = *ipart;
    Cell cell = particle.cell();

    // 1 - chercher à quel noeud appartient la particule
    //     pour le moment, on considère que la particule appartient au noeud le plus proche
    Real distance_min=1e30;
    Node node_proche;
    for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
      Real distance = (m_node_coord[inode] - m_particle_coord[ipart]).normL2();
      if (distance < distance_min){
        distance_min = distance;
        node_proche=*inode;
      }
    }

    // 2 - calcul de la force de traînée aux noeuds
    //     (ajout de la contribution de la particule à son noeud)
    Real Np=m_particle_weight[ipart];
    Real mp=4./3.*Pi*pow(m_particle_radius[ipart], 3.)*m_particle_density[ipart];
    Real3 up_chapo = m_particle_velocity[ipart] + dt*options()->getGravity();  // TODO: eviter de répéter l'option gravity qui est deja dans mahyco.axl
    Real fluid_node_density = m_node_mass[node_proche]/m_node_volume[node_proche];
    Real Cd=computeCd();

    Real3 Dp;
    Dp.x = 3./8.*fluid_node_density/m_particle_density[ipart]*Cd/m_particle_radius[ipart]*abs(m_velocity_n[node_proche].x-m_particle_velocity[ipart].x);
    Dp.y = 3./8.*fluid_node_density/m_particle_density[ipart]*Cd/m_particle_radius[ipart]*abs(m_velocity_n[node_proche].y-m_particle_velocity[ipart].y);
    Dp.z = 3./8.*fluid_node_density/m_particle_density[ipart]*Cd/m_particle_radius[ipart]*abs(m_velocity_n[node_proche].z-m_particle_velocity[ipart].z);

    Real3 one={1., 1., 1.};
    m_denom[node_proche] += dt/m_node_mass[node_proche]*Np*(mp*Dp)/(one + dt*Dp);
    m_force[node_proche] += dt/m_node_mass[node_proche]*Np*(mp*Dp)/(one + dt*Dp)*up_chapo;
  }

  // 3 - MaJ de la vitesse du fluide
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] += m_force[inode];
    m_velocity[inode] /= m_denom[inode];
  }
}

/*---------------------------------------------------------------------------*/
/** 
    mise à jour de la vitesse des particules
*/
/*---------------------------------------------------------------------------*/
void SprayFinService::updateParticleVelocity() {

  const Real dt (m_global_deltat());  //  dt^n

  ENUMERATE_PARTICLE (ipart, activeParticlesGroup) {

    Real3 up_chapo = m_particle_velocity[ipart] + dt*options()->getGravity();

    Particle particle = *ipart;
    Cell cell = particle.cell();
    
    // chercher à quel noeud appartient la particule
    // pour le moment, on considère que la particule appartient au noeud le plus proche
    Real distance_min=1e30;
    Node node_proche;
    for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
      Real distance = (m_node_coord[inode] - m_particle_coord[ipart]).normL2();
      if (distance < distance_min){
        distance_min = distance;
        node_proche=*inode;
      }
    }

    Real fluid_node_density = m_node_mass[node_proche]/m_node_volume[node_proche];
    Real Cd=computeCd();
    Real3 Dp;
    Dp.x = 3./8.*fluid_node_density/m_particle_density[ipart]*Cd/m_particle_radius[ipart]*abs(m_velocity_n[node_proche].x-m_particle_velocity[ipart].x);
    Dp.y = 3./8.*fluid_node_density/m_particle_density[ipart]*Cd/m_particle_radius[ipart]*abs(m_velocity_n[node_proche].y-m_particle_velocity[ipart].y);
    Dp.z = 3./8.*fluid_node_density/m_particle_density[ipart]*Cd/m_particle_radius[ipart]*abs(m_velocity_n[node_proche].z-m_particle_velocity[ipart].z);

    Real3 one={1., 1., 1.};
    m_particle_velocity[ipart] = ( up_chapo + dt*Dp*m_velocity[node_proche] )/( one + dt*Dp );
  }
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

Real SprayFinService::computeCd(){

  // TODO if necessary: add fluid viscosity for case Re_p < 1000
  
  // Re_p = 2.*m_density*m_particle_radius/fluid_viscosity*abs(m_velocity_n[]-m_particle_velocity);
  // if (Re_p < 1000.)
  //   return 24./Re_p * (1. + 1./6.* pow(Re_p, 2./3.));
  // else
  return 0.424;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_SPRAYFIN(SprayFin, SprayFinService);
