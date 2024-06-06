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

    ENUMERATE_PARTICLE (part_i, activeParticlesGroup) {  // TODO: on n'a qu'une seule particule, enumerate inutile ici ?
      m_particle_coord[part_i] = options()->getInitCoord();
      m_particle_velocity[part_i] = options()->getInitVelocity();
      m_particle_temperature[part_i] = options()->getInitTemperature();
      m_particle_weight[part_i] = options()->getInitWeight();
      m_particle_radius[part_i] = options()->getInitRadius();
    }


    // il faut affecter la particule à la cellule à laquelle elle appartient
    assignParticleToCell();
    
  }
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void OneParticleService::assignParticleToCell()
{
  IItemFamily* item_family = mesh()->findItemFamily (eItemKind::IK_Particle, "AllParticles");
  IParticleFamily* m_particles_family = item_family->toParticleFamily();
  ParticleGroup activeParticlesGroup = item_family->findGroup("activeItem");

  // on fait un truc lourdingue pour le moment en parcourant toutes les particules et toutes les cellules
  // mais comme il n'y a qu'une particule, c'est pas grave
  ENUMERATE_PARTICLE (part_i, activeParticlesGroup) {
    Particle particule = *part_i;
    ENUMERATE_CELL ( icell, allCells() ) {
      Cell cell = * icell;
      Real3 particule_coord = m_particle_coord[part_i] ;

      UniqueArray<Real3> nodes_coord;
      for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
        Real3 node_coord = m_node_coord[inode] ;
        nodes_coord.add(node_coord);
      }
        
      // if (isCoordInCell(particule_coord, cell))
      if (isCoordInCell(particule_coord, nodes_coord))
        m_particles_family->setParticleCell(particule, cell);
    }
  }

  // on vérifie que toutes les (la) particules sont dans une cellule.
  // les particules en dehors du domaine maillé n'en ont pas.
  ENUMERATE_PARTICLE (part_i, activeParticlesGroup) {
    if (!(part_i->hasCell()))
      info() << "WARNING: Particle " << part_i.localId() << " located in " << m_particle_coord[part_i] << " has no cell. ";
    else
      info() << "La particule " << part_i.localId() << " de coordonnées " << m_particle_coord[part_i] << " appartient à une cellule";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_ONEPARTICLE(OneParticle, OneParticleService);
