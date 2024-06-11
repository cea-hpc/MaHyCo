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
 * \brief Mise à jour des cellules auxquelles appartiennent les particules
 */
/*---------------------------------------------------------------------------*/

void IpgModule::
updateParticleCell()
{
    ENUMERATE_PARTICLE (ipart, activeParticlesGroup) {

    // on récupère les coordonnées de la particule et de sa cellule
    Particle particle = *ipart;
    Real3 particule_coord = m_particle_coord[ipart] ;
    
    Cell cell = particle.cell();
    UniqueArray<Real3> nodes_coord;
    for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
      Real3 node_coord = m_node_coord[inode] ;
      nodes_coord.add(node_coord);
    }

    // on vérifie si la particule est géométriquement dans la cellule qui lui est associée
    // si ce n'est pas le cas, on parcourt toutes les cellules pour lui associer la bonne
    if (!(isCoordInCell(particule_coord, nodes_coord))){

      // on doit également vérifier que la cellule n'es pas sortie du domaine
      bool particle_is_out = 1;
      
      ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        UniqueArray<Real3> nodes_coord_new;
        for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
          Real3 node_coord_new = m_node_coord[inode] ;
          nodes_coord_new.add(node_coord_new);
        }

        // si la particule est géométriquement dans la cellule, on associe les deux
        // et on précise que la particule n'est pas en dehors du domaine
        if (isCoordInCell(particule_coord, nodes_coord_new)){
          m_particles_family->setParticleCell(particle, cell);
          particle_is_out=0;
          info() << "Particle " << ipart.localId() << " with position " << m_particle_coord[ipart] << " has a new cell with nodes coordinates " << nodes_coord_new;
        }
      }

      // si la particules est sortie, on la retire du groupe des particules actives
      // mais il n'y a pas de removeParticles dans arcane pour supprimer la particule de m_particle_family
      UniqueArray<Int32> particle_to_remove;
      if (particle_is_out){
        info() << "Particle " << ipart.localId() << " with position " << m_particle_coord[ipart] << " is out of the domain and has been removed.";
        particle_to_remove.add(ipart.localId());
        activeParticlesGroup.removeItems(particle_to_remove);
      }

    }
  }

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
