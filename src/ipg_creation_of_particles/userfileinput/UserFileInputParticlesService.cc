// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "UserFileInputParticlesService.h"

using namespace Arcane;


/*---------------------------------------------------------------------------*/
/**
 * \brief Constructeur du service
 */
/*---------------------------------------------------------------------------*/
UserFileInputParticlesService::UserFileInputParticlesService(const ServiceBuildInfo & sbi) 
  : ArcaneUserFileInputParticlesObject(sbi)
{
  IItemFamily* item_family = mesh()->findItemFamily (eItemKind::IK_Particle, "AllParticles");

  // we create a new group for the particles read from the user file but not yet injected in the simulation
  toBeCreatedParticlesGroup = item_family->createGroup("toBeCreatedItem");
  // we create a new for the particles to be assigned to a cell
  toAssignCellParticlesGroup = item_family->createGroup("toAssignCellItem");
}


/*---------------------------------------------------------------------------*/
/* 
   \brief Initialisation de la création de particules dans start-init : lecture du fichier utilisateur et stockage des données d'entrée des particules

   Le fichier utilisateur est ainsi conçu :
   Chaque ligne contient une particule à créer.
   Les colonnes de données de ce fichier sont séparées par un espace et contiennent : 
   1) temps auquel la particule est injectée dans la simulation
   2) [Integer] poids initial de la particule (nombre de particules physiques représentées par la particule numérique)
   3-5) position initiale x y z
   6-8) vitesse initiale u_x u_y u_z
   9) rayon initial
   10) masse volumique initiale
   11) température initiale

   Le fichier peut contenir des lignes commentées commençant par #
   et des lignes vides.

*/
/*---------------------------------------------------------------------------*/

void UserFileInputParticlesService::initParticles()
{
  // ### we add particles to the particule family for each particle in the file
  initialize_particule_family();
  
  // ### we fill in the init data of the particules
  initialize_data_particule();

  // ### we find the time of injection of the first particle
  t_next_part = get_t_next_part();
}


/*---------------------------------------------------------------------------*/
/* 
   \brief Création des particules (dans compute-loop)

   Lorsque le temps de la simulation atteint le temps où des particules doivent
   être créées, on les change de groupe (toBeCreatedParticlesGroup -> activeParticlesGroup)
*/
/*---------------------------------------------------------------------------*/

void UserFileInputParticlesService::createParticles()
{
  // ### si au moins une particule doit être injectée à cette itération
  if ((t_next_part >= m_global_time()-m_global_deltat()) && (t_next_part < m_global_time())){
    // note: m_global_time() contient le temps t^{n+1} !!!

    IItemFamily* item_family = mesh()->findItemFamily (eItemKind::IK_Particle, "AllParticles");
    ParticleGroup activeParticlesGroup = item_family->findGroup("activeItem");

    // ### on récupère l'Id de toutes les particules à injecter
    UniqueArray<Int32> particles_to_move;
    ENUMERATE_PARTICLE (part_i, toBeCreatedParticlesGroup) {
      Real t_init = m_particle_init_time[part_i];
      if ((t_init >= m_global_time()-m_global_deltat()) && (t_init < m_global_time())) // note: m_global_time() contient t^{n+1} !
        particles_to_move.add(part_i.localId());
    }

    // ### On change de groupe les particules pour lesquelles c'est nécessaire
    if (!particles_to_move.empty()) {   // TODO: remove the if, I think it is not necessary
      info() << "We change the group of n= " << particles_to_move.size() << " particles to inject them in the simulation ";

      // on les ajoute au groupe des particules actives
      activeParticlesGroup.addItems(particles_to_move);
      // on les retire du groupe des particules à créer
      toBeCreatedParticlesGroup.removeItems(particles_to_move);
      // on les ajoute au groupe des particules pour lesquelles il faut affecter une cell
      toAssignCellParticlesGroup.addItems(particles_to_move);
    }

    // ### on assigne les nouvelles particules actives à leur cellule
    assignParticleToCell(item_family);

    // on vérifie que toutes les particules pour lesquelles il fallait affecter une cellule sont dans une cellule.
    ENUMERATE_PARTICLE (part_i, toAssignCellParticlesGroup) {
      if (!(part_i->hasCell())){
        info() << "WARNING: Particle " << part_i.localId() << " located in " << m_particle_coord[part_i] << " has no cell. We remove it from the group of active particles.";
        UniqueArray<Int32> particle_to_remove;
        particle_to_remove.add(part_i.localId());
        activeParticlesGroup.removeItems(particle_to_remove);
      }
      // else
      //   info() << "La particule " << part_i.localId() << " de coordonnées " << m_particle_coord[part_i] << " appartient à une cellule";
    }

    if (!particles_to_move.empty()) {
      // on vide le groupe des particules pour lesquelles il fallait affecter une cell
      toAssignCellParticlesGroup.removeItems(particles_to_move);
    }

    // ### we find the time of injection of the first particle
    t_next_part = get_t_next_part();
  }

}


/*---------------------------------------------------------------------------*/
/*
  Add particles to the particle family by reading them in the user-provided file
*/
/*---------------------------------------------------------------------------*/

void UserFileInputParticlesService::initialize_particule_family()
{

  // il faut ajouter la particule au groupe toBeCreatedParticlesGroup ET à la famille AllParticles
  IItemFamily* item_family = mesh()->findItemFamily (eItemKind::IK_Particle, "AllParticles");
  IParticleFamily* m_particles_family = item_family->toParticleFamily();
  
  // ### ouverture du fichier
  String  Sfilename = options()->getFilename();
  std::string filename = Sfilename.localstr();
  if (!std::filesystem::exists(filename)) {
    std::cout << "ERROR: Couldn't open the particles input data user file " << filename << std::endl;
    return;
  }
  std::ifstream user_file(filename);
  std::string one_particle;

  // ### pour chaque ligne du fichier, on crée une particule
  Integer N_particule = 0;
  while (std::getline(user_file, one_particle)) {
    if ( (!one_particle.empty()) && (one_particle[0] != '#')) {  // on vérifie que la ligne n'est pas vide ou commentée
      UniqueArray<Integer> lids({N_particule}); //local Id
      UniqueArray<Int64> uids({N_particule}); //unique Id
      info() << "Création des particules de localId " << lids.view();
      m_particles_family->addParticles(uids.view(), lids.view());
      toBeCreatedParticlesGroup.addItems(lids.view());
      N_particule++;
    }
  }
  user_file.close();

  m_particles_family->endUpdate();  // should not be forgotten

}



/*---------------------------------------------------------------------------*/
/*
  \brief Initialise les variables contenant les valeurs initiales des particules
*/
/*---------------------------------------------------------------------------*/

void UserFileInputParticlesService::initialize_data_particule()
{

  // ### On lit une deuxième fois le fichier d'input utilisateur, cette fois-ci
  //     pour en extraire toutes les valeurs.

  // ### ouverture du fichier
  String  Sfilename = options()->getFilename();
  std::string filename = Sfilename.localstr();
  if (!std::filesystem::exists(filename)) {
    std::cout << "ERROR: Couldn't open the particles input data user file " << filename << std::endl;
    return;
  }
  std::ifstream user_file(filename);
  std::string one_particle;

  // ### stockage des données initiales
  ENUMERATE_PARTICLE (part_i, toBeCreatedParticlesGroup) {

    // lecture des données d'une particule
    std::getline(user_file, one_particle);
    // on vérifie que la ligne n'est pas vide ou commentée, sinon on lit la suivante
    while ( one_particle.empty() || (one_particle[0] == '#')) {
      std::getline(user_file, one_particle);
    }

    // we want a stringstream to use the >> operator to extract the 11 data separated by a space
    std::stringstream s_one_particle(one_particle);

    // initial time, weight of the particle, positions, velocity, radius, dnesity, temperature
    std::string ti, wi, xi, yi, zi, uxi, uyi, uzi, ri, rhoi, Ti;
    s_one_particle >> ti >> wi >> xi >> yi >> zi >> uxi >> uyi >> uzi >> ri >> rhoi >> Ti;

    // stockage des valeurs initiales
    m_particle_init_time[part_i] = stod(ti);
    m_particle_weight[part_i] = stoi(wi);
    m_particle_coord[part_i] = Real3(stod(xi), stod(yi), stod(zi));
    m_particle_velocity[part_i] = Real3(stod(uxi), stod(uyi), stod(uzi));
    m_particle_radius[part_i] = stod(ri);
    m_particle_density[part_i] = stod(rhoi);
    m_particle_temperature[part_i] = stod(Ti);

    if (stod(ti) > options()->getTMaxInjection())
      info() << "WARNING : une particule ne sera pas injectée dans la simulation. Please increase option t-max-injection." ;
  }
  user_file.close();

}


/*---------------------------------------------------------------------------*/
/*
  \brief find the time of the injection of the next particle
*/
/*---------------------------------------------------------------------------*/

Real UserFileInputParticlesService::get_t_next_part()
{  

  Real tend = options()->getTMaxInjection();  // par defaut: 1e30
  
  // on parcours le groupe des particules à créer pour connaître le temps d'injection
  // de la prochaine particule.
  // pour rappel: les particules déjà créées n'appartiennent plus à ce groupe de particule
  ENUMERATE_PARTICLE (part_i, toBeCreatedParticlesGroup) {
    if (m_particle_init_time[part_i]<tend)
      tend = m_particle_init_time[part_i];    // fixme: reduction
  }

  if (tend < options()->getTMaxInjection() )
    info() << "Next particule(s) will be injected at time " << tend ;
  else
    info() << "All particles have been injected";
  
  return tend;
}



/*---------------------------------------------------------------------------*/
/*
  assigne les particules à la cellule qui les contient.
*/
/*---------------------------------------------------------------------------*/

void UserFileInputParticlesService::assignParticleToCell(IItemFamily* item_family)
{

  IParticleFamily* m_particles_family = item_family->toParticleFamily();
  ENUMERATE_PARTICLE (part_i, toAssignCellParticlesGroup) {

    // on récupère les coordonnées de la particule
    Particle particule = *part_i;
    Real3 particule_coord = m_particle_coord[part_i] ;

    // on récupère les coordonnées des noeuds de toutes le cellules
    ENUMERATE_CELL ( icell, allCells() ) {
      Cell cell = * icell;
      UniqueArray<Real3> nodes_coord;
      for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
        Real3 node_coord = m_node_coord[inode] ;
        nodes_coord.add(node_coord);
      }

      // si la particule est géométriquement dans la cellule, on associe les deux
      if (isCoordInCell(particule_coord, nodes_coord))
        m_particles_family->setParticleCell(particule, cell);
    }

  }



}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_USERFILEINPUTPARTICLES(UserFileInputParticles, UserFileInputParticlesService);
