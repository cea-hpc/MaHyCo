// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef USERFILEINPUTPARTICULESSERVICE_H
#define USERFILEINPUTPARTICULESSERVICE_H

#include <filesystem>

#include "arcane/IParticleFamily.h"
#include "arcane/IMesh.h"
#include "arcane/IItemFamily.h"

#include "ipg_creation_of_particles/ICreationOfParticles.h"
#include "ipg_creation_of_particles/userfileinput/UserFileInputParticles_axl.h"


#include "ipg_creation_of_particles/utilsIpg.h" // isCoordInCell


using namespace Arcane;


/**
Interface de la creation de particules.

*/


class UserFileInputParticlesService
: public ArcaneUserFileInputParticlesObject
{
public:
  /** Constructeur de la classe */
  UserFileInputParticlesService(const ServiceBuildInfo & sbi);

  /** Destructeur de la classe */
  virtual ~UserFileInputParticlesService() {};
  
  /** 
   *  Crée les particules et leurs propriétés (poids, position, vitesse, 
   *  rayon, température) à partir des données stockées lors de initParticles.
   */
  virtual void createParticles();

  /** 
   * Etapes à faire lors de start-init : lecture du fichier utilisateur
   * et stockage des données.
   */
  virtual void initParticles();



private:
  /*
  non-active particle group, contains the particules given in the user file 
  as long as they have not yet been injected in the simulation as active particule
  */
  ParticleGroup toBeCreatedParticlesGroup;    // fixme: pointeur ??

  /*
  Add particles to the particle family
  */
  void initialize_particule_family();

  /*
  init data of the particles
  */
  void initialize_data_particule();

  /*
  time of injection of the next particle, to avoid some loops over the particles
  */
  Real t_next_part;
  
  /*
  get the time of the next injection of particle
  */
  Real get_t_next_part();

  /*
  assigne les particules à la cellule qui les contient.
  */
  void assignParticleToCell(IItemFamily* item_family, UniqueArray<ParticleEnumerator> particules_to_move, Int32UniqueArray particles_to_move_Id);

};

#endif













