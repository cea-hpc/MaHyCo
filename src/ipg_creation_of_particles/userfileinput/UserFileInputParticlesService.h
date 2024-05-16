// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef USERFILEINPUTPARTICULESSERVICE_H
#define USERFILEINPUTPARTICULESSERVICE_H


#include "arcane/IParticleFamily.h"

#include "ipg_creation_of_particles/ICreationOfParticles.h"
#include "ipg_creation_of_particles/userfileinput/UserFileInputParticles_axl.h"



using namespace Arcane;


/**
Interface de la creation de particules.

*/


class UserFileInputParticlesService
: public ArcaneUserFileInputParticlesObject
{
public:
  /** Constructeur de la classe */
  UserFileInputParticlesService(const ServiceBuildInfo & sbi) 
    : ArcaneUserFileInputParticlesObject(sbi) {}
  /** Destructeur de la classe */
  virtual ~UserFileInputParticlesService() {};
  
  /** 
   *  Crée les particules et leurs propriétés (poids, position, vitesse, 
   *  rayon, température) à partir des données stockées lors de initParticles.
   */
  virtual void createParticles(IParticleFamily* m_particles_family);

  /** 
   * Etapes à faire lors de start-init : lecture du fichier utilisateur
   * et stockage des données.
   */
  virtual void initParticles();



// private:
//   m_particles_not_yet_created_family;
  
};

#endif













