// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef ONEPARTICLESERVICE_H
#define ONEPARTICLESERVICE_H


#include "arcane/IParticleFamily.h"
#include "arcane/utils/Real3.h"

#include "ipg_creation_of_particles/ICreationOfParticles.h"
#include "ipg_creation_of_particles/oneparticle/OneParticle_axl.h"



using namespace Arcane;


/**
Interface de la creation de particules.

*/


class OneParticleService
: public ArcaneOneParticleObject
{
public:
  /** Constructeur de la classe */
  OneParticleService(const ServiceBuildInfo & sbi) 
    : ArcaneOneParticleObject(sbi) {}
  /** Destructeur de la classe */
  virtual ~OneParticleService() {};
  
  /** 
   *  Crée les particules et leurs propriétés (poids, position, vitesse, 
   *  rayon, température). 
   */
  virtual void createParticles(IParticleFamily* m_particles_family);

  /** 
   * Rien à faire pour ce service.
   */
  void initParticles(){}
};

#endif













