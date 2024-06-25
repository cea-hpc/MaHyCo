// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#ifndef ONEPARTICLESERVICE_H
#define ONEPARTICLESERVICE_H


#include "arcane/IParticleFamily.h"
#include "arcane/IMesh.h"
#include "arcane/IItemFamily.h"
#include "arcane/utils/Real3.h"

#include "ipg_creation_of_particles/ICreationOfParticles.h"
#include "ipg_creation_of_particles/oneparticle/OneParticle_axl.h"


#include "ipg_creation_of_particles/utilsIpg.h" // isCoordInCell


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
  virtual void createParticles();

  /** 
   * Rien à faire pour ce service.
   */
  void initParticles(){}

private:

  /** 
   * affecte la particule à la cellule à laquelle elle appartient
   */
  virtual void assignParticleToCell();
  
};

#endif













