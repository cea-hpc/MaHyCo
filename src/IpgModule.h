// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef IPGMODULE_H
#define IPGMODULE_H

#include "arcane/IParticleFamily.h"
#include "ipg_output/IIpgOutput.h"
#include "ipg_creation_of_particles/ICreationOfParticles.h"

#include "Ipg_axl.h"

using namespace Arcane;

/**
 * Représente un module pour le traitement des particules
 * - création de particules avec une vitesse donnée
 * - calcul de la trajectoire des particules
 * 
 */


class IpgModule
: public ArcaneIpgObject
{
 public:
  /** Constructeur de la classe */
  IpgModule(const ModuleBuildInfo& mbi);
  
  /** Destructeur de la classe */
  ~IpgModule() {}

  /** Création des sorties pour les particules */
  virtual void initParticleOutput();
  
  /** Création des particules */
  virtual void injectParticles();

  /** Initialisations éventuellement nécessaires dans start-init pour créer des particules */
  virtual void initInjectParticles();

  /** Mise à jour de la position des particules */
  virtual void updateParticlePosition();

  /** Ecriture des sorties des particules */
  virtual void writeParticleOutput();

 private:
  // protected:
  IParticleFamily* m_particles_family;
};

#endif
