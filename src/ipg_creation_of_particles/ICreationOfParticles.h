// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef ICREATIONOFPARTICLES_H
#define ICREATIONOFPARTICLES_H




using namespace Arcane;


/**
Interface de la creation de particules.

*/


class ICreationOfParticles
{
public:
  /** Constructeur de la classe */
  ICreationOfParticles() {};
  /** Destructeur de la classe */
  virtual ~ICreationOfParticles() {};
  
  /** 
   *  Crée les particules et leurs propriétés (poids, position, vitesse, 
   *  rayon, température). 
   */
  virtual void createParticles(IParticleFamily* m_particles_family) = 0;  // TODO: arguments nécessaires ??

};

#endif













