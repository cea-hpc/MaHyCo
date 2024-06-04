// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef ISPRAYTYPE_H
#define ISPRAYTYPE_H




using namespace Arcane;


/**
Interface de la creation de particules.

*/


class ISprayType
{
public:
  /** Constructeur de la classe */
  ISprayType() {};
  /** Destructeur de la classe */
  virtual ~ISprayType() {};

  /** 
      initialisation du solver pour les particules
   */
  virtual void initSolverParticles() = 0;

  /** 
      correction de la vitesse du fluide (prise en compte de la traînée due aux particules).
   */
  virtual void correctFluidVelocity() = 0;

  /** 
      mise à jour de la vitesse des particules
   */
  virtual void updateParticleVelocity() = 0;

  /** 
      mise à jour de la position des particules
   */
  virtual void updateParticlePosition() = 0;
};

#endif













