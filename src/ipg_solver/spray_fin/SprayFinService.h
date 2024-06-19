// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef SPRAYFINSERVICE_H
#define SPRAYFINSERVICE_H


#include "arcane/IParticleFamily.h"
#include "arcane/IMesh.h"
#include "arcane/IItemFamily.h"
#include "arcane/utils/Real3.h"

#include "ipg_solver/ISprayType.h"
#include "ipg_solver/spray_fin/SprayFin_axl.h"

#include "mathUtils.h"  // Pi
#include "ipg_creation_of_particles/utilsIpg.h" // isCoordInCell


//#include <arcane/core/ITimeHistoryMng.h>
//#include <arcane/core/GlobalTimeHistoryAdder.h>


using namespace Arcane;


/**
Service de résolutoin des équations des spray fins.

*/


class SprayFinService
: public ArcaneSprayFinObject
{
public:

  /** Constructeur de la classe */
  SprayFinService(const ServiceBuildInfo & sbi)
    : ArcaneSprayFinObject(sbi) {}

  /** Destructeur de la classe */
  virtual ~SprayFinService() {};
    
  /** 
      initialisation du solver pour les particules
   */
  virtual void initSolverParticles();

  /** 
      correction de la vitesse du fluide (prise en compte de la traînée due aux particules).
   */
  virtual void correctFluidVelocity();

  /** 
      mise à jour de la vitesse des particules
   */
  virtual void updateParticleVelocity();

  /** 
      mise à jour de la position des particules
   */
  virtual void updateParticlePosition();

 
private:
  IParticleFamily* m_particles_family;
  ParticleGroup activeParticlesGroup;
  /** 
      calcul du coefficient de traînée
   */
  virtual Real computeCd();

};



#endif













