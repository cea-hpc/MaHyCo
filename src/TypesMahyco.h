#ifndef TYPESMAHYCO_H
#define TYPESMAHYCO_H

#include <arcane/ItemGroup.h>
#include "eos/IEquationOfState.h"
#include "casTest/IInitialisations.h"
#include "Remap/IRemap.h"

struct TypesMahyco
{
  enum eBoundaryCondition
  {
    VelocityX, //!< Vitesse X fixée
    VelocityY, //!< Vitesse Y fixée
    VelocityZ, //!< Vitesse Z fixée
    Unknown //!< Type inconnu
  }; 
};

const double Pi = 3.14159265359;

enum Test {
    // cas test
     UnitTestCase = 0,
     Sedov = 1,
     TriplePoint = 2,
     SodCaseX = 3,
     SodCaseY = 4,
     SodCaseZ = 5,
     NohTestCase = 6,
     AdvectionX = 7,
     AdvectionY = 8,
     AdvectionVitX = 9,
     AdvectionVitY = 10,
     BiSedov = 11,
     BiTriplePoint = 12,
     BiSodCaseX = 13,
     BiSodCaseY = 14,
     BiSodCaseZ = 15,
     BiShockBubble = 16,
     BiNohTestCase = 17,
     BiAdvectionX = 18,
     BiAdvectionY = 19,
     BiAdvectionVitX = 20,
     BiAdvectionVitY = 21,
     BiImplosion = 22,
     MonoRiderTx = 23,
     MonoRiderTy = 24,
     MonoRiderT45 = 25,
     MonoRiderRotation = 26,
     MonoRiderVortex = 27,
     MonoRiderDeformation = 28,
     MonoRiderVortexTimeReverse = 29,
     MonoRiderDeformationTimeReverse = 30,
     RiderTx = 31,
     RiderTy = 32,
     RiderT45 = 33,
     RiderRotation = 34,
     RiderVortex = 35,
     RiderDeformation = 36,
     RiderVortexTimeReverse = 37,
     RiderDeformationTimeReverse = 38,
  };
  
  enum limiteur {
    minmod, 
    superBee,
    vanLeer ,
    minmodG ,
    superBeeG,
    vanLeerG,
    arithmeticG,
    ultrabeeG,
  }; 
#endif

