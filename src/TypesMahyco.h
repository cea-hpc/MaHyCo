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
     MonoAdvectionTx = 7,
     MonoAdvectionTy = 8,
     MonoAdvectionT45 = 9,
     MonoAdvectionRotation = 10,
     BiSedov = 11,
     BiTriplePoint = 12,
     BiSodCaseX = 13,
     BiSodCaseY = 14,
     BiSodCaseZ = 15,
     BiShockBubble = 16,
     BiNohTestCase = 17,
     AdvectionTx = 18,
     AdvectionTy = 19,
     AdvectionT45 = 20,
     AdvectionRotation = 21,
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
     BiSodSph = 39,
  };
  
  enum limiteur {
    minmod = 0, 
    superBee = 1,
    vanLeer = 2,
    minmodG = 3,
    superBeeG = 4,
    vanLeerG = 5,
    arithmeticG = 6,
    ultrabeeG = 7,
  }; 
#endif

