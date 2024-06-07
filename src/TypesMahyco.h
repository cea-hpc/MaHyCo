#ifndef TYPESMAHYCO_H
#define TYPESMAHYCO_H

#include <arcane/ItemGroup.h>
#include "eos/IEquationOfState.h"
#include "elasto/IElasto.h"
#include "casTest/IInitialisations.h"
#include "Remap/IRemap.h"

#include "mathUtils.h"   //  Pi

struct TypesMahyco
{
  enum eBoundaryCondition
  {
    VelocityX, //!< Vitesse X fixée
    VelocityY, //!< Vitesse Y fixée
    VelocityZ, //!< Vitesse Z fixée
    OnFileVelocity, //!< Vitesse analytique
    Density, //!< Densité fixée
    Energy, //!< Energy fixée
    SuperGaussianEnergy, 
    Pressure, //!< Pression fixée
    GeometricPressure, //!< Pression fixée
    LinearPressure, //!< Pression linaire suivant X,Y,Z et t
    OnFilePressure, //!< Pression dans un fichier en fonction du temps
    SuperGaussianPressure, //!< Pression en supergaussienne suivant X,Y,Z
    ContactHerzPressure, //!< Pression contact de Herz (loi pression sqrt(1-r2))
    Unknown //!< Type inconnu
  }; 
  enum eEnergyDepot
  {
    DepotConstant, //!< Energy fixée constante
    DepotLineaire, //!< Energy linaire suivant X,Y,Z et t
    DepotSuperGaussian, //!< Energy en supergaussienne suivant X,Y,Z
    Inconnu //!< Type inconnu
  };
};

enum Test {
    // cas test
     OneCell = 0,
     Sedov = 1,
     MonoTriplePoint = 2,
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
     InstaRTX = 40,
     InstaRTY = 41,
     InstaRTZ = 42,
     BiInstaRTX = 43,
     BiInstaRTY = 44,
     BiInstaRTZ = 45,
     Vortex = 46,
     ChocBulle = 47,
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

