#ifndef TYPESMAHYCO_H
#define TYPESMAHYCO_H

#include <arcane/ItemGroup.h>
#include "eos/IEquationOfState.h"

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

#endif

