// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#ifndef TYPESIPG_H
#define TYPESIPG_H


struct TypesIpg
{
  enum eDrag
  {
   LinearDrag, //!< traînée linéaire en la vitesse
   QuadraticDrag  //!< traînée quadratique en la vitesse
  }; 
};

#endif

