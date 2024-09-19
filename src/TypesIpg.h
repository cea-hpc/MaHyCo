// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#ifndef TYPESIPG_H
#define TYPESIPG_H


struct TypesIpg
{
  enum eDrag
  {
   LinearDrag, //!< tra�n�e lin�aire en la vitesse
   QuadraticDrag  //!< tra�n�e quadratique en la vitesse
  }; 
};

#endif

