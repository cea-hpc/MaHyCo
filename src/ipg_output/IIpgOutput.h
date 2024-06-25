// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#ifndef IIPG_OUTPUT_H
#define IIPG_OUTPUT_H


#include "arcane/ItemTypes.h"


using namespace Arcane;


/**
 * Interface du service d'écriture des sorties aux particules
 */
class IIpgOutput
{
public:
  /** Constructeur de la classe */
  IIpgOutput() {};
  /** Destructeur de la classe */
  virtual ~IIpgOutput() {};
public:
  /** Initialisation des sorties */
  virtual void initOutput() = 0;
  /** Ecriture des sorties */
  virtual void writeOutput(ParticleGroup particles) = 0;
};

#endif // IIPG_OUTPUT_H
