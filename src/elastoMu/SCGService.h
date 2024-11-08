// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef SCGSERVICE_H
#define SCGSERVICE_H

#include "elastoMu/IElastoMu.h"
#include "elastoMu/SCG_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Représente le modèle d'élastop-plasticité
 */
class SCGService 
: public ArcaneSCGObject
{
public:
  /** Constructeur de la classe */
  SCGService(const ServiceBuildInfo & sbi)
    : ArcaneSCGObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~SCGService() {};

public:
   /** 
   *  Renvoie le module de cisaillement (Mu) de l'environnement. 
   */
  virtual Real getShearModulus(IMeshEnvironment* env, EnvCell ev);
   /** 
   *  Renvoie la dérivé (par rapport à la température) du module de cisaillement (Mu) de l'environnement. 
   */
  virtual Real getShearModulusDerivate(IMeshEnvironment* env, EnvCell ev);
};

#endif
