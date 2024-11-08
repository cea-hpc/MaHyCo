// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef EPPSERVICE_H
#define EPPSERVICE_H

#include "elastoMu/IElastoMu.h"
#include "elastoMu/EPP_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Représente le modèle d'élastop-plasticité
 */
class EPPService 
: public ArcaneEPPObject
{
public:
  /** Constructeur de la classe */
  EPPService(const ServiceBuildInfo & sbi)
    : ArcaneEPPObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~EPPService() {};

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
