// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef JCSERVICE_H
#define JCSERVICE_H

#include "elastoMu/IElastoMu.h"
#include "elastoMu/JC_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Représente le modèle d'élastop-plasticité
 */
class JCService 
: public ArcaneJCObject
{
public:
  /** Constructeur de la classe */
  JCService(const ServiceBuildInfo & sbi)
    : ArcaneJCObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~JCService() {};

public:
   /** 
   *  Renvoie la constante Mu de l'environnement. 
   */
  virtual Real getShearModulus(IMeshEnvironment* env, EnvCell ev);
};

#endif
