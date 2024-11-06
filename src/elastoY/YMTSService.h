// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#ifndef YMTSSERVICE_H
#define YMTSSERVICE_H

#include "elastoY/IElastoY.h"
#include "elastoY/YMTS_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Représente le modèle d'élastop-plasticité
 */
class YMTSService 
: public ArcaneYMTSObject
{
public:
  /** Constructeur de la classe */
  YMTSService(const ServiceBuildInfo & sbi)
    : ArcaneYMTSObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~YMTSService() {};

public:
  /** 
   *  Renvoie la La limite de l'environnement. 
   */
  virtual Real getElasticLimit(IMeshEnvironment* env, EnvCell ev);
};

#endif
