// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef NEWHYPOMODELSERVICE_H
#define NEWHYPOMODELSERVICE_H

#include "elasto/IElasto.h"
#include "yandg/IYandG.h"
#include "elasto/NewHypoModel_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Représente le modèle d'élastop-plasticité
 */
class NewHypoModelService 
: public ArcaneNewHypoModelObject
{
public:
  /** Constructeur de la classe */
  NewHypoModelService(const ServiceBuildInfo & sbi)
    : ArcaneNewHypoModelObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~NewHypoModelService() {};

public:
  /** 
   *  Initialise de l'elasto-plasticité au groupe de mailles passé en argument
   */
  virtual void initElasto(IMeshEnvironment* env);
  /** 
   *  Calcul des gradients de vitesse
   */
  virtual void ComputeVelocityGradient(Real delta_t);
    /** 
   *  Calcul du tenseur de déformation et de rotation
   */
  virtual void ComputeDeformationAndRotation();
  /** 
   *  Calcul du deviateur elasto-plastique au groupe de mailles passé en argument
   */
  virtual void ComputeElasticity(IMeshEnvironment* env, Real delta_t, Integer dim);
  /** 
   *  Calcul de la plasticité au groupe de mailles passé en argument
   */
  virtual void ComputePlasticity(IMeshEnvironment* env, Real delta_t, Integer dim);
  /** 
   *  Calcul du travail elasto-plastique
   */
  virtual void ComputeElastoEnergie(IMeshEnvironment* env, Real delta_t);
  
};

#endif
