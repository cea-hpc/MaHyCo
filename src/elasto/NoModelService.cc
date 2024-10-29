// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#include "NoModelService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/* Mise à zero de tous les tenseurs ou vecteurs liés à l'elasto-plasticité   */
/*---------------------------------------------------------------------------*/

void NoModelService::initElasto(IMeshEnvironment* env)
{
     m_velocity_gradient.fill(Real3x3::zero());
     m_spin_rate.fill(Real3::zero());
     m_deformation_rate.fill(Real3x3::zero());
     m_strain_tensor.fill(Real3x3::zero());
     m_strain_tensor_n.fill(Real3x3::zero());
    
}

/*---------------------------------------------------------------------------*/
/* Calcul des gradients de vitesses à la cell */
/*---------------------------------------------------------------------------*/

void NoModelService::ComputeVelocityGradient(Real delta_t)
{

}
/*---------------------------------------------------------------------------*/
/* Calcul du taux de déformations et du vecteur de rotation                  */
/*---------------------------------------------------------------------------*/

void NoModelService::ComputeDeformationAndRotation()
{

}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NoModelService::ComputeElasticity(IMeshEnvironment* env, Real delta_t, Integer dim)
{
 
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NoModelService::ComputePlasticity(IMeshEnvironment* env, Real delta_t, Integer dim)
{
 
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NoModelService::ComputeElastoEnergie(IMeshEnvironment* env, Real delta_t)
{

}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_NOMODEL(NoModel, NoModelService);
