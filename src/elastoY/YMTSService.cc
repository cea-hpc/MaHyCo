// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#include "elastoY/IElastoY.h"
#include "YMTSService.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real YMTSService::getElasticLimit(IMeshEnvironment* env, EnvCell ev) {
  Real D = options()->D;
  Real T0 = options()->T0;
  Real mu0 = options()->Mu0;
  Real mu = mu0 - D / (exp( T0 /  m_temperature[ev]) -1);
  Real Sigma_a = options()->SigmaA;
  Real Epsilon_0 = options()->Epsilon0;
  Real g_0 = options()->g0;
  Real kb =  options()->kb;
  Real Pi = options()->Pi;
  Real Qi = options()->Qi;
  Real Sigma_i = options()->SigmaISurMu0;
  Real b =  options()->b;
  Real facteur = kb * m_temperature[ev] / (g_0 * pow(b,3) * mu) * std::log ( Epsilon_0 / m_plastic_deformation[ev]);
  Real hardering_reponse = pow((1.- pow(facteur , 1./Qi)), 1./Pi);
  Real yelas =  Sigma_a + hardering_reponse * Sigma_i * mu;
  return yelas;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_YMTS(YMTS, YMTSService);
