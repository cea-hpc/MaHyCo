// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#include "YJCService.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real YJCService::getElasticLimit(IMeshEnvironment* env, EnvCell ev) {
    Real a = options()->A;
    Real b = options()->B;
    Real c = options()->C;
    Real n = options()->n;
    Real m = options()->m;
    Real tm = options()->Tm;
    Real eps = options()->Epsilon_init;
    Real Teff = std::min(std::max(0., (m_temperature[ev] - 300.)/(tm - 300.)), 1.);

    Real yelas = a + b * pow(m_plastic_deformation[ev], n)
        * (1 + c * std::log( std::max(m_plastic_deformation_velocity[ev],eps)/eps))
        * (1-Teff);
    return yelas;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_YJC(YJC, YJCService);
