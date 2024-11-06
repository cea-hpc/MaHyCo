// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#include "elastoMu/IElastoMu.h"
#include "EPPService.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real EPPService::getShearModulus(IMeshEnvironment* env, EnvCell ev) {
    Real mu =  options()->elasticCst();
    // adoucissement thermique en energie interne
    /* 
    Real enerfus = options()->Eint_fus;
    Real ener_sol = m_internal_energy_0[ev];
    Real x = (m_internal_energy[ev] - ener_sol)/(enerfus - ener_sol);
    mu *=f_ram(x); 
    */
    return mu;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_EPP(EPP, EPPService);
