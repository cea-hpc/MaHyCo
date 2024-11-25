// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#include "elastoMu/IElastoMu.h"
#include "SCGService.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real SCGService::getShearModulus(IMeshEnvironment* env, EnvCell ev) {
    Real gpp = options()->GPP;
    Real gpt = options()->GPT;
    Real mu0 = options()->Mu0;
    Real eta = m_density[ev] / m_density_0[ev];
    Real enerfus = options()->Eint_fus;
    Real ener_sol = m_internal_energy_0[ev];
    Real KoneOverThree = 1./3.;
    Real Coeff_press = gpp * std::max( 0., m_pressure[ev]) / (mu0*pow(eta, KoneOverThree));
    Real Coeff_temp = gpt * std::max( 0., m_temperature[ev] - 300.) / mu0;
    Real Coeff = (1 + Coeff_press + Coeff_temp);
    // adoucissement thermique en energie interne
    Real x = (m_internal_energy[ev] - ener_sol)/(enerfus - ener_sol);
    Real mu = f_ram(x) *Coeff * mu0;
    if (f_ram(x) == 0.)  info() << "mu" << mu << " x " << x << "f_ram(x)" << f_ram(x) ;
    // " temp " << m_temperature[ev] << " press " << m_pressure[ev] << " densité " <<  m_density[ev];
    // info() << " Coeff_press " <<  Coeff_press << " Coeff_temp " << Coeff_temp;
    return mu;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real SCGService::getShearModulusDerivate(IMeshEnvironment* env, EnvCell ev) {
    Real gpt = options()->GPT;
    // terme d'adoucissement thermique dans la dérivée ?
    return gpt;
}
ARCANE_REGISTER_SERVICE_SCG(SCG, SCGService);
