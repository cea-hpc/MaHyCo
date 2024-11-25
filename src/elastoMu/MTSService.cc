// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#include "elastoMu/IElastoMu.h"
#include "MTSService.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real MTSService::getShearModulus(IMeshEnvironment* env, EnvCell ev) {
    Real D = options()->D;
    Real T0 = options()->T0;
    Real mu0 = options()->Mu0;
    Real mu = mu0 - D / (exp( T0 /  m_temperature[ev]) -1);
    // adoucissement thermique en energie interne
    Real enerfus = options()->Eint_fus;
    Real ener_sol = m_internal_energy_0[ev];
    Real x = (m_internal_energy[ev] - ener_sol)/(enerfus - ener_sol);
    mu *=f_ram(x);
    // info() << "mu" << mu << " temp " << m_temperature[ev] << " press " << m_pressure[ev] << " densité " <<  m_density[ev];
    // info() << " Coeff_press " <<  Coeff_press << " Coeff_temp " << Coeff_temp;
    return mu;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real MTSService::getShearModulusDerivate(IMeshEnvironment* env, EnvCell ev) {
    Real D = options()->D;
    Real T0 = options()->T0;
    Real muprime = -D*T0*(exp(T0 /  m_temperature[ev])) /  pow(m_temperature[ev],2.) / pow((exp(T0 /  m_temperature[ev]) -1),2.);
    return muprime;
}
ARCANE_REGISTER_SERVICE_MTS(MTS, MTSService);
