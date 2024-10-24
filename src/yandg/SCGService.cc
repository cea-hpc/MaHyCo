// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#include "yandg/IYandG.h"
#include "SCGService.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real SCGService::getShearModulus(IMeshEnvironment* env) { 
    Real value(0.);
    info() << "Y0" << options()->Y0();
    return value;}
Real SCGService::getElasticLimit(IMeshEnvironment* env) { 
    Real value(0.);
    return value;}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_SCG(SCG, SCGService);
