// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#include "yandg/IYandG.h"
#include "EPPService.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real EPPService::getShearModulus(IMeshEnvironment* env, EnvCell ev) {
    // info() << options()->elasticCst();
    return options()->elasticCst();}
Real EPPService::getElasticLimit(IMeshEnvironment* env, EnvCell ev) {
    // info() << options()->limitElastic();
    return options()->limitElastic();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_EPP(EPP, EPPService);
