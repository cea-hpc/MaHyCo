// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "IpgModule.h"

#include "arcane/ITimeLoopMng.h"
#include "arcane/IMesh.h"

using namespace Arcane;
using namespace Arcane::Materials;


/*---------------------------------------------------------------------------*/
/**
 * \brief Création des particules
 */
/*---------------------------------------------------------------------------*/

void IpgModule::
createParticles()
{
    info() << "*********************************";
    info() << "coucou 1";
    info() << "***********************************";
}

/*---------------------------------------------------------------------------*/
/**
 * \brief Mise à jour de la position des particule
 */
/*---------------------------------------------------------------------------*/

void IpgModule::
updateParticlePosition()
{
    info() << "*********************************";
    info() << "coucou 2";
    info() << "***********************************";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_IPG(IpgModule);
