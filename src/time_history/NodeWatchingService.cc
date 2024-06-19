// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "NodeWatchingService.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/* Initialisation de la surveillance temporelle du noeud                     */
/*---------------------------------------------------------------------------*/
void NodeWatchingService::init()
{
   
}

/*---------------------------------------------------------------------------*/
/* Ecriture des sorties temporelles pour la surveillance du noeud            */
/*---------------------------------------------------------------------------*/
void NodeWatchingService::write()
{

}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_NODEWATCHING(NodeWatching, NodeWatchingService);
