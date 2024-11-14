// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

/*
Copyright 2000-2024 CEA

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include "NodeWatchingService.h"
#include "arcane/IParallelMng.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Initialisation de la surveillance temporelle du noeud                     */
/*---------------------------------------------------------------------------*/
void NodeWatchingService::init()
{
  Integer rank = subDomain()->parallelMng()->commRank();
  Real3 borneSupNoeud = options()->borneSup; 
  Real3 borneInfNoeud = options()->borneInf; 
  // Sortie Time History
  bool trouve = false;

  m_noeud_th = -1;
  ENUMERATE_NODE(inode, allNodes()){
      Node node = *inode;
      if ((borneInfNoeud.x < m_node_coord[node].x) && (m_node_coord[node].x < borneSupNoeud.x) 
      && (borneInfNoeud.y < m_node_coord[node].y) && (m_node_coord[node].y < borneSupNoeud.y)
      && (borneInfNoeud.z < m_node_coord[node].z) && (m_node_coord[node].z < borneSupNoeud.z)
      && node.isOwn() ) {
      m_noeud_th = node.localId();
      trouve = true;
      }
  }
  if (!trouve) {
    std::cout << "Pas de noeud trouvé pour le time history dans ce sous-domaine ( m_noeud_th = " << m_noeud_th() << " )" << std::endl;
  } else {
    std::cout << " Noeud " << m_noeud_th() << " trouvé dans le sous-domaine " << rank << std::endl;
  }
  std::ofstream fichier_noeud("time-history_noeud.cvs", std::ofstream::app );
  fichier_noeud << "Temps,noeud,x,y,z,vx,vy,vz"  << std::endl; 
}

/*---------------------------------------------------------------------------*/
/* Ecriture des sorties temporelles pour la surveillance du noeud            */
/*---------------------------------------------------------------------------*/
void NodeWatchingService::write()
{
 
  std::ofstream fichier_noeud("time-history_noeud.cvs", std::ofstream::app );

  Integer period = options()->periode;
  if (m_global_iteration()%period !=0) return;

  if (fichier_noeud.is_open()) {
    ENUMERATE_NODE(inode, allNodes()){
      Node node = *inode;
      if ( m_noeud_th() == node.localId()) {
          fichier_noeud << m_global_time() << ",";
          fichier_noeud << m_noeud_th() << ",";
          fichier_noeud << m_node_coord[node].x << ",";
          fichier_noeud << m_node_coord[node].y << ",";
          fichier_noeud << m_node_coord[node].z << ",";
          fichier_noeud << m_velocity[node].x << ",";
          fichier_noeud << m_velocity[node].y << ",";
          fichier_noeud << m_velocity[node].z << std::endl;
      }
    }
  }
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_NODEWATCHING(NodeWatching, NodeWatchingService);
