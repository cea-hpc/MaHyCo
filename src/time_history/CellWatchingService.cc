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

#include "CellWatchingService.h"
#include "arcane/IParallelMng.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/* Initialisation de la surveillance temporelle du noeud                     */
/*---------------------------------------------------------------------------*/
void CellWatchingService::init()
{
  mm = IMeshMaterialMng::getReference ( mesh() );
  Integer rank = subDomain()->parallelMng()->commRank();

  Real3 borneSup = options()->borneSup; 
  Real3 borneInf = options()->borneInf; 
  // Sortie Time History
  m_maille_th = -1;
  bool trouve = false;
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    if ((borneInf.x < m_cell_coord[cell].x) && (m_cell_coord[cell].x < borneSup.x) 
      && (borneInf.y < m_cell_coord[cell].y) && (m_cell_coord[cell].y < borneSup.y)
      && (borneInf.z < m_cell_coord[cell].z) && (m_cell_coord[cell].z < borneSup.z)
      && cell.isOwn() ) {
      m_maille_th = cell.localId();
      trouve = true;
    }
  }
  if (!trouve) {
      std::cout << "Pas de maille trouvé pour le time history dans ce sous-domaine ( m_maille_th = " << m_maille_th() << " )" << std::endl;
  } else {
      std::cout << " Maille " << m_maille_th() << " trouvé dans le sous-domaine " << rank << std::endl;
  }   
  
}

/*---------------------------------------------------------------------------*/
/* Ecriture des sorties temporelles pour la surveillance du noeud            */
/*---------------------------------------------------------------------------*/
void CellWatchingService::write()
{
  std::ofstream fichier("time-history.csv", std::ofstream::app );
  Integer period = options()->periode;
  if (m_global_iteration()%period !=0) return;

  if (fichier.is_open()) {
    fichier << " Temps, densité,energie_interne,pression,temperature,vitesse_son,fraction_ph1,fraction_ph2,fraction_ph3,norm_s,sxx,syy,sxy,syz,sxz"  << std::endl;
    ENUMERATE_ENV(ienv,mm){
      IMeshEnvironment* env = *ienv;
      ENUMERATE_ENVCELL(ienvcell,env) {
        EnvCell ev = *ienvcell;   
        Cell cell = ev.globalCell();
        if ( m_maille_th() == cell.localId()) {
          fichier << m_global_time() << ",";
          fichier << m_maille_th() << ",";
          // fichier << m_cell_coord[cell] << " ";
          fichier << m_density[ev] << ",";
          fichier << m_internal_energy[ev] << ",";
          fichier << m_pressure[ev] << ",";
          fichier << m_temperature[ev] << ",";
          fichier << m_sound_speed[ev] << ",";
          fichier << m_frac_phase1[ev] << ",";
          fichier << m_frac_phase2[ev] << ",";
          fichier << m_frac_phase3[ev] << ",";
          fichier << m_strain_norm[ev] << ",";
          fichier << m_strain_tensor_xx[ev] << ",";
          fichier << m_strain_tensor_yy[ev] << ",";
          fichier << m_strain_tensor_xy[ev] << ",";
          fichier << m_strain_tensor_yz[ev] << ",";
          fichier << m_strain_tensor_xz[ev] << std::endl;
        }
      }
    }   
  }
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_CELLWATCHING(CellWatching, CellWatchingService);
