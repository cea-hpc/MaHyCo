﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapADIService.h"
/**
 *******************************************************************************
 * \file remapVariables()
 * \brief Calcul des variables aux mailles et aux noeuds qui ont ete projetees
 *
 * \param  varlp->Uremap2
 * \return m_fracvol_env, m_mass_fraction_env
 *         varlp->mixte, varlp->pure
 *         m_density_nplus1, m_density_env_nplus1
 *         m_internal_energy_nplus1, m_internal_energy_env_nplus1
 *         m_cell_mass, m_cell_mass_env
 *         m_node_velocity_nplus1, m_x_velocity, m_y_velocity
 *******************************************************************************
 */
void RemapADIService::remapVariables(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env) {
  
  debug() << " Entree dans remapVariables";
  mm = IMeshMaterialMng::getReference(mesh());
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ConstArrayView<IMeshEnvironment*> envs = mm->environments();
  Int32UniqueArray cells_to_add;
  Int32UniqueArray cells_to_remove;
  for (Integer index_env=0; index_env < nb_env ; index_env++) { 
    IMeshEnvironment* env = envs[index_env];
    CellGroup env_cells = env->cells();
    cells_to_add.clear();
    cells_to_remove.clear();
    Integer max_id = allCells().itemFamily()->maxLocalId();
    Int32UniqueArray cells_marker(max_id,-1);
    ENUMERATE_CELL(icell, env_cells) {
      cells_marker[icell.localId()] = 0;
    }
    ENUMERATE_CELL(icell, allCells()){
      if (((m_u_lagrange[icell][index_env] / m_euler_volume[icell]) < options()->threshold ||(m_u_lagrange[icell][nb_env + index_env] == 0. && withDualProjection))
          && cells_marker[icell.localId()] == 0) {
        // m_u_lagrange[icell][nb_env + index_env] == 0. && !withDualProjection signifie 
        // masse nulle et cas avec lagrange i.e. pas dans le cas d'une advection pure dans le vide
        cells_to_remove.add(icell.localId());
        debug() << " cell " << icell.localId() << " ( " << icell->uniqueId() << " ) " << " retirée dans l'env " << env->name();
      } else if ((m_u_lagrange[icell][index_env] / m_euler_volume[icell]) > options()->threshold 
          && cells_marker[icell.localId()] == -1) {  
        // verification que le volume normaliseé fournit une fraction de volume au-dessus du threshold             
        Real volt(0.), vol_ev_apres_normalisation(0.) ;
        for (Integer index_env_loc=0; index_env_loc < nb_env; index_env_loc++) { 
          // somme des volumes
          volt += m_u_lagrange[icell][index_env_loc] ;
        }
        vol_ev_apres_normalisation = m_u_lagrange[icell][index_env] *  m_euler_volume[icell] / volt;
        if ((vol_ev_apres_normalisation / m_euler_volume[icell]) > options()->threshold  
            && m_u_lagrange[icell][nb_env + index_env] != 0.) {
          cells_to_add.add(icell.localId());
          debug() << " cell " << icell.localId() << " ( " << icell->uniqueId() << " ) " 
            << " ajoutée dans l'env apres normalisation" << env->name() << " d'index " << index_env;
          debug() << " volume : " <<  m_u_lagrange[icell][index_env] << " fracvol " << m_u_lagrange[icell][index_env] / m_euler_volume[icell];
          debug() << " volume-apres_nomalisation " << vol_ev_apres_normalisation <<  " fracvol " << vol_ev_apres_normalisation / m_euler_volume[icell];
          debug() << " masse projetée " << m_u_lagrange[icell][nb_env + index_env];
        }
      }
    }
    
    if (!cells_to_add.empty()) {
      debug() << "ADD_CELLS to env " << env->name() << " n=" << cells_to_add.size();
      env_cells.addItems(cells_to_add);
    }
    if (!cells_to_remove.empty()){
      debug() << "REMOVE_CELLS to env " << env->name() << " n=" << cells_to_remove.size();
      env_cells.removeItems(cells_to_remove);
    }
    
  }
  mm->forceRecompute();
  UniqueArray<Real> vol_nplus1(nb_env);
  UniqueArray<Real> density_env_nplus1(nb_env);
  UniqueArray<Real> internal_energy_env_nplus1(nb_env);
  m_cell_mass.fill(0.0);
  Integer index_env;
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = * icell;   
    Real vol = m_euler_volume[cell];  // volume euler   
    m_cell_volume[cell] = m_euler_volume[cell]; // retour à la grille euler
    Real volt = 0.;
    Real masset = 0.;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    
    // info() << " cell " << cell.localId() << " calcul des masses et volumes totales " 
    //        << m_u_lagrange[cell][0] << " et " << m_u_lagrange[cell][1];
    for (Integer index_env=0; index_env < nb_env; index_env++) { 
      vol_nplus1[index_env] = m_u_lagrange[cell][index_env];
      // somme des volumes
      volt += vol_nplus1[index_env];
      // somme des masses
      masset += m_u_lagrange[cell][nb_env + index_env];
    }
    
    double volt_normalise = 0.;   
    Real unsurvolt = 1./ volt;
    // normalisation des volumes + somme 
    for (Integer index_env=0; index_env < nb_env ; index_env++) { 
      // vol_nplus1[index_env] *= vol / volt;
      vol_nplus1[index_env] *= vol * unsurvolt;
      volt_normalise += vol_nplus1[index_env];
    }
    double somme_frac = 0.;
    Real unsurvol = 1. / vol;
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;     
      index_env = ev.environmentId();        
    //      m_fracvol[ev] = vol_nplus1[index_env] / volt_normalise;
       // m_fracvol[ev] = vol_nplus1[index_env] / vol;
       m_fracvol[ev] = vol_nplus1[index_env] * unsurvol;
      if (m_fracvol[ev] < options()->threshold)
        m_fracvol[ev] = 0.;
      somme_frac += m_fracvol[ev];
    }
    // apres normamisation
    Integer matcell(0);
    Integer imatpure(-1);  
    Real unsursomme_frac(0.);
    if (somme_frac !=0.)  unsursomme_frac = 1. / somme_frac;
    index_env = 0;  
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      index_env = ev.environmentId();  
      // m_fracvol[ev] /= somme_frac; 
      m_fracvol[ev] *= unsursomme_frac;
      if (m_fracvol[ev] > 0.) {
        matcell++;
        imatpure = index_env;
      }
    }
    if (matcell > 1) {
      m_est_mixte[cell] = 1;
      m_est_pure[cell] = -1;
    } else {
      m_est_mixte[cell] = 0;
      m_est_pure[cell] = imatpure;
    }
    // on ne recalcule par les mailles à masses nulles < threshold au carrée - cas advection
    // on enleve les petits fractions de volume aussi sur la fraction
    // massique et on normalise
    Real fmasset(0.);
    if (masset > options()->threshold*options()->threshold) {
      Real unsurmasset = 1./  masset;
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;  
        index_env = ev.environmentId();  
        // m_mass_fraction[ev] = m_u_lagrange[cell][nb_env + index_env] / masset;
        m_mass_fraction[ev] = m_u_lagrange[cell][nb_env + index_env] * unsurmasset;
        if (m_fracvol[ev] < options()->threshold) {
          m_mass_fraction[ev] = 0.;
        }
        fmasset += m_mass_fraction[ev];
      }
      if (fmasset != 0.) {
        Real unsurfmasset = 1. / fmasset;
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;  
          // m_mass_fraction[ev] /= fmasset;
          m_mass_fraction[ev] *= unsurfmasset;
        }
      }
    }
    Real density_nplus1 = 0.;   
    
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell; 
      index_env = ev.environmentId();  
      density_env_nplus1[index_env] = 0.;
      if (m_fracvol[ev] > options()->threshold && vol_nplus1[index_env] != 0.)
          density_env_nplus1[index_env] = m_u_lagrange[cell][nb_env + index_env] 
            / vol_nplus1[index_env];
      density_nplus1 += m_fracvol[ev] * density_env_nplus1[index_env];
      if (cell.localId() == -1) pinfo() << cell.localId() << " calcul density moy env " << index_env << " et density " << density_env_nplus1[index_env] << " fraction " << m_fracvol[ev] << " m_u_lagrange[cell][nb_env + index_env] " <<  m_u_lagrange[cell][nb_env + index_env] << " et vol " << vol_nplus1[index_env];
      // 1/density_nplus1 += m_mass_fraction_env(cCells)[imat] / density_env_nplus1[imat];  
    }
    Real energie_nplus1 = 0.;
    Real pseudo_nplus1 = 0.;
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell; 
      index_env = ev.environmentId();  
      internal_energy_env_nplus1[index_env] = 0.;
      
      if (m_fracvol[ev] > options()->threshold && m_u_lagrange[cell][nb_env + index_env] != 0.) {
        internal_energy_env_nplus1[index_env] =
          m_u_lagrange[cell][2 * nb_env + index_env] / m_u_lagrange[cell][nb_env + index_env];
         /* indice de phase 
          m_frac_phase1[ev] = m_u_lagrange[cell][3 * nb_env + index_env] / m_u_lagrange[cell][nb_env + index_env];
          m_frac_phase2[ev] = m_u_lagrange[cell][4 * nb_env + index_env] / m_u_lagrange[cell][nb_env + index_env];
          m_frac_phase3[ev] = m_u_lagrange[cell][5 * nb_env + index_env] / m_u_lagrange[cell][nb_env + index_env]; */
      }
    }
    // mise à jour des valeurs moyennes aux allCells
    // densite
    m_density[cell] = density_nplus1;
    if (cell.localId() == -1) pinfo() << cell.localId() << " apres proj density " << m_density[cell];
    // pinfo() << cell.localId() << " apres proj " << m_u_lagrange[cell];
    // recalcul de la masse
    m_cell_mass[cell] = m_euler_volume[cell] * density_nplus1;
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell; 
      index_env = ev.environmentId();  
      m_cell_mass[ev] = m_mass_fraction[ev] * m_cell_mass[cell];
      // recuperation de la pseudo projetee
      // m_pseudo_viscosity[ev] = m_u_lagrange[cell][3 * nb_env + 4] / vol;
      m_pseudo_viscosity[ev] = m_u_lagrange[cell][7 * nb_env ] * unsurvol;
      // recuperation de la densite
      m_density[ev] = density_env_nplus1[index_env];
      // recuperation de l'energie
      m_internal_energy[ev] = internal_energy_env_nplus1[index_env];
      // conservation energie totale
      // delta_ec : energie specifique
      // m_internal_energy_env[ev] += delta_ec;
      // energie interne totale
      energie_nplus1 += m_mass_fraction[ev] * m_internal_energy[ev];
      pseudo_nplus1 += m_fracvol[ev] * m_pseudo_viscosity[ev];
      if (cell.localId() == -1) pinfo() << cell.localId() << " apres proj env " << index_env << " et density " << m_density[ev];
    }
    // energie interne
    m_internal_energy[cell] = energie_nplus1;
    // pseudoviscosité
    m_pseudo_viscosity[cell] = pseudo_nplus1;
    
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;
      index_env = ev.environmentId();  

      if (m_density[ev] < 0. || m_internal_energy[ev] < 0.) {
        
        pinfo() << " cell " << cell.localId() << " --energy ou masse negative pour l'environnement "
               << index_env;
        pinfo() << " energie interne env " << m_internal_energy[ev] 
               << " cell " << m_internal_energy[cell] ;
        pinfo() << " densite env " << m_density[ev] 
               << " cell " << m_density[cell] ;
        pinfo() << " fraction vol env " << m_fracvol[ev];
        pinfo() << " fraction massique env " <<  m_mass_fraction[ev];
        m_internal_energy[ev] = 0.;
        m_density[ev] = 0.;
        m_fracvol[ev] = 0.;
        m_mass_fraction[ev]=0.;
      }
      if (m_density[ev] != m_density[ev] || m_internal_energy[ev] != m_internal_energy[ev]) {
        pinfo() << " cell NAN " << cell.localId() << " --energy ou masse NAN "
               << index_env;
        pinfo() << " energie interne env " << m_internal_energy[ev] 
               << " cell " << m_internal_energy[cell] ;
        pinfo() << " densite env " << m_density[ev] 
               << " cell " << m_density[cell] ;
        pinfo() << " fraction vol env " << m_fracvol[ev];
        pinfo() << " fraction massique env " <<  m_mass_fraction[ev];
        exit(1);
      }
    } 
  }   
  if (withDualProjection) {
  // variables aux noeuds
    ENUMERATE_NODE(inode, allNodes()){
      Real mass_nodale_proj = m_u_dual_lagrange[inode][3];
      if (mass_nodale_proj != 0.) {
        m_velocity[inode].x = m_u_dual_lagrange[inode][0] / mass_nodale_proj;
        m_velocity[inode].y = m_u_dual_lagrange[inode][1] / mass_nodale_proj;
        m_velocity[inode].z = m_u_dual_lagrange[inode][2] / mass_nodale_proj;
      } 
    }
  }
  // recalcule de la masse mass_nodale et 
  m_node_mass.fill(0.);
  Real one_over_nbnode = dimension == 2 ? .25  : .125 ;
  ENUMERATE_CELL(icell, allCells()){    
    Cell cell = * icell;
    Real contrib_node_mass = one_over_nbnode * m_cell_mass[cell];
    for( NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode){
      m_node_mass[inode] += contrib_node_mass; 
    }
    // des volumes matieres
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;        
        Cell cell = ev.globalCell();
        m_cell_volume[ev] = m_fracvol[ev] * m_cell_volume[cell];
      }

    }
  } 
  if (isEuler()) {
      // info() << " recopie des valeurs initiales pour les sorties ";
      m_node_coord.copy(m_node_coord_0);
  }
  m_node_mass.synchronize();
  // conservation energie totale lors du remap
  if (hasConservationEnergieTotale()) {
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;
      Real delta_ec(0.);
      Real ec_proj(0.);
      Real ec_reconst(0.);
      AllEnvCell all_env_cell = all_env_cell_converter[cell];
      ENUMERATE_NODE(inode, cell->nodes()) {
        if (m_u_dual_lagrange[inode][3] != 0.) {
          ec_proj = m_u_dual_lagrange[inode][4] / m_u_dual_lagrange[inode][3];
          ec_reconst = 0.5 * m_velocity[inode].squareNormL2();
          delta_ec += 0.25 * ( ec_proj - ec_reconst);
        }
      }
      delta_ec = std::max(0., delta_ec);
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell; 
        m_internal_energy[ev] += m_mass_fraction[ev] * delta_ec;
      }
      if (all_env_cell.nbEnvironment() !=1) {
        // depot sur la valeur moyenne
        m_internal_energy[cell] += delta_ec;
      }
    }
  }
}
/**
 *******************************************************************************
 * \file computeVariablesGlobalesT()
 * \brief Calcul de l'energie totale et la masse initiale du systeme apres projection
 *
 * \param  m_cell_velocity_nplus, m_density_nplus, m_euler_volume
 * \return m_total_energy_T, m_global_masse_T,
 *         m_global_total_energy_T, m_total_masse_T
 *
 *******************************************************************************
 */
// void Vnr::computeVariablesGlobalesT() noexcept {
//   m_global_total_energy_T = 0.;
//   int nbmat = options->nbmat;
//   Kokkos::parallel_for(
//       "remapVariables", nbCells, KOKKOS_LAMBDA(const int& cCells) {
// 	const Id cId(cCells);
// 	const auto nodesOfCellC(mesh->getNodesOfCell(cId));
// 	const size_t nbNodesOfCellC(nodesOfCellC.size());
// 	double ec_reconst(0.);
// 	for (size_t pNodesOfCellC = 0; pNodesOfCellC < nbNodesOfCellC;
// 	     pNodesOfCellC++) {
// 	  const Id pId(nodesOfCellC[pNodesOfCellC]);
// 	  const size_t pNodes(pId);
// 	  ec_reconst += 0.25 * 0.5 *
// 	    (m_node_velocity_nplus1(pNodes)[0] * m_node_velocity_nplus1(pNodes)[0] +
// 	     m_node_velocity_nplus1(pNodes)[1] * m_node_velocity_nplus1(pNodes)[1]);
// 	}
//         m_total_energy_T(cCells) = m_density_nplus1(cCells) * m_euler_volume(cCells) *
// 	  (m_internal_energy_nplus1(cCells) + ec_reconst);
//         m_total_masse_T(cCells) = 0.;
//         for (int imat = 0; imat < nbmat; imat++)
//           m_total_masse_T(cCells) += m_density_env_nplus1(cCells)[imat] *
//                                      m_euler_volume(cCells) *
//                                      m_fracvol_env(cCells)[imat];
//         // m_mass_fraction_env(cCells)[imat] * (density_nplus1 * vol) ; //
//         // m_density_env_nplus1(cCells)[imat] * vol_nplus1[imat];
//       });
//   double reductionE(0.), reductionM(0.);
//   {
//     Kokkos::Sum<double> reducerE(reductionE);
//     Kokkos::parallel_reduce("reductionE", nbCells,
//                             KOKKOS_LAMBDA(const int& cCells, double& x) {
//                               reducerE.join(x, m_total_energy_T(cCells));
//                             },
//                             reducerE);
//     Kokkos::Sum<double> reducerM(reductionM);
//     Kokkos::parallel_reduce("reductionM", nbCells,
//                             KOKKOS_LAMBDA(const int& cCells, double& x) {
//                               reducerM.join(x, m_total_masse_T(cCells));
//                             },
//                             reducerM);
//   }
//   m_global_total_energy_T = reductionE;
//   m_global_total_masse_T = reductionM;
// }
