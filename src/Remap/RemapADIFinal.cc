// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapADIService.h"
#include "accenv/AcceleratorUtils.h"

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
  
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans remapVariables";
  mm = IMeshMaterialMng::getReference(mesh());
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  PROF_ACC_BEGIN("cellStatus");

  bool to_add_rm_cells = false;
#if 0
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
      if ((m_u_lagrange[icell][index_env] / m_euler_volume[icell]) < options()->threshold 
          && cells_marker[icell.localId()] == 0) {
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
            << " ajoutée dans l'env apres normalisation" << env->name();
          debug() << " volume : " <<  m_u_lagrange[icell][index_env] << " fracvol " << m_u_lagrange[icell][index_env] / m_euler_volume[icell];
          debug() << " volume-apres_nomalisation " << vol_ev_apres_normalisation <<  " fracvol " << vol_ev_apres_normalisation / m_euler_volume[icell];
          debug() << " masse projetée " << m_u_lagrange[icell][nb_env + index_env];
        }
      }
    }
    
    if (!cells_to_add.empty()) {
      pinfo() << "ADD_CELLS to env " << env->name() << " n=" << cells_to_add.size() 
        << " ITERATION " << globalIteration();
      env_cells.addItems(cells_to_add);
      to_add_rm_cells = true;
    }
    if (!cells_to_remove.empty()){
      pinfo() << "REMOVE_CELLS to env " << env->name() << " n=" << cells_to_remove.size()
        << " ITERATION " << globalIteration();
      env_cells.removeItems(cells_to_remove);
      to_add_rm_cells = true;
    }
    
  }
#else
  ParallelLoopOptions mtopt;
  mtopt.setPartitioner(ParallelLoopOptions::Partitioner::Static);

  Integer max_nb_task = TaskFactory::nbAllowedThread();
  UniqueArray< Int32UniqueArray > cid_to_add_per_task(max_nb_task);
  UniqueArray< Int32UniqueArray > cid_to_rm_per_task(max_nb_task);

  // Par environnement, on détermine si une maille :
  //  - doit être ajouté à l'env : +1
  //  - doit être retiré de l'env : -1
  //  - ne change pas de status
  const Real threshold = options()->threshold;
  ConstArrayView<IMeshEnvironment*> envs = mm->environments();
  auto queue_arm = m_acc_env->newQueue();

  for (Integer index_env=0; index_env < nb_env ; index_env++) 
  { 
    auto command = makeCommand(queue_arm);
    
    auto in_euler_volume = ax::viewIn(command, m_euler_volume);
    auto in_u_lagrange   = ax::viewIn(command, m_u_lagrange);

    auto out_cell_status = ax::viewOut(command, m_cell_status); // var tempo
 
    // Pour décrire l'accés multi-env sur GPU
    auto in_menv_cell(m_acc_env->multiEnvCellStorage()->viewIn(command));

    command.addKernelName("add_rm") << RUNCOMMAND_ENUMERATE(Cell,cid,allCells())
    {
      Real fvol = in_u_lagrange[cid][index_env] / in_euler_volume[cid];

      // Recherche de l'appartenance de la maille à l'env index_env
      bool cell_in_env = false;
      for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
        auto evi = in_menv_cell.envCell(cid,ienv);
        Integer index_env_loc = in_menv_cell.envId(cid,ienv);
        if (index_env_loc == index_env)
          cell_in_env = true;
      }

      Integer cell_status = 0; // par défaut, le status de la maille ne change pas
      if (fvol < threshold && cell_in_env) {
        cell_status = -1; // maille à retirer de l'env
      } else if (fvol > threshold && !cell_in_env) {
        // verification que le volume normaliseé fournit une fraction de volume au-dessus du threshold             
        Real volt(0.);
        for (Integer index_env_loc=0; index_env_loc < nb_env; index_env_loc++) { 
          // somme des volumes
          volt += in_u_lagrange[cid][index_env_loc];
        }
        // Simplification :
        // (m_u_lagrange[icell][index_env]* m_euler_volume[icell] /volt)/m_euler_volume[icell]
        //  == m_u_lagrange[icell][index_env]/volt
        Real fvol_ev_apres_normalisation = in_u_lagrange[cid][index_env] / volt;
        if (fvol_ev_apres_normalisation > threshold
            && in_u_lagrange[cid][nb_env + index_env] != 0.) {
          cell_status = +1; // maille à ajouter à l'env
        }
      }
      out_cell_status[cid] = cell_status;
    };

    // On revient sur CPU pour lire m_cell_status et créer les listes de mailles à ajouter/supprimer
    // On multithread le traitement
    IMeshEnvironment* env = envs[index_env];

    arcaneParallelForeach(allCells(), mtopt, [&](CellVectorView cells) {
      // Chaque tache met à jour ses listes spécifiques
      Integer task_id = TaskFactory::currentTaskIndex();
      Int32UniqueArray& cid_to_add    = cid_to_add_per_task[task_id];
      Int32UniqueArray& cid_to_remove = cid_to_rm_per_task [task_id];
  
      ENUMERATE_CELL(icell, cells){
        Integer cell_status = m_cell_status[icell];
        if (cell_status < 0) 
        {
          cid_to_remove.add(icell.localId());
          debug() << " cell " << icell.localId() << " ( " << icell->uniqueId() << " ) " << " retirée dans l'env " << env->name();
        } 
        else if (cell_status > 0) 
        {
          cid_to_add.add(icell.localId());
          debug() << " cell " << icell.localId() << " ( " << icell->uniqueId() << " ) " 
            << " ajoutée dans l'env apres normalisation" << env->name();
          debug() << " volume : " <<  m_u_lagrange[icell][index_env] << " fracvol " << m_u_lagrange[icell][index_env] / m_euler_volume[icell];
          /* debug() << " volume-apres_nomalisation " << vol_ev_apres_normalisation <<  " fracvol " << vol_ev_apres_normalisation / m_euler_volume[icell]; */
          debug() << " masse projetée " << m_u_lagrange[icell][nb_env + index_env];
        }
      }
    });

    // En séquentiel, on va ajouter les listes de toutes les taches
    // On précalcule les tailles pour réduire le nb d'alloc dynamiques
    Integer ncells_to_add=0, ncells_to_rm=0; 
    for(Integer itask=0 ; itask<max_nb_task ; ++itask) {
      ncells_to_add += cid_to_add_per_task[itask].size();
      ncells_to_rm  += cid_to_rm_per_task [itask].size();
    }
    if (ncells_to_add) {
      Int32UniqueArray cells_to_add;
      cells_to_add.reserve(ncells_to_add); // but : ne faire qu'une allocation
      for(Integer itask=0 ; itask<max_nb_task ; ++itask) {
        cells_to_add.addRange(cid_to_add_per_task[itask]);
        cid_to_add_per_task[itask].clear(); // Nettoyage pour le prochain environnement
      }
      // On fait un tri pour avoir la même liste quelles que soient les contributions de chaque tache
      std::sort(cells_to_add.data(), cells_to_add.data()+cells_to_add.size());

      // On ajoute réellement les items à l'environnement
      CellGroup env_cells = env->cells();
      pinfo() << "ADD_CELLS to env " << env->name() << " n=" << cells_to_add.size(); 
      env_cells.addItems(cells_to_add);
      to_add_rm_cells = true;
    }
    if (ncells_to_rm) {
      Int32UniqueArray cells_to_remove;
      cells_to_remove.reserve(ncells_to_rm); // but : ne faire qu'une allocation
      for(Integer itask=0 ; itask<max_nb_task ; ++itask) {
        cells_to_remove.addRange(cid_to_rm_per_task[itask]);
        cid_to_rm_per_task[itask].clear(); // Nettoyage pour le prochain environnement
      }
      // On fait un tri pour avoir la même liste quelles que soient les contributions de chaque tache
      std::sort(cells_to_remove.data(), cells_to_remove.data()+cells_to_remove.size());

      // On retire réellement les items de l'environnement
      CellGroup env_cells = env->cells();
      pinfo() << "REMOVE_CELLS to env " << env->name() << " n=" << cells_to_remove.size();
      env_cells.removeItems(cells_to_remove);
      to_add_rm_cells = true;
    }
  }
#endif
  PROF_ACC_END;

  if (to_add_rm_cells) 
  {
    PROF_ACC_BEGIN("forceRecompute");
    // finalisation avant remplissage des variables
    mm->forceRecompute();
    PROF_ACC_END;

    // Ici, la carte des environnements a changé
    m_acc_env->updateMultiEnv(mm);
  }

#if 0
  UniqueArray<Real> vol_nplus1(nb_env);
  UniqueArray<Real> density_env_nplus1(nb_env);
  UniqueArray<Real> internal_energy_env_nplus1(nb_env);
  m_cell_mass.fill(0.0);
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
    /*
    info() << " cell " << cell.localId() << " fin des masses et volumes " << volt;*/
    double volt_normalise = 0.;   
    Real unsurvolt = 1./ volt;
    // normalisation des volumes + somme 
    for (Integer index_env=0; index_env < nb_env ; index_env++) { 
      // vol_nplus1[index_env] *= vol / volt;
      vol_nplus1[index_env] *= vol * unsurvolt;
      volt_normalise += vol_nplus1[index_env];
    }
    // info() << " cell " << cell.localId() << " fin des masses et volumes normalisées ";
    double somme_frac = 0.;
    Real unsurvol = 1. / vol;
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;     
      Integer index_env = ev.environmentId();        
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
    Real unsursomme_frac = 1. / somme_frac;
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();  
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
    
    // on ne recalcule par les mailles à masses nulles - cas advection
    // on enleve les petits fractions de volume aussi sur la fraction
    // massique et on normalise
    Real fmasset = 0.;
    if (masset != 0.) {
      Real unsurmasset = 1./  masset;
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;  
        Integer index_env = ev.environmentId();  
        // m_mass_fraction[ev] = m_u_lagrange[cell][nb_env + index_env] / masset;
        m_mass_fraction[ev] = m_u_lagrange[cell][nb_env + index_env] * unsurmasset;
        if (m_fracvol[ev] < options()->threshold) {
          m_mass_fraction[ev] = 0.;
        }
        fmasset += m_mass_fraction[ev];
      }
      if (fmasset!= 0.) {
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
      Integer index_env = ev.environmentId();  
      density_env_nplus1[index_env] = 0.;
      if (m_fracvol[ev] > options()->threshold)
        density_env_nplus1[index_env] = m_u_lagrange[cell][nb_env + index_env] 
                / vol_nplus1[index_env];
      density_nplus1 += m_fracvol[ev] * density_env_nplus1[index_env];
      // 1/density_nplus1 += m_mass_fraction_env(cCells)[imat] / density_env_nplus1[imat];  
    }
    Real energie_nplus1 = 0.;
    Real pseudo_nplus1 = 0.;
    // Boucle A
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell; 
      Integer index_env = ev.environmentId();  
      internal_energy_env_nplus1[index_env] = 0.;
      
      if (m_fracvol[ev] > options()->threshold && m_u_lagrange[cell][nb_env + index_env] != 0.) {
        internal_energy_env_nplus1[index_env] =
          m_u_lagrange[cell][2 * nb_env + index_env] / m_u_lagrange[cell][nb_env + index_env];
      }
    }
    // mise à jour des valeurs moyennes aux allCells
    // densite
    m_density[cell] = density_nplus1;
    // pinfo() << cell.localId() << " apres proj " << m_u_lagrange[cell];
    // recalcul de la masse
    // Boucle B
    m_cell_mass[cell] = m_euler_volume[cell] * density_nplus1;
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell; 
      Integer index_env = ev.environmentId();  
      m_cell_mass[ev] = m_mass_fraction[ev] * m_cell_mass[cell];
      // recuperation de la pseudo projetee
      // m_pseudo_viscosity[ev] = m_u_lagrange[cell][3 * nb_env + 4] / vol;
      m_pseudo_viscosity[ev] = m_u_lagrange[cell][3 * nb_env + 4] * unsurvol;
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
    }
    // energie interne
    m_internal_energy[cell] = energie_nplus1;
    // pseudoviscosité
    m_pseudo_viscosity[cell] = pseudo_nplus1;
    
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;
      Integer index_env = ev.environmentId();  

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
#else
  auto queue = m_acc_env->newQueue();
  {
    auto command = makeCommand(queue);

    const Real threshold = options()->threshold;

    auto in_euler_volume = ax::viewIn(command, m_euler_volume);
    auto in_u_lagrange   = ax::viewIn(command, m_u_lagrange);
    
    auto out_cell_volume_g      = ax::viewOut(command, m_cell_volume.globalVariable());
    auto out_density_g          = ax::viewOut(command, m_density.globalVariable());
    auto out_pseudo_viscosity_g = ax::viewOut(command, m_pseudo_viscosity.globalVariable());
    auto out_internal_energy_g  = ax::viewOut(command, m_internal_energy.globalVariable());
    auto out_est_pure           = ax::viewOut(command, m_est_pure);

    auto inout_est_mixte     = ax::viewInOut(command, m_est_mixte);
    auto inout_cell_mass_g   = ax::viewInOut(command, m_cell_mass.globalVariable());

    // Variables multi-environnement
    MultiEnvVar<Real> menv_fracvol(m_fracvol, mm);
    auto inout_menv_fracvol(menv_fracvol.span());

    MultiEnvVar<Real> menv_mass_fraction(m_mass_fraction, mm);
    auto inout_menv_mass_fraction(menv_mass_fraction.span());

    MultiEnvVar<Real> menv_cell_mass(m_cell_mass, mm);
    auto inout_menv_cell_mass(menv_cell_mass.span());

    MultiEnvVar<Real> menv_pseudo_viscosity(m_pseudo_viscosity, mm);
    auto inout_menv_pseudo_viscosity(menv_pseudo_viscosity.span());

    MultiEnvVar<Real> menv_density(m_density, mm);
    auto inout_menv_density(menv_density.span());

    MultiEnvVar<Real> menv_internal_energy(m_internal_energy, mm);
    auto inout_menv_internal_energy(menv_internal_energy.span());

    // Pour décrire l'accés multi-env sur GPU
    auto in_menv_cell(m_acc_env->multiEnvCellStorage()->viewIn(command));

    command.addKernelName("moy") << RUNCOMMAND_ENUMERATE(Cell,cid,allCells())
    {
      Real vol = in_euler_volume[cid];  // volume euler   
      out_cell_volume_g[cid] = vol; // retour à la grille euler
      Real volt = 0.;
      Real masset = 0.;

      // info() << " cell " << cell.localId() << " calcul des masses et volumes totales " 
      //        << m_u_lagrange[cell][0] << " et " << m_u_lagrange[cell][1];
      for (Integer index_env=0; index_env < nb_env; index_env++) { 
        Real vol_nplus1 = in_u_lagrange[cid][index_env];
        // somme des volumes
        volt += vol_nplus1;
        // somme des masses
        masset += in_u_lagrange[cid][nb_env + index_env];
      }
      /*
         info() << " cell " << cell.localId() << " fin des masses et volumes " << volt;*/
      double volt_normalise = 0.;   
      Real unsurvolt = 1./ volt;
      // normalisation des volumes + somme 
      for (Integer index_env=0; index_env < nb_env ; index_env++) { 
        // vol_nplus1[index_env] *= vol / volt;
        Real vol_nplus1_norm = in_u_lagrange[cid][index_env] * vol * unsurvolt;
        volt_normalise += vol_nplus1_norm;
      }
      // info() << " cell " << cell.localId() << " fin des masses et volumes normalisées ";
      double somme_frac = 0.;
      Real unsurvol = 1. / vol;
      for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
        auto evi = in_menv_cell.envCell(cid,ienv);
        Integer index_env = in_menv_cell.envId(cid,ienv);

        // Simplification : 
        // in_u_lagrange[cid][index_env] *vol*unsurvolt*unsurvol == in_u_lagrange[cid][index_env] *unsurvolt
        Real fvol = in_u_lagrange[cid][index_env] * unsurvolt;
        if (fvol < threshold)
          fvol = 0.;
        inout_menv_fracvol.setValue(evi, fvol);
        somme_frac += fvol;
      }
      // apres normamisation
      Integer matcell(0);
      Integer imatpure(-1);  
      Real unsursomme_frac = 1. / somme_frac;
      for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
        auto evi = in_menv_cell.envCell(cid,ienv);
        Integer index_env = in_menv_cell.envId(cid,ienv);

        Real fvol = inout_menv_fracvol[evi] * unsursomme_frac;
        if (fvol > 0.) {
          matcell++;
          imatpure = index_env;
        }
        inout_menv_fracvol.setValue(evi, fvol);
      }
      if (matcell > 1) {
        inout_est_mixte[cid] = 1;
        out_est_pure[cid] = -1;
      } else {
        inout_est_mixte[cid] = 0;
        out_est_pure[cid] = imatpure;
      }

      // on ne recalcule par les mailles à masses nulles - cas advection
      // on enleve les petits fractions de volume aussi sur la fraction
      // massique et on normalise
      Real fmasset = 0.;
      if (masset != 0.) {
        Real unsurmasset = 1./  masset;
        for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
          auto evi = in_menv_cell.envCell(cid,ienv);
          Integer index_env = in_menv_cell.envId(cid,ienv);

          // m_mass_fraction[ev] = m_u_lagrange[cell][nb_env + index_env] / masset;
          Real mass_frac = in_u_lagrange[cid][nb_env + index_env] * unsurmasset;
          if (inout_menv_fracvol[evi] < threshold) {
            mass_frac = 0.;
          }
          fmasset += mass_frac;
          inout_menv_mass_fraction.setValue(evi, mass_frac);
        }
        if (fmasset!= 0.) {
          Real unsurfmasset = 1. / fmasset;
          for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
            auto evi = in_menv_cell.envCell(cid,ienv);

            // m_mass_fraction[ev] /= fmasset;
            Real mass_frac = inout_menv_mass_fraction[evi] * unsurfmasset;
            inout_menv_mass_fraction.setValue(evi, mass_frac);
          }
        }
      }
      Real density_nplus1 = 0.;
      for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
        auto evi = in_menv_cell.envCell(cid,ienv);
        Integer index_env = in_menv_cell.envId(cid,ienv);

        Real density_env_nplus1 = 0.;
        Real fvol = inout_menv_fracvol[evi];
        if (fvol > threshold) {
          Real vol_nplus1_norm = in_u_lagrange[cid][index_env] * vol * unsurvolt;
          density_env_nplus1 = in_u_lagrange[cid][nb_env + index_env] 
            / vol_nplus1_norm;
        }
        // recuperation de la densite, se trouvait avant dans boucle B
        inout_menv_density.setValue(evi, density_env_nplus1);

        density_nplus1 += fvol * density_env_nplus1;
        // 1/density_nplus1 += m_mass_fraction_env(cCells)[imat] / density_env_nplus1[imat];  
      }
      // Ici, on fusionne 2 boucles pour ne pas utiliser un tableau internal_energy_env_nplus1[*]
      //
      // mise à jour des valeurs moyennes aux allCells
      // densite
      if (inout_est_mixte[cid])
        out_density_g[cid] = density_nplus1; 
      // density_nplus1 est la valeur moyenne,
      // ne pas affecter dans le cas d'une maille pure sinon on n'aurait pas la valeur partielle

      // pinfo() << cell.localId() << " apres proj " << m_u_lagrange[cell];
      // recalcul de la masse
      // Rem : on n'a plus besoin de m_cell_mass.fill(0.0)
      inout_cell_mass_g[cid] = in_euler_volume[cid] * density_nplus1;

      Real energie_nplus1 = 0.;
      Real pseudo_nplus1 = 0.;
      for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
        auto evi = in_menv_cell.envCell(cid,ienv);
        Integer index_env = in_menv_cell.envId(cid,ienv);

        Real fvol      = inout_menv_fracvol[evi];
        Real mass_frac = inout_menv_mass_fraction[evi];

        // coeur boucle A
        Real internal_energy_env_nplus1 = 0.;
        if (fvol > threshold && in_u_lagrange[cid][nb_env + index_env] != 0.) {
          internal_energy_env_nplus1 =
            in_u_lagrange[cid][2 * nb_env + index_env] / in_u_lagrange[cid][nb_env + index_env];
        }
        // recuperation de l'energie, se trouvait avant dans boucle B
        inout_menv_internal_energy.setValue(evi, internal_energy_env_nplus1);
        // conservation energie totale
        // delta_ec : energie specifique
        // m_internal_energy_env[ev] += delta_ec;
        // energie interne totale, se trouvait avant dans boucle B
        energie_nplus1 += mass_frac * internal_energy_env_nplus1;

        // coeur boucle B
        inout_menv_cell_mass.setValue(evi, mass_frac * inout_cell_mass_g[cid]);
        // recuperation de la pseudo projetee
        // m_pseudo_viscosity[ev] = m_u_lagrange[cell][3 * nb_env + 4] / vol;
        Real pseudo_visco = in_u_lagrange[cid][3 * nb_env + 4] * unsurvol;
        inout_menv_pseudo_viscosity.setValue(evi, pseudo_visco);
        pseudo_nplus1 += fvol * pseudo_visco;
      }
      // energie interne
      out_internal_energy_g[cid] = energie_nplus1;
      // pseudoviscosité
      out_pseudo_viscosity_g[cid] = pseudo_nplus1;

      // Boucle de "blindage" pas totalement portée
      for(Integer ienv=0 ; ienv<in_menv_cell.nbEnv(cid) ; ++ienv) {
        auto evi = in_menv_cell.envCell(cid,ienv);

        if (inout_menv_density[evi] < 0. || inout_menv_internal_energy[evi] < 0.) {
          // Comment gérer des messages d'erreur ou d'avertissement ?
          /*
          pinfo() << " cell " << cell.localId() << " --energy ou masse negative pour l'environnement "
            << index_env;
          pinfo() << " energie interne env " << m_internal_energy[ev] 
            << " cell " << m_internal_energy[cell] ;
          pinfo() << " densite env " << m_density[ev] 
            << " cell " << m_density[cell] ;
          pinfo() << " fraction vol env " << m_fracvol[ev];
          pinfo() << " fraction massique env " <<  m_mass_fraction[ev];
          */
          inout_menv_internal_energy.setValue(evi, 0.);
          inout_menv_density.setValue(evi, 0.);
          inout_menv_fracvol.setValue(evi, 0.);
          inout_menv_mass_fraction.setValue(evi, 0.);
        }
        // NAN sur GPU ? comment gérer l'arrêt du calcul ?
        /*
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
        */
      }
    };
  }
#endif

  auto queue_v = m_acc_env->newQueue();
  if (withDualProjection) {
    // variables aux noeuds
#if 0
    ENUMERATE_NODE(inode, allNodes()){
      Real mass_nodale_proj = m_u_dual_lagrange[inode][3];
      if (mass_nodale_proj != 0.) {
        m_velocity[inode].x = m_u_dual_lagrange[inode][0] / mass_nodale_proj;
        m_velocity[inode].y = m_u_dual_lagrange[inode][1] / mass_nodale_proj;
        m_velocity[inode].z = m_u_dual_lagrange[inode][2] / mass_nodale_proj;
      } 
    }
#else
    auto command_v = makeCommand(queue_v);

    auto in_u_dual_lagrange = ax::viewIn(command_v, m_u_dual_lagrange);
    auto out_velocity       = ax::viewOut(command_v, m_velocity);

    command_v.addKernelName("vel") << RUNCOMMAND_ENUMERATE(Node,nid,allNodes()) {
      Real mass_nodale_proj = in_u_dual_lagrange[nid][3];
      if (mass_nodale_proj != 0.) {
        out_velocity[nid].setX(in_u_dual_lagrange[nid][0] / mass_nodale_proj);
        out_velocity[nid].setY(in_u_dual_lagrange[nid][1] / mass_nodale_proj);
        out_velocity[nid].setZ(in_u_dual_lagrange[nid][2] / mass_nodale_proj);
      }
    };
#endif
  }
  // recalcule de la masse mass_nodale et 
#if 0
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
  m_node_mass.synchronize();
#else
  // On peut calculer simultanément m_node_mass sur les noeuds 
  // et m_cell_volume sur les mailles mixtes
  auto queue_nm = m_acc_env->newQueue();
  queue_nm.setAsync(true);
  {
    auto command_nm = makeCommand(queue_nm);

    auto in_cell_mass  = ax::viewIn(command_nm, m_cell_mass.globalVariable());
    auto out_node_mass = ax::viewOut(command_nm, m_node_mass);

    auto nc_cty = m_acc_env->connectivityView().nodeCell();
    
    Real one_over_nbnode = dimension == 2 ? .25  : .125 ;

    // On a inversé la boucle Cell<->Node pour être parallèle
    command_nm.addKernelName("nodem") << RUNCOMMAND_ENUMERATE(Node,nid,allNodes()) {
      Real contrib_cell_mass=0.;
      for( CellLocalId cid : nc_cty.cells(nid) ){
        contrib_cell_mass += in_cell_mass[cid];
      }
      out_node_mass[nid] = one_over_nbnode * contrib_cell_mass;
    }; // asynchrone et non bloquant vis-à-vis du CPU et des autres queues
  }

  // Mise à jour des m_cell_volume sur les mailles mixtes
  auto menv_queue = m_acc_env->multiEnvQueue();
  ENUMERATE_ENV(ienv, mm) {
    IMeshEnvironment* env = *ienv;

    // Nombre de mailles impures (mixtes) de l'environnement
    Integer nb_imp = env->impureEnvItems().nbItem();

    // Des sortes de vues sur les valeurs impures pour l'environnement env
    Span<Real> out_cell_volume(envView(m_cell_volume, env));
    Span<const Real> in_fracvol(envView(m_fracvol, env));
    Span<const Integer> in_global_cell(envView(m_global_cell, env));

    // Les kernels sont lancés de manière asynchrone environnement par environnement
    auto command = makeCommand(menv_queue->queue(env->id()));

    auto in_cell_volume_g  = ax::viewIn(command,m_cell_volume.globalVariable()); 

    command << RUNCOMMAND_LOOP1(iter,nb_imp) {
      auto [imix] = iter(); // imix \in [0,nb_imp[
      CellLocalId cid(in_global_cell[imix]); // on récupère l'identifiant de la maille globale
      out_cell_volume[imix] = in_fracvol[imix] * in_cell_volume_g[cid];
    };
  }

  queue_nm.barrier();
  m_node_mass.synchronize();
  menv_queue->waitAllQueues();
#endif
 
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
          ec_reconst = 0.5 * m_velocity[inode].abs2();
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
  PROF_ACC_END;
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
