// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "MahycoModule.h"

#include <arcane/MathUtils.h>
#include <arcane/IParallelMng.h>
#include <arcane/ITimeLoopMng.h>

#include <arcane/geometry/IGeometry.h>
#include <arcane/mesh/GhostLayerMng.h>

#include <arcane/AcceleratorRuntimeInitialisationInfo.h>

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/* Constructeur */
/*---------------------------------------------------------------------------*/
MahycoModule::
  MahycoModule(const ModuleBuildInfo& mbi)
: ArcaneMahycoObject(mbi), 
  m_node_index_in_cells(platform::getAcceleratorHostMemoryAllocator())
{}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
accBuild()
{
  PROF_ACC_BEGIN(__FUNCTION__);

  info() << "Using MaHyCo with accelerator";
  IApplication* app = subDomain()->application();
  initializeRunner(m_runner,traceMng(),app->acceleratorRuntimeInitialisationInfo());

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MahycoModule::
_computeNodeIndexInCells() {
  debug() << "_computeNodeIndexInCells";
  // Un noeud est connecté au maximum à MAX_NODE_CELL mailles
  // Calcul pour chaque noeud son index dans chacune des
  // mailles à laquelle il est connecté.
  NodeGroup nodes = allNodes();
  Integer nb_node = nodes.size();
  m_node_index_in_cells.resize(MAX_NODE_CELL*nb_node);
  m_node_index_in_cells.fill(-1);
  auto node_cell_cty = m_connectivity_view.nodeCell();
  auto cell_node_cty = m_connectivity_view.cellNode();
  ENUMERATE_NODE(inode,nodes){
    NodeLocalId node = *inode;
    Int32 index = 0; 
    Int32 first_pos = node.localId() * MAX_NODE_CELL;
    for( CellLocalId cell : node_cell_cty.cells(node) ){
      Int16 node_index_in_cell = 0; 
      for( NodeLocalId cell_node : cell_node_cty.nodes(cell) ){
        if (cell_node==node)
          break;
        ++node_index_in_cell;
      }    
      m_node_index_in_cells[first_pos + index] = node_index_in_cell;
      ++index;
    }    
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
hydroStartInit()
{
  PROF_ACC_BEGIN(__FUNCTION__);

  m_connectivity_view.setMesh(this->mesh());

   IParallelMng* m_parallel_mng = subDomain()->parallelMng();
   my_rank = m_parallel_mng->commRank();
   
  
  pinfo() <<  "Mon rang " << my_rank << " et mon nombre de mailles " << allCells().size();
  pinfo() <<  " Mes mailles pures : " << ownCells().size();
  pinfo() <<  " Mes mailles frantomes : " << allCells().size() - ownCells().size();
  
  info() << " Check donnees ";
   if ((options()->remap()->getOrdreProjection() == 3) && (mesh()->ghostLayerMng()->nbGhostLayer() != 3) && (m_parallel_mng->isParallel() == true)) {
       pinfo() << " mode parallele : " << m_parallel_mng->isParallel();
       pinfo() << " nombre de couches de mailles fantomes : " << mesh()->ghostLayerMng()->nbGhostLayer();
       pinfo() << " incompatible avec la projection d'ordre " << options()->remap()->getOrdreProjection();
       pinfo() << " ----------------------------- fin du calcul à la fin de l'init ---------------------------------------------";
       subDomain()->timeLoopMng()->stopComputeLoop(true);
  }
  if ((options()->withProjection == true) && (mesh()->ghostLayerMng()->nbGhostLayer() < 2) && (m_parallel_mng->isParallel() == true)) {
      pinfo() << " mode parallele : " << m_parallel_mng->isParallel();
      pinfo() << " nombre de couches de mailles fantomes : " << mesh()->ghostLayerMng()->nbGhostLayer();
      pinfo() << " incompatible avec la projection ";
      pinfo() << " ----------------------------- fin du calcul à la fin de l'init ---------------------------------------------";
      subDomain()->timeLoopMng()->stopComputeLoop(true);
  }
  
  m_cartesian_mesh = ICartesianMesh::getReference(mesh());
  m_dimension = mesh()->dimension(); 
  m_cartesian_mesh->computeDirections();
  
  // Dimensionne les variables tableaux
  m_cell_cqs.resize(4*(m_dimension-1));
  m_cell_cqs_n.resize(4*(m_dimension-1));

  // Permet la lecture des cqs quand on boucle sur les noeuds
  _computeNodeIndexInCells();
  
    // Initialise le delta-t
  Real deltat_init = options()->deltatInit();
  m_global_deltat = deltat_init;

  info() << " Initialisation des environnements";
  hydroStartInitEnvAndMat();
 
  // Initialise les données géométriques: volume, cqs, longueurs caractéristiques
  computeGeometricValues(); 
  
  info() << " Initialisation des groupes de faces";
  PrepareFaceGroup();
  
  info() << " Initialisation des variables";
  // Initialises les variables (surcharge l'init d'arcane)
  options()->casModel()->initVar(m_dimension);
  
  if (!options()->sansLagrange) {
    for( Integer i=0,n=options()->environment().size(); i<n; ++i ) {
        IMeshEnvironment* ienv = mm->environments()[i];
        // Initialise l'énergie et la vitesse du son
        options()->environment[i].eosModel()->initEOS(ienv);
    }

    CellToAllEnvCellConverter all_env_cell_converter(mm);

    // valeur moyenne
    ENUMERATE_CELL(icell, allCells()){
        Cell cell = * icell;
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        if (all_env_cell.nbEnvironment() !=1) {
            m_internal_energy[cell] = 0.;
            ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
                EnvCell ev = *ienvcell;        
                m_internal_energy[cell] += m_internal_energy[ev] * m_mass_fraction[ev];
                m_sound_speed[cell] = std::max(m_sound_speed[ev], m_sound_speed[cell]);
            }
        }
    }
  }
  // initialisation du volume Euler
  ENUMERATE_CELL(icell,allCells()){
    Cell cell = *icell;
    // initialisation du volume Euler
    m_euler_volume[cell] = m_cell_volume[cell];
  }
  
  Arcane::Numerics::IGeometryMng* gm = options()->geometry();
  gm->init();
  auto g = gm->geometry();
  bool is_verbose = false;
  if (is_verbose){
    ENUMERATE_CELL(icell,allCells()){
      Cell c = *icell;
      // initialisation du volume Euler
      info() << "Volume = " << g->computeMeasure(c);  
      for (NodeEnumerator inode(c.nodes()); inode.hasNext(); ++inode)
        info() << c.localId() << " avec " << inode.localId();
    }
  }
  
  InitGeometricValues();
  
  // deplacé de hydrocontinueInit
  // calcul des volumes via les cqs, differents deu calcul dans computeGeometricValues ?
  {
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    auto in_node_coord = ax::viewIn(command,m_node_coord);
    auto in_cell_cqs   = ax::viewIn(command,m_cell_cqs);
    auto out_cell_volume_g  = ax::viewOut(command,m_cell_volume.globalVariable()); 

    auto cnc = m_connectivity_view.cellNode();

    // NOTE : on ne peut pas utiliser un membre sur accélérateur (ex : m_dimension), 
    // cela revient à utiliser this->m_dimension avec this pointeur illicite
    Real inv_dim = 1./m_dimension;

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()){

      // Calcule le volume de la maille
      Span<const Real3> cell_cqs = in_cell_cqs[cid];
      Real volume = 0;

      Int64 index=0;
      for( NodeLocalId nid : cnc.nodes(cid) ){
        volume += math::dot(in_node_coord[nid],  cell_cqs[index]);
        ++index;
      }
      volume *= inv_dim;
      out_cell_volume_g[cid] = volume;
    };
  }
 
  // Remplit la structure contenant les informations sur les conditions aux limites
  // Cela permet de garantir avec les accélérateurs qu'on pourra accéder
  // de manière concurrente aux données.
  {
    m_boundary_conditions.clear();
    for (Integer i = 0, nb = options()->boundaryCondition.size(); i < nb; ++i){
      String NomBC = options()->boundaryCondition[i]->surface;
      FaceGroup face_group = mesh()->faceFamily()->findGroup(NomBC);
      Real value = options()->boundaryCondition[i]->value();
      TypesMahyco::eBoundaryCondition type = options()->boundaryCondition[i]->type();

      BoundaryCondition bcn;
      bcn.nodes = face_group.nodeGroup();
      // attention, cette vue sera à reconstruire si bcn.nodes est modifié
      bcn.boundary_nodes = bcn.nodes.view(); 
      bcn.value = value;
      bcn.type = type;
      m_boundary_conditions.add(bcn);
    }
  }

  PROF_ACC_END;
}

/**
 *******************************************************************************
 * \file computeCellMass()
 * \brief Calcul de la masse des mailles
 *
 * \param  m_cell_volume, m_density, m_mass_fraction_env
 * \return m_cell_mass, m_cell_mass_env
 *******************************************************************************
 */
void MahycoModule::
computeCellMass()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans computeCellMass()";
  {
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    auto in_cell_volume_g = ax::viewIn(command, m_cell_volume.globalVariable());
    auto in_density_g     = ax::viewIn(command, m_density.globalVariable());
    auto out_cell_mass_g  = ax::viewOut(command, m_cell_mass.globalVariable());

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()){
      out_cell_mass_g[cid] = in_density_g[cid] * in_cell_volume_g[cid];
    };
  }
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      Cell cell = ev.globalCell();
      m_cell_mass[ev] = m_mass_fraction[ev] * m_cell_mass[cell];
    }
  }
  PROF_ACC_END;
}
/**
 *******************************************************************************
 * \file computeNodeMass()
 * \brief Calcul de la masse nodale
 *
 * \param  m_cell_mass
 * \return m_node_mass
 *******************************************************************************
 */
void MahycoModule::
computeNodeMass() 
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans computeNodeMass()";
   // Initialisation ou reinitialisation de la masse nodale
  {
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    Real one_over_nbnode = m_dimension == 2 ? .25  : .125 ;
    auto nc_cty = m_connectivity_view.nodeCell();

    auto in_cell_mass_g = ax::viewIn(command,m_cell_mass.globalVariable());
    auto out_node_mass = ax::viewOut(command, m_node_mass);

    // Le calcul sur les noeuds fantômes externes n'est pas correct
    // de toute façon (d'où m_node_mass.synchronize) autant ne calculer
    // que sur ownNodes()
    //NodeGroup node_group = allNodes();
    NodeGroup node_group = ownNodes();

    command << RUNCOMMAND_ENUMERATE(Node,nid,node_group) {
      Real sum_mass = 0.;
      for( CellLocalId cid : nc_cty.cells(nid) )
        sum_mass += in_cell_mass_g[cid];
      out_node_mass[nid] = one_over_nbnode * sum_mass;
    };
  }
  m_node_mass.synchronize();
  PROF_ACC_END;
}
/**
 *******************************************************************************
 * \file hydroContinueInit()
 * \brief Initialisation suite à une reprise
 *
 * \param  
 * \return m_nb_env, m_nb_vars_to_project, m_sens_projection
 *******************************************************************************
 */

void MahycoModule::
hydroContinueInit()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  if (subDomain()->isContinue()) {
    
    debug() << " Entree dans hydroContinueInit()";
    // en reprise 
    m_cartesian_mesh = ICartesianMesh::getReference(mesh());
    m_dimension = mesh()->dimension(); 
    m_cartesian_mesh->computeDirections();
    
    mm = IMeshMaterialMng::getReference(defaultMesh());
  
    mm->recreateFromDump();
    m_nb_env = mm->environments().size();
    m_nb_vars_to_project = 3 * m_nb_env + 3 + 1 + 1;
    
    m_global_old_deltat = m_old_deltat;
    // mise a jour nombre iteration 
    m_global_iteration = m_global_iteration() +1;
  }
  PROF_ACC_END;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MahycoModule::
saveValuesAtN()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans saveValuesAtN()";

  // le pas de temps a été mis a jour a la fin dunpas de temps precedent et arcanne met dans m_global_old_deltat ce pas de temps ?
  // donc nous ont remet le bon old pas de temps
  m_global_old_deltat = m_old_deltat;
  
  // synchronisation debut de pas de temps (avec projection nécéssaire ?)
  m_pseudo_viscosity.synchronize();
  m_density.synchronize();
  m_internal_energy.synchronize();
  m_cell_volume.synchronize();
  m_pressure.synchronize();
  m_cell_cqs.synchronize();
  m_velocity.synchronize();
  
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      m_pseudo_viscosity_nmoins1[ev] = m_pseudo_viscosity_n[ev];
      m_pseudo_viscosity_n[ev] = m_pseudo_viscosity[ev];
      m_pressure_n[ev] = m_pressure[ev];
      m_cell_volume_n[ev] = m_cell_volume[ev];
      m_density_n[ev] = m_density[ev];
      m_internal_energy_n[ev] = m_internal_energy[ev];
    }
  }
  auto queue = makeQueue(m_runner);
  // on va recopier de façon asynchrone et concurrente les grandeurs aux mailles et aux noeuds
  queue.setAsync(true); 

  {
    auto command = makeCommand(queue);

    auto in_pseudo_viscosity = ax::viewIn(command,m_pseudo_viscosity.globalVariable());
    auto in_pressure         = ax::viewIn(command,m_pressure.globalVariable());
    auto in_cell_volume      = ax::viewIn(command,m_cell_volume.globalVariable());
    auto in_density          = ax::viewIn(command,m_density.globalVariable());
    auto in_internal_energy  = ax::viewIn(command,m_internal_energy.globalVariable());
    auto in_cell_cqs         = ax::viewIn(command,m_cell_cqs);

    auto inout_pseudo_viscosity_n = ax::viewInOut(command,m_pseudo_viscosity_n.globalVariable());

    auto out_pseudo_viscosity_nmoins1 = ax::viewOut(command,m_pseudo_viscosity_nmoins1.globalVariable());
    auto out_pressure_n         = ax::viewOut(command,m_pressure_n.globalVariable());
    auto out_cell_volume_n      = ax::viewOut(command,m_cell_volume_n.globalVariable());
    auto out_density_n          = ax::viewOut(command,m_density_n.globalVariable());
    auto out_internal_energy_n  = ax::viewOut(command,m_internal_energy_n.globalVariable());
    auto out_cell_cqs_n         = ax::viewInOut(command,m_cell_cqs_n);

    const Integer nb_node_in_cell = m_cell_cqs.arraySize();

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()){
      out_pseudo_viscosity_nmoins1[cid] = inout_pseudo_viscosity_n[cid];
      inout_pseudo_viscosity_n[cid] = in_pseudo_viscosity[cid];
      out_pressure_n[cid] = in_pressure[cid];
      out_cell_volume_n[cid] = in_cell_volume[cid];
      out_density_n[cid] = in_density[cid];
      out_internal_energy_n[cid] = in_internal_energy[cid];

      // anciennement .copy(...)
      for(Integer i=0 ; i<nb_node_in_cell ; ++i) {
        out_cell_cqs_n[cid][i] = in_cell_cqs[cid][i];
      }
    }; // asynchrone
  }

  bool copy_velocity = !options()->sansLagrange;
  bool copy_node_coord = options()->withProjection && options()->remap()->isEuler();

  if (copy_velocity || copy_node_coord) {
    auto command = makeCommand(queue);

    // if (!options()->sansLagrange) m_velocity_n.copy(m_velocity)
    auto in_velocity = ax::viewIn(command, m_velocity);
    auto out_velocity_n = ax::viewOut(command, m_velocity_n);

    // if (options()->withProjection && options()->remap()->isEuler()) m_node_coord.copy(m_node_coord_0)
    auto in_node_coord_0 = ax::viewIn(command, m_node_coord_0);
    auto out_node_coord = ax::viewOut(command, m_node_coord);

    command << RUNCOMMAND_ENUMERATE(Node,nid,allNodes()) {
      if (copy_velocity)
        out_velocity_n[nid] = in_velocity[nid];
      if (copy_node_coord)
        out_node_coord[nid] = in_node_coord_0[nid];
    }; // asynchrone
  }

  queue.barrier();
 
  PROF_ACC_END;
}    
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
computeArtificialViscosity()
{
  if (options()->sansLagrange) return;
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans computeArtificialViscosity()";
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst(env);
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      Cell cell = ev.globalCell();
      m_pseudo_viscosity[ev] = 0.;     
      if (m_div_u[cell] < 0.0) {
        m_pseudo_viscosity[ev] = 1. / m_tau_density[ev]
          * (-0.5 * m_caracteristic_length[cell] * m_sound_speed[cell] * m_div_u[cell]
             + (adiabatic_cst + 1) / 2.0 * m_caracteristic_length[cell] * m_caracteristic_length[cell]
             * m_div_u[cell] * m_div_u[cell]);
      }
    }
  }
  

  // maille mixte
  // moyenne sur la maille
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      m_pseudo_viscosity[cell] = 0.;     
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;        
        m_pseudo_viscosity[cell] += m_pseudo_viscosity[ev] * m_fracvol[ev];
      }
    }
  }
  PROF_ACC_END;
}
/**
 *******************************************************************************
 * \file updateForceAndVelocity()
 * \brief Calcul de la force et de la vitesse 
 *
 * \param  dt, velocity_in
 *         v_pressure, v_pseudo_viscosity, v_cell_cqs
 * \return m_force et v_velocity_out
 *******************************************************************************
 */
void MahycoModule::
updateForceAndVelocity(Real dt,
    const MaterialVariableCellReal& v_pressure,
    const MaterialVariableCellReal& v_pseudo_viscosity,
    const VariableCellArrayReal3& v_cell_cqs,
    const VariableNodeReal3& v_velocity_in,
    VariableNodeReal3& v_velocity_out) 
{
  debug() << " Entree dans updateForceAndVelocity()";
 
#if 0 
  // Remise à zéro du vecteur des forces.
  m_force.fill(Real3::zero());

  // Calcul pour chaque noeud de chaque maille la contribution
  // des forces de pression
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real pressure = v_pressure[icell] + v_pseudo_viscosity[icell];
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode) {
      m_force[inode] += pressure * v_cell_cqs[icell] [inode.index()];
     }
  }

  VariableNodeReal3InView in_force(viewIn(m_force));
  VariableNodeReal3InView in_velocity(viewIn(v_velocity_in));
  VariableNodeRealInView  in_mass(viewIn(m_node_mass));
  VariableNodeReal3OutView out_velocity(viewOut(v_velocity_out));

  // Calcule l'impulsion aux noeuds
  PRAGMA_IVDEP
  ENUMERATE_SIMD_NODE(inode, allNodes()){
    SimdNode snode=*inode;
    out_velocity[snode] = in_velocity[snode] + ( dt / in_mass[snode]) * in_force[snode];;
  }
#else
  {
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    auto in_pressure         = ax::viewIn(command, v_pressure.globalVariable());
    auto in_pseudo_viscosity = ax::viewIn(command, v_pseudo_viscosity.globalVariable());
    auto in_cell_cqs         = ax::viewIn(command, v_cell_cqs);
    // TODO : supprimer m_force, qui ne devient qu'une variable temporaire de travail
    auto out_force           = ax::viewOut(command, m_force);

    auto node_index_in_cells = m_node_index_in_cells.constSpan();
    auto nc_cty = m_connectivity_view.nodeCell();

    auto in_mass      = ax::viewIn(command, m_node_mass);
    auto in_velocity  = ax::viewIn(command, v_velocity_in);
    auto out_velocity = ax::viewOut(command, v_velocity_out);
    
    command << RUNCOMMAND_ENUMERATE(Node,nid,allNodes()) {
      Int32 first_pos = nid.localId() * MAX_NODE_CELL;
      Integer index = 0;
      Real3 node_force = Real3::zero();
      for( CellLocalId cid : nc_cty.cells(nid) ){
        Int16 node_index = node_index_in_cells[first_pos + index];
        node_force += (in_pressure[cid]+in_pseudo_viscosity[cid]) 
          * in_cell_cqs[cid][node_index];
        ++index;
      }
      out_force[nid] = node_force;

      // On peut mettre la vitesse à jour dans la foulée
      out_velocity[nid] = in_velocity[nid] + ( dt / in_mass[nid]) * node_force;
    };
  }
#endif

  v_velocity_out.synchronize();
}

/**
 *******************************************************************************
 * \file updateVelocity()
 * \brief Calcul de la vitesse de n-1/2 a n+1/2
 *
 * \param  m_global_old_deltat(), m_global_deltat(), m_velocity_n
 *         m_pressure_n, m_pseudo_viscosity_n, m_cqs_n
 * \return m_velocity
 *******************************************************************************
 */

void MahycoModule::
updateVelocity()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  if (options()->sansLagrange) {
    updateVelocityWithoutLagrange();
    PROF_ACC_END;
    return;
  }
  
  // passage des vitesse de n à n-1/2 si projection 
  // la vitesse m_velocity_n 
  // (qui dans le cas de vnr-csts est la vitesse apres projection du pas de temps précédent donc n),
  // est replacee à n-1/2 pour vnr-csts.
  // Dans le cas de vnr (pas csts) elle est deja en n-1/2
  if (options()->schemaCsts() && options()->withProjection) updateVelocityBackward();
      
      
  debug() << " Entree dans updateVelocity()";
#if 0
  // Remise à zéro du vecteur des forces.
  m_force.fill(Real3::zero());

  // Calcul pour chaque noeud de chaque maille la contribution
  // des forces de pression
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real pressure = m_pressure_n[icell] + m_pseudo_viscosity_n[icell];
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode) {
      m_force[inode] += pressure * m_cell_cqs_n[icell] [inode.index()];
     }
  }
  
  VariableNodeReal3InView in_force(viewIn(m_force));
  VariableNodeReal3InView in_velocity(viewIn(m_velocity_n));
  VariableNodeRealInView  in_mass(viewIn(m_node_mass));
  VariableNodeReal3OutView out_velocity(viewOut(m_velocity));

  const Real dt(0.5 * (m_global_old_deltat() + m_global_deltat()));
  // Calcule l'impulsion aux noeuds
  PRAGMA_IVDEP
  ENUMERATE_SIMD_NODE(inode, allNodes()){
    SimdNode snode=*inode;
    out_velocity[snode] = in_velocity[snode] + ( dt / in_mass[snode]) * in_force[snode];;
  }

  m_velocity.synchronize();
#else
  updateForceAndVelocity(0.5 * (m_global_old_deltat() + m_global_deltat()),
        /* calcul m_force : */ m_pressure_n, m_pseudo_viscosity_n, m_cell_cqs_n,
        /* calcul m_velocity : */ m_velocity_n, m_velocity);
#endif
  PROF_ACC_END;
}
/**
 *******************************************************************************
 * \file updateVelocityBackward()
 * \brief Calcul de la vitesse de n a n-1/2
 *
 * \param  gt->deltat_n, m_velocity_n
 *         m_pressure_n, m_pseudo_viscosity_n, m_cqs_n
 * \return m_velocity_n
 *******************************************************************************
 */
void MahycoModule::
updateVelocityBackward()
{
  if (options()->sansLagrange) return;
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans updateVelocityBackward()";
#if 0
  // Remise à zéro du vecteur des forces.
  m_force.fill(Real3::null());

  // Calcul pour chaque noeud de chaque maille la contribution
  // des forces de pression
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real pressure = m_pressure_n[icell] + m_pseudo_viscosity_n[icell];
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      m_force[inode] += pressure * m_cell_cqs_n[icell] [inode.index()];
  }
  const Real dt(-0.5 * m_global_old_deltat());  
  
  VariableNodeReal3InView in_force(viewIn(m_force));
  VariableNodeReal3InView in_velocity(viewIn(m_velocity_n));
  VariableNodeRealInView  in_mass(viewIn(m_node_mass));
  VariableNodeReal3OutView out_velocity(viewOut(m_velocity_n));

  // Calcule l'impulsion aux noeuds
  PRAGMA_IVDEP
  ENUMERATE_SIMD_NODE(inode, allNodes()){
    SimdNode snode=*inode;
    out_velocity[snode] = in_velocity[snode] + ( dt / in_mass[snode]) * in_force[snode];;
  }

  m_velocity_n.synchronize();
#else
  updateForceAndVelocity(-0.5 * m_global_old_deltat(),
      /* calcul m_force : */ m_pressure_n, m_pseudo_viscosity_n, m_cell_cqs_n,
      /* calcul m_velocity_n : */ m_velocity_n, m_velocity_n);
#endif
  PROF_ACC_END;
}
/*******************************************************************************
 * \file updateVelocityForward()
 * \brief Calcul de la vitesse de n+1/2 a n+1
 *
 * \param  m_global_deltat, m_velocity
 *         m_pressure, m_pseudo_viscosity, m_cqs (calcule juste avant)
 * \return m_velocity
 *******************************************************************************
 */
void MahycoModule::
updateVelocityForward()
{
  if (options()->sansLagrange) return;
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans updateVelocityForward()";
#if 0
  // Remise à zéro du vecteur des forces.
  m_force.fill(Real3::null());

  // Calcul pour chaque noeud de chaque maille la contribution
  // des forces de pression
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real pressure = m_pressure[icell] + m_pseudo_viscosity[icell];
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      m_force[inode] += pressure * m_cell_cqs[icell] [inode.index()];
  }
  const Real dt(0.5 * m_global_deltat());
  
  VariableNodeReal3InView in_force(viewIn(m_force));
  VariableNodeReal3InView in_velocity(viewIn(m_velocity));
  VariableNodeRealInView  in_mass(viewIn(m_node_mass));
  VariableNodeReal3OutView out_velocity(viewOut(m_velocity));

  // Calcule l'impulsion aux noeuds
  PRAGMA_IVDEP
  ENUMERATE_SIMD_NODE(inode, allNodes()){
    SimdNode snode=*inode;
    out_velocity[snode] = in_velocity[snode] + ( dt / in_mass[snode]) * in_force[snode];;
  }
  m_velocity.synchronize();
#else
  updateForceAndVelocity(0.5 * m_global_deltat(),
      /* calcul m_force : */ m_pressure, m_pseudo_viscosity, m_cell_cqs,
      /* calcul m_velocity : */ m_velocity, m_velocity);
#endif
  PROF_ACC_END;
}
/**
 *******************************************************************************
 * \file updateVelocityWithoutLagrange()
 * \brief Calcul de la vitesse pour les cas d'advection pure
 *
 * \param  gt->t_nplus1, m_node_velocity_n0, m_node_velocity_n, m_node_coord_n
 * \return m_node_velocity_nplus1
 *******************************************************************************
 */
void MahycoModule::
updateVelocityWithoutLagrange()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  bool option(options()->casModel()->hasReverseOption());
  Real factor(options()->casModel()->getReverseParameter());
  
  Real option_real( (Real) option);
  
  VariableNodeReal3InView in_velocity(viewIn(m_velocity_n));
  VariableNodeReal3OutView out_velocity(viewOut(m_velocity));

  PRAGMA_IVDEP
  ENUMERATE_SIMD_NODE(inode, allNodes()){
    SimdNode snode=*inode;
    out_velocity[snode] = in_velocity[snode] * (1. -option_real)
          + option_real * in_velocity[snode] * cos(Pi * m_global_time() * factor);
  }

  m_velocity.synchronize();
  PROF_ACC_END;
}
/**
 *******************************************************************************
 * \file updatePosition()
 * \brief Calcul de la position
 *
 * \param  gt->t_nplus1, m_node_velocity_nplus1, m_node_coord_n
 * \return m_node_coord_nplus1
 *******************************************************************************
 */
void MahycoModule::
updatePosition()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans updatePosition()";
  Real deltat = m_global_deltat();
#if 0
  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    if (((options()->sansLagrange) && (node.nbCell() == 4)) || (!options()->sansLagrange))
        m_node_coord[inode] += deltat * m_velocity[inode];
  }
  Real one_over_nbnode = m_dimension == 2 ? .25  : .125 ;
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real3 somme = {0. , 0. , 0.};
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      somme += m_node_coord[inode];
    m_cell_coord[cell] = one_over_nbnode * somme;
  }
#else
  auto queue = makeQueue(m_runner);
  if (!options()->sansLagrange)
  {
    auto command = makeCommand(queue);

    auto in_velocity    = ax::viewIn(command,m_velocity);
    auto out_node_coord = ax::viewOut(command,m_node_coord);

    command << RUNCOMMAND_ENUMERATE(Node,nid,allNodes()) {
      out_node_coord[nid] += deltat * in_velocity[nid];
    };
  }
  else // options()->sansLagrange == true
  {
    ENUMERATE_NODE(inode, allNodes()){
      Node node = *inode;
      if (node.nbCell() == 4) {
        m_node_coord[inode] += deltat * m_velocity[inode];
      }
    }
  }
  {
    Real one_over_nbnode = m_dimension == 2 ? .25  : .125 ;

    auto command = makeCommand(queue);

    auto in_node_coord  = ax::viewIn(command,m_node_coord);
    auto out_cell_coord = ax::viewOut(command,m_cell_coord);
    auto cnc = m_connectivity_view.cellNode();

    command << RUNCOMMAND_ENUMERATE(Cell,cid,allCells()) {
      Real3 somme = {0. , 0. , 0.};
      for( NodeLocalId nid : cnc.nodes(cid) ){
        somme += in_node_coord[nid];
      }
      out_cell_coord[cid] = one_over_nbnode * somme;
    };
  }
#endif
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
applyBoundaryCondition()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans applyBoundaryCondition()";
#if 0
  for (Integer i = 0, nb = options()->boundaryCondition.size(); i < nb; ++i){
    String NomBC = options()->boundaryCondition[i]->surface;
    FaceGroup face_group = mesh()->faceFamily()->findGroup(NomBC);
    Real value = options()->boundaryCondition[i]->value();
    TypesMahyco::eBoundaryCondition type = options()->boundaryCondition[i]->type();

    // boucle sur les faces de la surface
    ENUMERATE_FACE(j, face_group){
      Face face = * j;
      Integer nb_node = face.nbNode();

      // boucle sur les noeuds de la face
      for (Integer k = 0; k < nb_node; ++k){
        Node node = face.node(k);
        Real3& velocity = m_velocity[node];

        switch (type){
        case TypesMahyco::VelocityX:
          velocity.x = value;
          break;
        case TypesMahyco::VelocityY:
          velocity.y = value;
          break;
        case TypesMahyco::VelocityZ:
          velocity.z = value;
          break;
        case TypesMahyco::Unknown:
          break;
        }
      }
    }
  }
#else
  auto queue = makeQueue(m_runner);

  // Pour cette méthode, comme les conditions aux limites sont sur des groupes
  // indépendants (ou alors avec la même valeur si c'est sur les mêmes noeuds),
  // on peut exécuter les noyaux en asynchrone.
  queue.setAsync(true);

  for( auto bc : m_boundary_conditions ) {
    Real value = bc.value;
    TypesMahyco::eBoundaryCondition type = bc.type;

    auto command = makeCommand(queue);
    auto inout_velocity = ax::viewInOut(command,m_velocity);

    command << RUNCOMMAND_ENUMERATE(Node,nid,bc.boundary_nodes) {
      Real3 velocity = inout_velocity[nid];
      switch(type) {
        case TypesMahyco::VelocityX: velocity.x = value; break;
        case TypesMahyco::VelocityY: velocity.y = value; break;
        case TypesMahyco::VelocityZ: velocity.z = value; break;
        case TypesMahyco::Unknown: break;
      }
      inout_velocity[nid] = velocity;
    }; // non-bloquant
  }

  queue.barrier();
#endif
  PROF_ACC_END;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
InitGeometricValues()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans InitGeometricValues() ";
  ENUMERATE_NODE(inode, allNodes()){
      m_node_coord_0[inode] = m_node_coord[inode];
  }
  Real one_over_nbnode = m_dimension == 2 ? .5  : .25 ;
  if ( m_dimension == 3) {
    ENUMERATE_FACE (iFace, allFaces()) {
        Face face = *iFace;
        Real3 face_vec1 = m_node_coord[face.node(2)] - m_node_coord[face.node(0)]; 
        Real3 face_vec2 = m_node_coord[face.node(3)] - m_node_coord[face.node(1)];
        m_face_normal[iFace].x = produit(face_vec1.y, face_vec2.z, face_vec1.z, face_vec2.y);
        m_face_normal[iFace].y = - produit(face_vec2.x, face_vec1.z, face_vec2.z, face_vec1.x);
        m_face_normal[iFace].z = produit(face_vec1.x, face_vec2.y, face_vec1.y, face_vec2.x);
        m_face_normal[iFace] /= m_face_normal[iFace].abs();
    }
  } else {
    ENUMERATE_FACE (iFace, allFaces()) {
      Face face = *iFace;
      m_face_normal[iFace].x = (m_node_coord[face.node(1)].y - m_node_coord[face.node(0)].y); 
      m_face_normal[iFace].y = (m_node_coord[face.node(1)].x - m_node_coord[face.node(0)].x); 
      m_face_normal[iFace] /= m_face_normal[iFace].abs();
    }
  }
  ENUMERATE_FACE (iFace, allFaces()) {
      Face face = *iFace;
      m_face_coord[face] = 0.;
      for (Integer inode = 0; inode < face.nbNode(); ++inode) 
        m_face_coord[face] +=  one_over_nbnode * m_node_coord[face.node(inode)];
  }
  ENUMERATE_CELL(icell,allCells()) {
    ENUMERATE_FACE(iface, (*icell).faces()){
      const Face& face = *iface;
      Integer index = iface.index(); 
        m_outer_face_normal[icell][index] = (m_face_coord[face]-m_cell_coord[icell]) 
            / (m_face_coord[face]-m_cell_coord[icell]).abs();
    }
  }
  PROF_ACC_END;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
computeGeometricValues()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << my_rank << " : " << " Entree dans computeGeometricValues() ";
  
  m_node_coord.synchronize();
  if ( m_dimension == 3) {
    {
      auto queue = makeQueue(m_runner);
      auto command = makeCommand(queue);

      auto in_node_coord = ax::viewIn(command,m_node_coord);
      auto out_cell_cqs = ax::viewInOut(command,m_cell_cqs);

      auto cnc = m_connectivity_view.cellNode();

      command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()){
        // Recopie les coordonnées locales (pour le cache)
        Real3 coord[8];
        Int32 index=0;
        for( NodeLocalId nid : cnc.nodes(cid) ){
          coord[index]=in_node_coord[nid];
          ++index;
        }
        // Calcul les coordonnées des centres des faces
        Real3 face_coord[6] = {
          0.25 * (coord[0] + coord[3] + coord[2] + coord[1]),
          0.25 * (coord[0] + coord[4] + coord[7] + coord[3]),
          0.25 * (coord[0] + coord[1] + coord[5] + coord[4]),
          0.25 * (coord[4] + coord[5] + coord[6] + coord[7]),
          0.25 * (coord[1] + coord[2] + coord[6] + coord[5]),
          0.25 * (coord[2] + coord[3] + coord[7] + coord[6])
        };

        // Calcule les résultantes aux sommets
        computeCQs(coord, face_coord, out_cell_cqs[cid]);
      };
    }
  } else {
    Real3 npc[5];
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;
      // Recopie les coordonnées locales (pour le cache)
      Real3 coord[8];
      for (NodeEnumerator inode(cell.nodes()); inode.index() < cell.nbNode(); ++inode) {
        coord[inode.index()] = m_node_coord[inode];
      }
      coord[4] = coord[0];
      for (NodeEnumerator inode(cell.nodes()); inode.index() < cell.nbNode(); ++inode) {
        npc[inode.index()+1].x = 0.5 * (coord[inode.index()+1].y -  coord[inode.index()].y);
        npc[inode.index()+1].y = 0.5 * (coord[inode.index()].x -  coord[inode.index()+1].x);
        // npc[inode.index()+1] = npc[inode.index()+1] / npc[inode.index()+1].abs();
      }
      npc[0] = npc[4];
      for (Integer ii = 0; ii < 4; ++ii) {
        m_cell_cqs[icell][ii] = npc[ii+1] + npc[ii]; 
      }
    }  
 }
  m_cell_cqs.synchronize(); // TODO : pourquoi synchronize ?
 
  if (options()->longueurCaracteristique() == "faces-opposees")
  {
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;    
      // Calcule le volume de la maille
      {
        Real volume = 0.;

        for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
          volume += math::dot(m_node_coord[cell.node(inode)], m_cell_cqs[icell] [inode]);
          // pinfo() << cell.localId() << " coor " << m_node_coord[cell.node(inode)] << " et " << m_cell_cqs[icell] [inode];
           
          // TODO : bien s'assurer que m_node_volume ne sert à rien
          // m_node_volume[cell.node(inode)] += volume;
        }
        volume /= m_dimension;

        m_cell_volume[cell] = volume;
      }
      // Calcule la longueur caractéristique de la maille.
      {
        // Recopie les coordonnées locales (pour le cache)
        Real3 coord[8];
        for (NodeEnumerator inode(cell.nodes()); inode.index() < 8; ++inode) {
          coord[inode.index()] = m_node_coord[inode];
        }
        // Calcul les coordonnées des centres des faces
        Real3 face_coord[6];
        face_coord[0] = 0.25 * (coord[0] + coord[3] + coord[2] + coord[1]);
        face_coord[1] = 0.25 * (coord[0] + coord[4] + coord[7] + coord[3]);
        face_coord[2] = 0.25 * (coord[0] + coord[1] + coord[5] + coord[4]);
        face_coord[3] = 0.25 * (coord[4] + coord[5] + coord[6] + coord[7]);
        face_coord[4] = 0.25 * (coord[1] + coord[2] + coord[6] + coord[5]);
        face_coord[5] = 0.25 * (coord[2] + coord[3] + coord[7] + coord[6]);

        Real3 median1 = face_coord[0] - face_coord[3];
        Real3 median2 = face_coord[2] - face_coord[5];
        Real3 median3 = face_coord[1] - face_coord[4];
        Real d1 = median1.abs();
        Real d2 = median2.abs();
        Real d3 = median3.abs();

        Real dx_numerator = d1 * d2 * d3;
        Real dx_denominator = d1 * d2 + d1 * d3 + d2 * d3;
        m_caracteristic_length[icell] = dx_numerator / dx_denominator;
      }
    } 
  }
  else if (options()->longueurCaracteristique() == "racine-cubique-volume")
  {
    Real racine = m_dimension == 2 ? .5  : 1./3. ;
    // Calcul des volumes aux mailles puis longueur caractéristique
    // Attention, m_node_volume n'est plus calculé car n'est pas utilisé
    {
      auto queue = makeQueue(m_runner);
      auto command = makeCommand(queue);

      auto in_node_coord = ax::viewIn(command,m_node_coord);
      auto in_cell_cqs   = ax::viewIn(command,m_cell_cqs);
      auto out_cell_volume_g        = ax::viewOut(command,m_cell_volume.globalVariable()); 
      auto out_caracteristic_length = ax::viewOut(command,m_caracteristic_length);
      
      auto cnc = m_connectivity_view.cellNode();

      // NOTE : on ne peut pas utiliser un membre sur accélérateur (ex : m_dimension), 
      // cela revient à utiliser this->m_dimension avec this pointeur illicite
      Real inv_dim = 1./m_dimension;

      command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()){

        // Calcule le volume de la maille
        Span<const Real3> cell_cqs = in_cell_cqs[cid];
        Real volume = 0;

        Int64 index=0;
        for( NodeLocalId nid : cnc.nodes(cid) ){
          volume += math::dot(in_node_coord[nid],  cell_cqs[index]);
          ++index;
        }
        volume *= inv_dim;
        out_cell_volume_g[cid] = volume;

        // Calcule la longueur caractéristique de la maille.
        out_caracteristic_length[cid] = math::pow(volume, racine);
      };
    }
  }
  else
  {
    info() << " pas de longeur caractéritique definie dans le .arc " << options()->longueurCaracteristique(); 
    subDomain()->timeLoopMng()->stopComputeLoop(true);
  }

  // maille mixte
  // moyenne sur la maille
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;        
        Cell cell = ev.globalCell();
        m_cell_volume[ev] = m_fracvol[ev] * m_cell_volume[cell];
      }
    }
  }
  PROF_ACC_END;
}
/**
 *******************************************************************************
 * \file updateDensity()
 * \brief Calcul de la densité
 * \brief Calcul de la divergence de la vitesse
 * \brief Calcul de la variation du volume specifique (variation de 1/densite)
 *
 * pour les grandeurs moyennes et pour chaque materiau
 *
 * \param  m_density_n, ,m_cell_mass, m_cell_volume, m_global_deltat()
 * \return m_density, m_tau_density, m_divu
 *******************************************************************************
 */
void MahycoModule::
updateDensity()
{
  if (options()->sansLagrange) return;
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << my_rank << " : " << " Entree dans updateDensity() ";
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      Cell cell = ev.globalCell();
       // pinfo() << my_rank << " : " << cell.uniqueId() << " " << m_cell_volume[ev];
       Real new_density = m_cell_mass[ev] / m_cell_volume[ev];
       // nouvelle density
       m_density[ev] = new_density;
       // volume specifique de l'environnement au temps n+1/2
       m_tau_density[ev] = 
        0.5 * (1.0 / m_density_n[ev] + 1.0 / m_density[ev]);
        
    }
  }
  {
    auto queue = makeQueue(m_runner);
    auto command = makeCommand(queue);

    Real inv_deltat = 1.0/m_global_deltat(); // ne pas appeler de méthodes de this dans le kernel

    auto in_cell_mass_g   = ax::viewIn(command, m_cell_mass.globalVariable());
    auto in_cell_volume_g = ax::viewIn(command, m_cell_volume.globalVariable());
    auto in_density_n_g   = ax::viewIn(command, m_density_n.globalVariable());

    auto iou_density_g     = ax::viewInOut(command, m_density.globalVariable());
    auto iou_tau_density_g = ax::viewInOut(command, m_tau_density.globalVariable());

    auto out_div_u = ax::viewOut(command, m_div_u);

    command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()){
      Real new_density = in_cell_mass_g[cid] / in_cell_volume_g[cid];
      // nouvelle density
      iou_density_g[cid] = new_density;
      // volume specifique moyen au temps n+1/2
      iou_tau_density_g[cid] = 
        0.5 * (1.0 / in_density_n_g[cid] + 1.0 / iou_density_g[cid]);
      // divergence de la vitesse mode A1
      out_div_u[cid] =
        inv_deltat  * ( 1.0 / iou_density_g[cid] - 1.0 / in_density_n_g[cid] )
        / iou_tau_density_g[cid];
    };
  }
  
  m_density.synchronize();
  m_tau_density.synchronize();
  m_div_u.synchronize();
  PROF_ACC_END;
}
/**
 *******************************************************************************
 * \file updateEnergy()
 * \brief Calcul de l'energie interne ( cas du gaz parfait ou methode de newton)
 *
 * \param  m_density_env_nplus1, m_density_env_n,
 *         m_pseudo_viscosity_env_nplus1, m_pseudo_viscosity_env_n
 *         m_pressure_env_n, m_mass_fraction_env
 *
 * \return m_internal_energy_env_nplus1, m_internal_energy_nplus1
 *******************************************************************************
 */
void MahycoModule::
updateEnergyAndPressure()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  if (options()->withNewton) 
    updateEnergyAndPressurebyNewton();
  else
    updateEnergyAndPressureforGP();  
    
  // maille mixte
  // moyenne sur la maille
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      m_internal_energy[cell] = 0.;
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;
        m_internal_energy[cell] += m_mass_fraction[ev] * m_internal_energy[ev];
      }
    }
  }
  if (! options()->withProjection) {
    // Calcul de la Pression si on ne fait pas de projection 
    for( Integer i=0,n=options()->environment().size(); i<n; ++i ) {
        IMeshEnvironment* ienv = mm->environments()[i];
        // Calcul de la pression et de la vitesse du son
        options()->environment[i].eosModel()->applyEOS(ienv);
    }
    computePressionMoyenne();
  }
  PROF_ACC_END;
}
/*
 *******************************************************************************
*/
void MahycoModule::updateEnergyAndPressurebyNewton()  {  
    
  if (options()->sansLagrange) return;
  PROF_ACC_BEGIN(__FUNCTION__);
    debug() << " Entree dans updateEnergyAndPressure()";
    bool csts = options()->schemaCsts();
    bool pseudo_centree = options()->pseudoCentree();
    // Calcul de l'énergie interne
    if (!csts) {
      ENUMERATE_ENV(ienv,mm){
        IMeshEnvironment* env = *ienv;
        Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst(env);
        Real tension_limit = options()->environment[env->id()].eosModel()->getTensionLimitCst(env);
        ENUMERATE_ENVCELL(ienvcell,env){
          EnvCell ev = *ienvcell;
          Real pseudo(0.);
          if (pseudo_centree &&
            ((m_pseudo_viscosity_n[ev] + m_pseudo_viscosity[ev]) * (1.0 / m_density[ev] - 1.0 / m_density_n[ev]) < 0.))
          pseudo = 0.5 * (m_pseudo_viscosity[ev] + m_pseudo_viscosity_n[ev]);
          if (!pseudo_centree &&
                (m_pseudo_viscosity[ev] * (1.0 / m_density[ev] - 1.0 / m_density_n[ev]) < 0.))
            pseudo = m_pseudo_viscosity[ev];
            
          double rn  = m_density_n[ev];
          double pn  = m_pressure_n[ev];
          double qnn1 = pseudo;
          double m   = m_cell_mass[ev];
          double rn1 = m_density[ev];
          double en  = m_internal_energy_n[ev];
          double g = adiabatic_cst;
          double t = tension_limit;
          // les iterations denewton
          double epsilon = options()->threshold;
          double itermax = 50;
          double enew=0, e=en, p, c, dpde;
          int i = 0;
        
          while(i<itermax && abs(fvnr(e, p, dpde, en, qnn1, pn, rn1, rn))>=epsilon)
            {
              m_internal_energy[ev] = e;
              options()->environment[env->id()].eosModel()->applyOneCellEOS(env, ev);
              p = m_pressure[ev];
              c = m_sound_speed[ev];
              dpde = m_dpde[ev];
              enew = e - fvnr(e, p, dpde, en, qnn1, pn, rn1, rn) / fvnrderiv(e, dpde, rn1, rn);
              e = enew;
              i++;
            }
          m_internal_energy[ev] = e;
	      m_sound_speed[ev] = c;
	      m_pressure[ev] = p;
        }
      }
    } else {
      ENUMERATE_ENV(ienv,mm){
        IMeshEnvironment* env = *ienv;
        Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst(env);
        Real tension_limit = options()->environment[env->id()].eosModel()->getTensionLimitCst(env);
        ENUMERATE_ENVCELL(ienvcell,env){
          EnvCell ev = *ienvcell;
          Cell cell=ev.globalCell();
          Real cqs_v_nplus1(0.);
          Real cqs_v_n(0.);
          Real cqs_delta_v(0.);
          Real cqs_v_old_n(0.);
          for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
            cqs_v_nplus1 += math::dot(m_velocity[cell.node(inode)], m_cell_cqs[cell] [inode])
              * m_global_deltat();
            cqs_v_n += math::dot(m_velocity[cell.node(inode)], m_cell_cqs_n[cell] [inode])
              * m_global_deltat();
            cqs_delta_v +=  math::dot(m_velocity[cell.node(inode)] - m_velocity_n[cell.node(inode)], m_cell_cqs_n[cell] [inode])
            * (m_global_deltat() - m_global_old_deltat());
            cqs_v_old_n += math::dot(m_velocity_n[cell.node(inode)], m_cell_cqs_n[cell] [inode])
            * m_global_old_deltat();
          }
          double rn  = m_density_n[ev];
          double pn  = m_pressure_n[ev];
          double qn = m_pseudo_viscosity_n[ev];
          double qn1 = m_pseudo_viscosity[ev];
          double m   = m_cell_mass[ev];
          double rn1 = m_density[ev];
          double en  = m_internal_energy_n[ev];
          double g = adiabatic_cst;
          double t = tension_limit;
          double cn1 = cqs_v_nplus1;
          double cn = cqs_v_n;
          double cdn = cqs_delta_v;
          double cdon = 0. ;
          double qnm1 = m_pseudo_viscosity_nmoins1[ev]; 
          // les iterations de newton
          double epsilon = options()->threshold;
          double itermax = 50;
          double enew=0, e=en, p, c, dpde;
          int i = 0;
        
          while(i<itermax && abs(f(e, p, dpde, en, qn, pn, cn1, cn, m, qn1, cdn, cdon, qnm1))>=epsilon)
	        {
              m_internal_energy[ev] = e;
              options()->environment[env->id()].eosModel()->applyOneCellEOS(env, ev);
              p = m_pressure[ev];
              c = m_sound_speed[ev];
              dpde = m_dpde[ev];
              enew = e - f(e, p, dpde, en, qn, pn, cn1, cn, m, qn1, cdn, cdon, qnm1) / fderiv(e, p, dpde, cn1, m);
              e = enew;
              i++;
            }
          m_internal_energy[ev] = e;
	      m_sound_speed[ev] = c;
	      m_pressure[ev] = p;
        }
      }
    }
  PROF_ACC_END;
}

/*
 *******************************************************************************
 */
void MahycoModule::
updateEnergyAndPressureforGP()
{
  if (options()->sansLagrange) return;
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans updateEnergyAndPressure()";
  bool csts = options()->schemaCsts();
  bool pseudo_centree = options()->pseudoCentree();
  // Calcul de l'énergie interne
  if (!csts) {
    ENUMERATE_ENV(ienv,mm){
      IMeshEnvironment* env = *ienv;
      Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst(env);
      ENUMERATE_ENVCELL(ienvcell,env){
        EnvCell ev = *ienvcell;
        Real pseudo(0.);
        Cell cell = ev.globalCell();
        if (pseudo_centree &&
            ((m_pseudo_viscosity_n[ev] + m_pseudo_viscosity[ev]) * (1.0 / m_density[ev] - 1.0 / m_density_n[ev]) < 0.))
          pseudo = 0.5 * (m_pseudo_viscosity[ev] + m_pseudo_viscosity_n[ev]);
        if (!pseudo_centree &&
            (m_pseudo_viscosity[ev] * (1.0 / m_density[ev] - 1.0 / m_density_n[ev]) < 0.))
          pseudo = m_pseudo_viscosity[ev];
          
        Real denom_accrois_nrj(1 + 0.5 * (adiabatic_cst - 1.0) *
                               m_density[ev] *
                               (1.0 / m_density[ev] -
                                1.0 / m_density_n[ev]));
        Real numer_accrois_nrj(m_internal_energy_n[ev] -
                               (0.5 * m_pressure[ev] + pseudo) *
                               (1.0 / m_density[ev] - 1.0 / m_density_n[ev]));
        m_internal_energy[ev] = numer_accrois_nrj / denom_accrois_nrj;
      }
    }
  } else {
    ENUMERATE_ENV(ienv,mm){
      IMeshEnvironment* env = *ienv;
      Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst(env);
      ENUMERATE_ENVCELL(ienvcell,env){
        EnvCell ev = *ienvcell;
        Cell cell = ev.globalCell();
        Real cqs_v_nplus1(0.);
        Real cqs_v_n(0.);
        Real cqs_delta_v(0.);
        Real cqs_v_old_n(0.);
        for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
          cqs_v_nplus1 += math::dot(m_velocity[cell.node(inode)], m_cell_cqs[cell] [inode])
            * m_global_deltat();
          cqs_v_n += math::dot(m_velocity[cell.node(inode)], m_cell_cqs_n[cell] [inode])
            * m_global_deltat();
          cqs_delta_v +=  math::dot(m_velocity[cell.node(inode)] - m_velocity_n[cell.node(inode)], m_cell_cqs_n[cell] [inode])
            * (m_global_deltat() - m_global_old_deltat());
          cqs_v_old_n += math::dot(m_velocity_n[cell.node(inode)], m_cell_cqs_n[cell] [inode])
            * m_global_old_deltat();
        }
        Real denom_accrois_nrj(1 + 0.5 * (adiabatic_cst - 1.0)
                               * m_density[ev] * cqs_v_nplus1 / m_cell_mass[ev]) ;
        Real numer_accrois_nrj(m_internal_energy[ev]
                               - (0.5 * (m_pressure_n[ev] + m_pseudo_viscosity_n[ev])
                                  * cqs_v_n / m_cell_mass[ev])
                               - (0.5 * m_pseudo_viscosity[ev]
                                  * cqs_v_nplus1 / m_cell_mass[ev])
                               + (0.25 * (m_pressure_n[ev] + m_pseudo_viscosity_n[ev]) 
                                  * cqs_delta_v / m_cell_mass[ev])
                               - (0.5 * (m_pseudo_viscosity_n[ev] - m_pseudo_viscosity_nmoins1[ev]) 
                                  * cqs_v_old_n / m_cell_mass[ev])
                               );

        m_internal_energy[ev] = numer_accrois_nrj / denom_accrois_nrj;
      }
    }
  }
  PROF_ACC_END;
}   
/**
 *******************************************************************************
 * \file computePressionMoyenne()
 * \brief Calcul de la pression moyenne
 *
 * \param  m_fracvol_env, m_pressure_env_nplus1, m_speed_velocity_env_nplus1
 * \return m_pressure_nplus1, m_speed_velocity_nplus1
 *******************************************************************************
 */
void MahycoModule::
computePressionMoyenne()
{
  if (options()->sansLagrange) return;
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans computePressionMoyenne() ";
  // maille mixte
  // moyenne sur la maille
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      m_pressure[cell] = 0.;
      m_sound_speed[icell] = 1.e-20;
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;        
        m_pressure[cell] += m_fracvol[ev] * m_pressure[ev];
        m_sound_speed[cell] = std::max(m_sound_speed[ev], m_sound_speed[cell]);
      }
    }
  }
  PROF_ACC_END;
}     
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
computeDeltaT()
{
  PROF_ACC_BEGIN(__FUNCTION__);
  debug() << " Entree dans compute DT avec " << m_global_old_deltat()
         << " et " << options()->deltatInit()
         << " et " << m_global_deltat()
         << " et " << options()->deltatMax();
         
  m_global_old_deltat = m_global_deltat;
  m_old_deltat = m_global_old_deltat();
         
  Real new_dt = FloatInfo < Real >::maxValue();
  if (options()->sansLagrange) {
    // on garde le meme pas de temps
    new_dt = options()->deltatInit();
    
  } else {
    CellToAllEnvCellConverter all_env_cell_converter(mm);

    // Calcul du pas de temps pour le respect du critère de CFL
    Real minimum_aux = FloatInfo < Real >::maxValue();
    Integer cell_id(-1), nbenvcell(-1);
    Real cc(0.), ll(0.);

    ENUMERATE_CELL(icell, allCells()){
        Cell cell = * icell;
        Real cell_dx = m_caracteristic_length[icell];
        Real sound_speed = m_sound_speed[icell];
        Real vmax(0.);
        if (options()->withProjection)
          for (NodeEnumerator inode(cell.nodes()); inode.index() < cell.nbNode(); ++inode) {
            vmax = math::max(m_velocity[inode].abs(), vmax);
          }
        Real dx_sound = cell_dx / (sound_speed + vmax);
        minimum_aux = math::min(minimum_aux, dx_sound);
        if (minimum_aux == dx_sound) {
            cell_id = icell.localId();
            cc = m_sound_speed[icell];
            ll = m_caracteristic_length[icell];
            AllEnvCell all_env_cell = all_env_cell_converter[cell];
            nbenvcell = all_env_cell.nbEnvironment();
        }
    }

    new_dt = options()->cfl() * minimum_aux;
    new_dt = parallelMng()->reduce(Parallel::ReduceMin, new_dt);
    
    // respect de taux de croissance max
    new_dt = math::min(new_dt, 1.05 * m_global_old_deltat());
    // respect de la valeur max imposée par le fichier de données .plt
    debug() << " nouveau pas de temps " << new_dt << " par " << cell_id << 
        " (avec " << nbenvcell << " envs) avec " << cc << " " << ll << " et min " << minimum_aux;
    new_dt = math::min(new_dt, options()->deltatMax());
    // respect du pas de temps minimum
    if (new_dt < options()->deltatMin()) {
        info() << " pas de temps minimum ";
        info() << " nouveau pas de temps " << new_dt << " par " << cell_id 
            << " (avec " << nbenvcell << " envs) avec vitson = " << cc << " et longeur  =  " << ll << " et min " << minimum_aux;
        exit(1);
    }
    
    debug() << " nouveau pas de temps2 " << new_dt;
        
  }
  
  // Mise à jour du pas de temps
  m_global_deltat = new_dt;
  
  // Le dernier calcul se fait exactement au temps stopTime()
  Real stop_time  = options()->finalTime();
  bool not_yet_finish = (m_global_time() < stop_time);
  bool finish = (m_global_time() > stop_time);
  bool too_much = ((m_global_time()+new_dt) > stop_time);

  
  debug() << " nouveau pas de temps " << new_dt;
  if ( not_yet_finish && too_much){
    new_dt = stop_time - m_global_time();
    subDomain()->timeLoopMng()->stopComputeLoop(true);
  }
  if (finish) {
    subDomain()->timeLoopMng()->stopComputeLoop(true);
  }
    
  debug() << " time " << m_global_time() << " et fin à " << stop_time;
  debug() << " not_yet_finish " << not_yet_finish;
  debug() << " finish " << finish;
  debug() << " too_much " << too_much;
  
  PROF_ACC_END;
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline void MahycoModule::
computeCQs(Real3 node_coord[8], Real3 face_coord[6], Span<Real3> out_cqs)
{
  const Real3 c0 = face_coord[0];
  const Real3 c1 = face_coord[1];
  const Real3 c2 = face_coord[2];
  const Real3 c3 = face_coord[3];
  const Real3 c4 = face_coord[4];
  const Real3 c5 = face_coord[5];

  // Calcul des normales face 1 :
  const Real3 n1a04 = 0.5 * math::vecMul(node_coord[0] - c0, node_coord[3] - c0);
  const Real3 n1a03 = 0.5 * math::vecMul(node_coord[3] - c0, node_coord[2] - c0);
  const Real3 n1a02 = 0.5 * math::vecMul(node_coord[2] - c0, node_coord[1] - c0);
  const Real3 n1a01 = 0.5 * math::vecMul(node_coord[1] - c0, node_coord[0] - c0);

  // Calcul des normales face 2 :
  const Real3 n2a05 = 0.5 * math::vecMul(node_coord[0] - c1, node_coord[4] - c1);
  const Real3 n2a12 = 0.5 * math::vecMul(node_coord[4] - c1, node_coord[7] - c1);
  const Real3 n2a08 = 0.5 * math::vecMul(node_coord[7] - c1, node_coord[3] - c1);
  const Real3 n2a04 = 0.5 * math::vecMul(node_coord[3] - c1, node_coord[0] - c1);

  // Calcul des normales face 3 :
  const Real3 n3a01 = 0.5 * math::vecMul(node_coord[0] - c2, node_coord[1] - c2);
  const Real3 n3a06 = 0.5 * math::vecMul(node_coord[1] - c2, node_coord[5] - c2);
  const Real3 n3a09 = 0.5 * math::vecMul(node_coord[5] - c2, node_coord[4] - c2);
  const Real3 n3a05 = 0.5 * math::vecMul(node_coord[4] - c2, node_coord[0] - c2);

  // Calcul des normales face 4 :
  const Real3 n4a09 = 0.5 * math::vecMul(node_coord[4] - c3, node_coord[5] - c3);
  const Real3 n4a10 = 0.5 * math::vecMul(node_coord[5] - c3, node_coord[6] - c3);
  const Real3 n4a11 = 0.5 * math::vecMul(node_coord[6] - c3, node_coord[7] - c3);
  const Real3 n4a12 = 0.5 * math::vecMul(node_coord[7] - c3, node_coord[4] - c3);

  // Calcul des normales face 5 :
  const Real3 n5a02 = 0.5 * math::vecMul(node_coord[1] - c4, node_coord[2] - c4);
  const Real3 n5a07 = 0.5 * math::vecMul(node_coord[2] - c4, node_coord[6] - c4);
  const Real3 n5a10 = 0.5 * math::vecMul(node_coord[6] - c4, node_coord[5] - c4);
  const Real3 n5a06 = 0.5 * math::vecMul(node_coord[5] - c4, node_coord[1] - c4);

  // Calcul des normales face 6 :
  const Real3 n6a03 = 0.5 * math::vecMul(node_coord[2] - c5, node_coord[3] - c5);
  const Real3 n6a08 = 0.5 * math::vecMul(node_coord[3] - c5, node_coord[7] - c5);
  const Real3 n6a11 = 0.5 * math::vecMul(node_coord[7] - c5, node_coord[6] - c5);
  const Real3 n6a07 = 0.5 * math::vecMul(node_coord[6] - c5, node_coord[2] - c5);

  // Calcul des résultantes aux sommets :
  out_cqs[0] = (5. * (n1a01 + n1a04 + n2a04 + n2a05 + n3a05 + n3a01) +
                          (n1a02 + n1a03 + n2a08 + n2a12 + n3a06 + n3a09)) * (1. / 12.);
  out_cqs[1] = (5. * (n1a01 + n1a02 + n3a01 + n3a06 + n5a06 + n5a02) +
                          (n1a04 + n1a03 + n3a09 + n3a05 + n5a10 + n5a07)) * (1. / 12.);
  out_cqs[2] = (5. * (n1a02 + n1a03 + n5a07 + n5a02 + n6a07 + n6a03) +
                          (n1a01 + n1a04 + n5a06 + n5a10 + n6a11 + n6a08)) * (1. / 12.);
  out_cqs[3] = (5. * (n1a03 + n1a04 + n2a08 + n2a04 + n6a08 + n6a03) +
                          (n1a01 + n1a02 + n2a05 + n2a12 + n6a07 + n6a11)) * (1. / 12.);
  out_cqs[4] = (5. * (n2a05 + n2a12 + n3a05 + n3a09 + n4a09 + n4a12) +
                          (n2a08 + n2a04 + n3a01 + n3a06 + n4a10 + n4a11)) * (1. / 12.);
  out_cqs[5] = (5. * (n3a06 + n3a09 + n4a09 + n4a10 + n5a10 + n5a06) +
                          (n3a01 + n3a05 + n4a12 + n4a11 + n5a07 + n5a02)) * (1. / 12.);
  out_cqs[6] = (5. * (n4a11 + n4a10 + n5a10 + n5a07 + n6a07 + n6a11) +
                          (n4a12 + n4a09 + n5a06 + n5a02 + n6a03 + n6a08)) * (1. / 12.);
  out_cqs[7] = (5. * (n2a08 + n2a12 + n4a12 + n4a11 + n6a11 + n6a08) +
                          (n2a04 + n2a05 + n4a09 + n4a10 + n6a07 + n6a03)) * (1. / 12.);
}

/*---------------------------------------------------------------------------*/
Real MahycoModule::produit(Real A, Real B, Real C, Real D)
  {
  return (A*B-C*D);
  }
/*---------------------------------------------------------------------------*/
inline Real MahycoModule::norme(Real E, Real F, Real G)
  {
  return (sqrt((E*E)+(F*F)+(G*G)));
  }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_MAHYCO(MahycoModule);
