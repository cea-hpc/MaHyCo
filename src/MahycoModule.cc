// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "MahycoModule.h"

#include <arcane/MathUtils.h>
#include <arcane/IParallelMng.h>
#include <arcane/ITimeLoopMng.h>

#include <arcane/geometry/IGeometry.h>
#include <arcane/mesh/GhostLayerMng.h>

using namespace Arcane;
using namespace Arcane::Materials;


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
hydroStartInit()
{
    IParallelMng* m_parallel_mng = subDomain()->parallelMng();
    my_rank = m_parallel_mng->commRank();


    info() <<  "Mon rang " << my_rank << " et mon nombre de mailles " << allCells().size();
    info() <<  " Mes mailles pures : " << ownCells().size();
    info() <<  " Mes mailles frantomes : " << allCells().size() - ownCells().size();

    info() << " Check donnees ";
    if ( ( options()->remap()->getOrdreProjection() == 3 ) && ( mesh()->ghostLayerMng()->nbGhostLayer() != 3 ) && ( m_parallel_mng->isParallel() == true ) ) {
        info() << " mode parallele : " << m_parallel_mng->isParallel();
        info() << " nombre de couches de mailles fantomes : " << mesh()->ghostLayerMng()->nbGhostLayer();
        info() << " incompatible avec la projection d'ordre " << options()->remap()->getOrdreProjection();
        info() << " ----------------------------- fin du calcul à la fin de l'init ---------------------------------------------";
        subDomain()->timeLoopMng()->stopComputeLoop ( true );
    }
    if ( ( options()->withProjection == true ) && ( mesh()->ghostLayerMng()->nbGhostLayer() < 2 ) && ( m_parallel_mng->isParallel() == true ) ) {
        info() << " mode parallele : " << m_parallel_mng->isParallel();
        info() << " nombre de couches de mailles fantomes : " << mesh()->ghostLayerMng()->nbGhostLayer();
        info() << " incompatible avec la projection ";
        info() << " ----------------------------- fin du calcul à la fin de l'init ---------------------------------------------";
        subDomain()->timeLoopMng()->stopComputeLoop ( true );
    }
    if ( ( options()->remap()->hasProjectionSimplePente() == true ) && options()->remap()->hasProjectionPenteBorne() == false )  {
        info() << " Mode simple Pente activé que dans le cas ou Pente Borne est active";
        info() << " ----------------------------- fin du calcul à la fin de l'init ---------------------------------------------";
        subDomain()->timeLoopMng()->stopComputeLoop ( true );
    }
    m_dimension = mesh()->dimension();
    
    if (options()->withProjection == true ) {
        m_cartesian_mesh = ICartesianMesh::getReference ( mesh() );
        m_cartesian_mesh->computeDirections();
    } 
    
    // Dimensionne les variables tableaux
    m_cell_cqs.resize ( 4* ( m_dimension-1 ) );
    m_cell_cqs_n.resize ( 4* ( m_dimension-1 ) );

    // Initialise le delta-t
    Real deltat_init = options()->deltatInit();
    m_global_deltat = deltat_init;

    info() << " Initialisation des environnements";
    hydroStartInitEnvAndMat();
    
    info() << "Initialise les données géométriques: volume, cqs, longueurs caractéristiques";
    computeGeometricValues();
    
    info() << " Initialisation des groupes de faces";
    PrepareFaceGroup();

    for ( Integer i = 0, nb = options()->boundaryCondition.size(); i < nb; ++i ) {
        if ( options()->boundaryCondition[i]->type() == TypesMahyco::OnFilePressure ) {
            // String fichier =  options()->boundaryCondition[i]->fichier();
            const std::string& fichier = "Pression.data";
            if ( fichier != "RIEN" ) {
                lireFichierCDL ( fichier );
            }
        }
    }
    
    // on considère par défaut que tous les mailles sont pures à l'init
    // cela peut etre surcharger suivant les "casModel" dans initVar
    m_fracvol.fill ( 1.0 );
    m_mass_fraction.fill ( 1.0 );
    // fraction des phases
    m_frac_phase1.fill ( 0.0 );
    m_frac_phase2.fill ( 0.0 );
    m_frac_phase3.fill ( 0.0 );
    m_frac_phase4.fill ( 0.0 );
    m_frac_phase5.fill ( 0.0 );
    m_frac_phase6.fill ( 0.0 );
    // initialisation de la pseudo à 0.
    m_pseudo_viscosity.fill ( 0.0 );
    
    m_maille_endo.fill ( 0.0 );
    m_density_fracture.fill ( 0.0 );
    m_internal_energy_fracture.fill ( 0.0 );
    
    m_internal_energy_n.fill ( 0.0 );
    m_temperature_n.fill ( 0.0 );
    m_temperature.fill ( 0.0 );
    
    if (options()->casModel()->isInternalModel() == true) { 
        info() << " Initialisation des variables pour les cas test internes";
        Real* densite_initiale = (double *)malloc(sizeof(double) * m_nb_env);
        Real* pression_initiale = (double *)malloc(sizeof(double) * m_nb_env);
        Real* energie_initiale = (double *)malloc(sizeof(double) * m_nb_env);
        Real* temperature_initiale=  (double *)malloc(sizeof(double) * m_nb_env);
        Real3x3 vitesse_initiale;
        
        for ( Integer i=0,n=options()->environment().size(); i<n; ++i ) {
            densite_initiale[i] = options()->environment[i].densiteInitiale;
            pression_initiale[i] = options()->environment[i].pressionInitiale;
            energie_initiale[i] = options()->environment[i].energieInitiale;
            vitesse_initiale[i] = options()->environment[i].vitesseInitiale;
            temperature_initiale[i] = options()->environment[i].temperatureInitiale;
            info() << " L'environnement " << options()->environment() [i]->name() << " est initialisé à ";
            info() << " Densité = " << densite_initiale[i];
            info() << " Pression = " << pression_initiale[i];
        }
        
        // Initialises les variables (surcharge l'init d'arcane)
        options()->casModel()->initVar ( m_dimension, densite_initiale, energie_initiale, 
                                        pression_initiale, temperature_initiale, vitesse_initiale );
    } else {
        for ( Integer i=0,n=options()->environment().size(); i<n; ++i ) { 
            if (options()->environment[i].initialisationUtilisateur) {
                // Initialises les variables par programme (surcharge l'init d'arcane)
                options()->casModel()->initUtilisateur();
            }
        }
    }

    // initialisation du time-history
    if ( options()->timeHistory.size() >0 ) {
        initTH();
    }
      
    info() << " - Appel aux eos ";
    if ( !options()->sansLagrange ) {
        for ( Integer i=0,n=options()->environment().size(); i<n; ++i ) {
            IMeshEnvironment* ienv = mm->environments() [i];
            info() << " - Appel aux eos pour l'env " << ienv->name();
            // Initialise l'énergie et la vitesse du son
            options()->environment[i].eosModel()->initEOS ( ienv );
            info() << " - Appel modèle elasto pour l'env " << ienv->name();
            // Initialisation des varibles ElastoPlastic (mise à zéro)
            options()->environment[i].elastoModel()->initElasto ( ienv );
        }

        CellToAllEnvCellConverter all_env_cell_converter ( mm );
        info() << " Calcul Valeurs moyennes ";
        // valeur moyenne
        ENUMERATE_CELL ( icell, allCells() ) {
            Cell cell = * icell;
            AllEnvCell all_env_cell = all_env_cell_converter[cell];
            if ( all_env_cell.nbEnvironment() !=1 ) {
                m_internal_energy[cell] = 0.;
                ENUMERATE_CELL_ENVCELL ( ienvcell,all_env_cell ) {
                    EnvCell ev = *ienvcell;
                    m_internal_energy[cell] += m_internal_energy[ev] * m_mass_fraction[ev];
                    m_sound_speed[cell] = std::max ( m_sound_speed[ev], m_sound_speed[cell] );
                    m_density[cell] += m_density[ev] * m_fracvol[ev];
                    m_pressure[cell] += m_pressure[ev] * m_fracvol[ev];
                }
            }
        }
    }
    info() << " Apres passage dans EOS ";
    // initialisation du volume Euler
    ENUMERATE_CELL ( icell,allCells() ) {
        Cell cell = *icell;
        // initialisation du volume Euler
        m_euler_volume[cell] = m_cell_volume[cell];
    }

    Arcane::Numerics::IGeometryMng* gm = options()->geometry();
    gm->init();
    auto g = gm->geometry();
    bool is_verbose = false;
    if ( is_verbose ) {
        ENUMERATE_CELL ( icell,allCells() ) {
            Cell c = *icell;
            // initialisation du volume Euler
            info() << "Volume = " << g->computeMeasure ( c );
            for ( NodeEnumerator inode ( c.nodes() ); inode.hasNext(); ++inode ) {
                info() << c.localId() << " avec " << inode.localId();
            }
        }
    }

    InitGeometricValues();

    // deplacé de hydrocontinueInit
    // calcul des volumes via les cqs, differents deu calcul dans computeGeometricValues ?
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = *icell;
        // Calcule le volume de la maille et des noeuds
        Real volume = 0.0;
        for ( Integer inode = 0; inode < cell.nbNode(); ++inode ) {
            volume += math::dot ( m_node_coord[cell.node ( inode )], m_cell_cqs[icell] [inode] );
            m_node_volume[cell.node ( inode )] += volume;
        }
        volume /= m_dimension;
        m_cell_volume[icell] = volume;
    }
}

/**
 *******************************************************************************
 * \file computeCellMass()
 * \brief Calcul de la masse des mailles
 *
 * \param  m_cell_volume, m_density, m_mass_fraction_ev
 * \return m_cell_mass, m_cell_mass_env
 *******************************************************************************
 */
void MahycoModule::
computeCellMass()
{
    debug() << " Entree dans computeCellMass()";
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        m_cell_mass[cell] = m_density[cell] * m_cell_volume[cell];
    }
    ENUMERATE_ENV ( ienv,mm ) {
        IMeshEnvironment* env = *ienv;
        ENUMERATE_ENVCELL ( ienvcell,env ) {
            EnvCell ev = *ienvcell;
            Cell cell = ev.globalCell();
            m_cell_mass[ev] = m_mass_fraction[ev] * m_cell_mass[cell];
        }
    }
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
    debug() << " Entree dans computeNodeMass()";
    // Initialisation ou reinitialisation de la masse nodale
    m_node_mass.fill ( 0. );
    Real one_over_nbnode = m_dimension == 2 ? .25  : .125 ;
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        Real contrib_node_mass = one_over_nbnode * m_cell_mass[cell];
        for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
            m_node_mass[inode] += contrib_node_mass;
        }
    }
    m_node_mass.synchronize();
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
    if ( subDomain()->isContinue() ) {

        debug() << " Entree dans hydroContinueInit()";
        // en reprise
        m_cartesian_mesh = ICartesianMesh::getReference ( mesh() );
        m_dimension = mesh()->dimension();
        m_cartesian_mesh->computeDirections();

        mm = IMeshMaterialMng::getReference ( defaultMesh() );

        mm->recreateFromDump();
        m_nb_env = mm->environments().size();
        m_nb_vars_to_project = 3 * m_nb_env + 3 + 1 + 1;

        m_global_old_deltat = m_old_deltat;
        // mise a jour nombre iteration
        m_global_iteration = m_global_iteration() +1;

        if ( !options()->sansLagrange ) {

            debug() << " Nombre Mailles en reprise " << allCells().size();

            for ( Integer i=0,n=options()->environment().size(); i<n; ++i ) {
                IMeshEnvironment* ienv = mm->environments() [i];
                // Re-Initialise certains parametres pour l'EOS
                options()->environment[i].eosModel()->ReinitEOS ( ienv );
            }
        }
        
        IParallelMng* m_parallel_mng = subDomain()->parallelMng();
        my_rank = m_parallel_mng->commRank();
    }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MahycoModule::
saveValuesAtN()
{
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
    m_strain_tensor.synchronize();
    m_cell_cqs.synchronize();
    m_velocity.synchronize();
    // avec Elasto
    m_strain_tensor.synchronize();
    m_spin_rate.synchronize();
    m_deformation_rate.synchronize();
    m_velocity_gradient.synchronize();
    m_plastic_deformation_velocity.synchronize();
    m_plastic_deformation.synchronize();
    // pour les phases
    m_frac_phase1.synchronize();
    m_frac_phase2.synchronize();
    m_frac_phase3.synchronize();
    m_frac_phase4.synchronize();
    m_frac_phase5.synchronize();
    m_frac_phase6.synchronize();


    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = *icell;
        m_pseudo_viscosity_nmoins1[cell] = m_pseudo_viscosity_n[cell];
        m_pseudo_viscosity_n[cell] = m_pseudo_viscosity[cell];
        m_pressure_n[cell] = m_pressure[cell];
        m_strain_tensor_n[cell] =m_strain_tensor[cell];
        m_cell_volume_n[cell] = m_cell_volume[cell];
        m_density_n[cell] = m_density[cell];
        m_internal_energy_n[cell] = m_internal_energy[cell];
        m_temperature_n[cell] = m_temperature[cell];
    }
    ENUMERATE_ENV ( ienv,mm ) {
        IMeshEnvironment* env = *ienv;
        ENUMERATE_ENVCELL ( ienvcell,env ) {
            EnvCell ev = *ienvcell;
            m_pseudo_viscosity_nmoins1[ev] = m_pseudo_viscosity_n[ev];
            m_pseudo_viscosity_n[ev] = m_pseudo_viscosity[ev];
            m_pressure_n[ev] = m_pressure[ev];
            m_strain_tensor_n[ev] =m_strain_tensor[ev];
            m_cell_volume_n[ev] = m_cell_volume[ev];
            m_density_n[ev] = m_density[ev];
            m_internal_energy_n[ev] = m_internal_energy[ev];
            m_temperature_n[ev] = m_temperature[ev];
        }
    }

    m_cell_cqs_n.copy ( m_cell_cqs );

    if ( !options()->sansLagrange ) {
        m_velocity_n.copy ( m_velocity );
    }
//    FAIT en fin de projection
//    if ( options()->withProjection && options()->remap()->isEuler() != "NON" ) {
//        // info() << " recopie des valeurs initiales ";
//        m_node_coord.copy ( m_node_coord_0 );
//    }

}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
computeArtificialViscosity()
{
    if ( options()->sansLagrange ) {
        return;
    }
    debug() << " Entree dans computeArtificialViscosity()";
    ENUMERATE_ENV ( ienv,mm ) {
        IMeshEnvironment* env = *ienv;
        Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst ( env );
        Real a1 = options()->environment[env->id()].linearPseudoCoeff();
        Real a2 = options()->environment[env->id()].quadraticPseudoCoeff();
        ENUMERATE_ENVCELL ( ienvcell,env ) {
            EnvCell ev = *ienvcell;
            Cell cell = ev.globalCell();
            m_pseudo_viscosity[ev] = 0.;
            if ( m_div_u[cell] < 0.0 ) {
                m_pseudo_viscosity[ev] = 1. / m_tau_density[ev]
                                         * ( - a1 * m_caracteristic_length[cell] * m_sound_speed[cell] * m_div_u[cell]
                                             + a2 * m_caracteristic_length[cell] * m_caracteristic_length[cell]
                                             * m_div_u[cell] * m_div_u[cell] );
            }
        }
    }


    // maille mixte
    // moyenne sur la maille
    CellToAllEnvCellConverter all_env_cell_converter ( mm );
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        if ( all_env_cell.nbEnvironment() !=1 ) {
            m_pseudo_viscosity[cell] = 0.;
            ENUMERATE_CELL_ENVCELL ( ienvcell,all_env_cell ) {
                EnvCell ev = *ienvcell;
                m_pseudo_viscosity[cell] += m_pseudo_viscosity[ev] * m_fracvol[ev];
            }
        }
    }
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
    if ( options()->sansLagrange ) {
        updateVelocityWithoutLagrange();
        return;
    }

    // passage des vitesse de n à n-1/2 si projection
    // la vitesse m_velocity_n
    // (qui dans le cas de vnr-csts est la vitesse apres projection du pas de temps précédent donc n),
    // est replacee à n-1/2 pour vnr-csts.
    // Dans le cas de vnr (pas csts) elle est deja en n-1/2
    if ( options()->schemaCsts() && options()->withProjection ) {
        updateVelocityBackward();
    }
    
    // pinfo() << my_rank << " : " << " Entree dans updateVelocity()";
    // Remise à zéro du vecteur des forces.
    m_force.fill ( Real3::zero() );

    // Calcul pour chaque noeud de chaque maille la contribution
    // des forces de pression
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        Real pressure = m_pressure_n[icell] + m_pseudo_viscosity_n[icell];
        for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
            m_force[inode] += pressure * m_cell_cqs_n[icell] [inode.index()];
        }
    }
    // Calcul de la force Elasto-plastique
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        Real3x3 force_tensor = -  m_strain_tensor_n[icell];
        for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
            m_force[inode] += math::prodTensVec ( force_tensor, m_cell_cqs_n[icell] [inode.index()] );
        }
    }
    
    const Real dt ( 0.5 * ( m_global_old_deltat() + m_global_deltat() ) );
    
//     VariableNodeReal3InView in_force ( viewIn ( m_force ) );
//     VariableNodeReal3InView in_velocity ( viewIn ( m_velocity_n ) );
//     VariableNodeRealInView  in_mass ( viewIn ( m_node_mass ) );
//     VariableNodeReal3OutView out_velocity ( viewOut ( m_velocity ) );
// 
//     Calcul l'impulsion aux noeuds
//     PRAGMA_IVDEP
//     ENUMERATE_SIMD_NODE ( inode, allNodes() ) {
//         SimdNode snode=*inode;
//         out_velocity[snode] = in_velocity[snode] + ( dt / in_mass[snode] ) * in_force[snode];;
//     }
    
    // Calcul l'impulsion aux noeuds
    ENUMERATE_NODE ( inode, allNodes() ) {
        Node node = *inode;
        m_velocity[inode] += (dt / m_node_mass[inode]) * m_force[inode];
        /* if (node.uniqueId() ==  177437 || node.uniqueId() ==  177438 || node.uniqueId() == 178439  || node.uniqueId() == 178438 ) {
            pinfo() << my_rank << " " << node.uniqueId() << " Vitesse " << m_velocity[node];
            pinfo() << my_rank << " " << node.uniqueId() << " Masse " << m_node_mass[inode];
            pinfo() << my_rank << " " << node.uniqueId() << " Force " << m_force[inode];
        } */
    }
     
    Real3 g = options()->gravity;
    ENUMERATE_NODE ( inode, allNodes() ) {
        Node node = *inode;
        m_velocity[inode] +=  dt * g;
    }
    m_velocity.synchronize();
    
    // pinfo() << my_rank << " : " << " Fin dans updateVelocity()";
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
    if ( options()->sansLagrange ) {
        return;
    }
    debug() << " Entree dans updateVelocityBackward()";
    // Remise à zéro du vecteur des forces.
    m_force.fill ( Real3::null() );

    // Calcul pour chaque noeud de chaque maille la contribution
    // des forces de pression
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        Real pressure = m_pressure_n[icell] + m_pseudo_viscosity_n[icell];
        for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
            m_force[inode] += pressure * m_cell_cqs_n[icell] [inode.index()];
        }
    }

    // Calcul de la force Elasto-plastique
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        Real3x3 force_tensor = -  m_strain_tensor_n[icell];
        for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
            m_force[inode] += math::prodTensVec ( force_tensor, m_cell_cqs_n[icell] [inode.index()] );
        }
    }
    const Real dt ( -0.5 * m_global_old_deltat() );

//     VariableNodeReal3InView in_force ( viewIn ( m_force ) );
//     VariableNodeReal3InView in_velocity ( viewIn ( m_velocity_n ) );
//     VariableNodeRealInView  in_mass ( viewIn ( m_node_mass ) );
//     VariableNodeReal3OutView out_velocity ( viewOut ( m_velocity_n ) );
// 
//     Calcule l'impulsion aux noeuds
//     PRAGMA_IVDEP
//     ENUMERATE_SIMD_NODE ( inode, allNodes() ) {
//         SimdNode snode=*inode;
//         out_velocity[snode] = in_velocity[snode] + ( dt / in_mass[snode] ) * in_force[snode];;
//     }
    
    // Calcul l'impulsion aux noeuds
    ENUMERATE_NODE ( inode, allNodes() ) {
        Node node = *inode;
        m_velocity[inode] += (dt / m_node_mass[inode]) * m_force[inode];
    }
     
    Real3 g = options()->gravity;
    ENUMERATE_NODE ( inode, allNodes() ) {
        Node node = *inode;
        m_velocity[inode] +=  dt * g;
    }
    m_velocity_n.synchronize();
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
    if ( options()->sansLagrange ) {
        return;
    }
    debug() << " Entree dans updateVelocityForward()";
    // Remise à zéro du vecteur des forces.
    m_force.fill ( Real3::null() );

    // Calcul pour chaque noeud de chaque maille la contribution
    // des forces de pression
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        Real pressure = m_pressure[icell] + m_pseudo_viscosity[icell];
        for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
            m_force[inode] += pressure * m_cell_cqs[icell] [inode.index()];
        }
    }

    // Calcul de la force Elasto-plastique
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        Real3x3 force_tensor = -  m_strain_tensor[icell];
        for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
            m_force[inode] += math::prodTensVec ( force_tensor, m_cell_cqs[icell] [inode.index()] );
        }
    }
    const Real dt ( 0.5 * m_global_deltat() );

//     VariableNodeReal3InView in_force ( viewIn ( m_force ) );
//     VariableNodeReal3InView in_velocity ( viewIn ( m_velocity ) );
//     VariableNodeRealInView  in_mass ( viewIn ( m_node_mass ) );
//     VariableNodeReal3OutView out_velocity ( viewOut ( m_velocity ) );
// 
//     Calcule l'impulsion aux noeuds
//     PRAGMA_IVDEP
//     ENUMERATE_SIMD_NODE ( inode, allNodes() ) {
//         SimdNode snode=*inode;
//         out_velocity[snode] = in_velocity[snode] + ( dt / in_mass[snode] ) * in_force[snode];;
//     }
    
    // Calcul l'impulsion aux noeuds
    ENUMERATE_NODE ( inode, allNodes() ) {
        Node node = *inode;
        m_velocity[inode] += (dt / m_node_mass[inode]) * m_force[inode];
    }
    
    // Ajout du terme de gravité : TODO à transformer
    Real3 g = options()->gravity;
    ENUMERATE_NODE ( inode, allNodes() ) {
        Node node = *inode;
        m_velocity[inode] +=  dt * g;
    }
    m_velocity.synchronize();
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
    bool option ( options()->casModel()->hasReverseOption() );
    Real factor ( options()->casModel()->getReverseParameter() );

    Real option_real ( ( Real ) option );

    VariableNodeReal3InView in_velocity ( viewIn ( m_velocity_n ) );
    VariableNodeReal3OutView out_velocity ( viewOut ( m_velocity ) );

    PRAGMA_IVDEP
    ENUMERATE_SIMD_NODE ( inode, allNodes() ) {
        SimdNode snode=*inode;
        out_velocity[snode] = in_velocity[snode] * ( 1. -option_real )
                              + option_real * in_velocity[snode] * cos ( Pi * m_global_time() * factor );
    }
    m_velocity.synchronize();
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
    debug() << my_rank << " : " << " Entree dans updatePosition()";
    Real deltat = m_global_deltat();
    ENUMERATE_NODE ( inode, allNodes() ) {
        Node node = *inode;
        if ( ( ( options()->sansLagrange ) && ( node.nbCell() == 4 ) ) || ( !options()->sansLagrange ) ) {
            m_node_coord[inode] += deltat * m_velocity[inode];
        }
    }
    Real one_over_nbnode = m_dimension == 2 ? .25  : .125 ;
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        Real3 somme = {0., 0., 0.};
        for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
            somme += m_node_coord[inode];
        }
        m_cell_coord[cell] = one_over_nbnode * somme;
    }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
applyBoundaryCondition()
{
    debug() << my_rank << " : " << " Entree dans applyBoundaryCondition()";
    for ( Integer i = 0, nb = options()->boundaryCondition.size(); i < nb; ++i ) {
        String NomBC = options()->boundaryCondition[i]->surface;
        FaceGroup face_group = mesh()->faceFamily()->findGroup ( NomBC );
        Real value = options()->boundaryCondition[i]->value();
        TypesMahyco::eBoundaryCondition type = options()->boundaryCondition[i]->type();

        bool Novalue = false;

        if ( m_global_time() < options()->boundaryCondition[i]->Tdebut ) {
            Novalue=true;
        }
        if ( m_global_time() > options()->boundaryCondition[i]->Tfin ) {
            Novalue=true;
        }

        if ( Novalue == false ) {
            // boucle sur les faces de la surface
            ENUMERATE_FACE ( j, face_group ) {
                Face face = * j;
                Integer nb_node = face.nbNode();

                // boucle sur les noeuds de la face
                for ( Integer k = 0; k < nb_node; ++k ) {
                    Node node = face.node ( k );
                    Real3& velocity = m_velocity[node];
                    Real3& coord = m_node_coord[node];
                    switch ( type ) {
                    case TypesMahyco::VelocityX:
                        velocity.x = value;
                        break;
                    case TypesMahyco::VelocityY:
                        velocity.y = value;
                        break;
                    case TypesMahyco::VelocityZ:
                        velocity.z = value;
                        break;
                    case TypesMahyco::OnFileVelocity:
                        velocity.x = - coord.x / coord.normL2() ;
                        velocity.y = - coord.y / coord.normL2() ;
                    case TypesMahyco::Unknown:
                        break;
                    }
                }
            }
        }
    }
    applyBoundaryConditionForCellVariables();
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
applyBoundaryConditionForCellVariables()
{
    debug() << my_rank << " : " << " Entree dans applyBoundaryConditionForCellVariables()";

    Real one_over_nbnode = m_dimension == 2 ? .25  : .125 ;
    Real one_over_nbnodeface = m_dimension == 2 ? .5  : .25 ;

    const Real dt ( 0.5 * ( m_global_old_deltat() + m_global_deltat() ) );

    for ( Integer i = 0, nb = options()->boundaryCondition.size(); i < nb; ++i ) {
        String NomBC = options()->boundaryCondition[i]->surface;
        FaceGroup face_group = mesh()->faceFamily()->findGroup ( NomBC );
        Real value = options()->boundaryCondition[i]->value();
        TypesMahyco::eBoundaryCondition type = options()->boundaryCondition[i]->type();
        if ( m_global_time() < options()->boundaryCondition[i]->Tdebut ) {
            value = 0.;
        }
        if ( m_global_time() > options()->boundaryCondition[i]->Tfin ) {
            value = 0.;
        }

        // boucle sur les faces de la surface
        ENUMERATE_FACE ( j, face_group ) {
            Face face = * j;
            Integer nb_node = face.nbNode();

            if ( type == TypesMahyco::Density ) {
                Cell cell = face.backCell();
                if ( cell.localId() != -1 ) {
                    m_density[cell] = value;
                    m_cell_mass[cell] = m_density[cell] * m_cell_volume[cell];
                    Real contrib_node_mass = one_over_nbnode * m_cell_mass[cell];
                    for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
                        m_node_mass[inode] += contrib_node_mass;
                    }
                }
                cell = face.frontCell();
                if ( cell.localId() != -1 ) {
                    m_density[cell] = value;
                    m_cell_mass[cell] = m_density[cell] * m_cell_volume[cell];
                    Real contrib_node_mass = one_over_nbnode * m_cell_mass[cell];
                    for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
                        m_node_mass[inode] += contrib_node_mass;
                    }
                }
            }
            if ( type == TypesMahyco::Energy ) {
                value = options()->boundaryCondition[i]->value();
                Cell cell = face.backCell();
                if ( cell.localId() != -1 ) {
                    // pinfo() << " cell coor " <<  m_cell_coord[cell] << "  " << value << " cut " << options()->boundaryCondition[i]->cutoffY;
                    if (math::abs(m_cell_coord[cell].x)  > options()->boundaryCondition[i]->cutoffX) value=0.;
                    if (math::abs(m_cell_coord[cell].y)  > options()->boundaryCondition[i]->cutoffY) value=0.;
                    // pinfo() << " cell coor " <<  m_cell_coord[cell] << "  " << value << " cut " << options()->boundaryCondition[i]->cutoffY;
                    m_internal_energy[cell] += value;
                }
                cell = face.frontCell();
                if ( cell.localId() != -1 ) {
                    // pinfo() << " cell coor " <<  m_cell_coord[cell] << "  " << value << " cut " << options()->boundaryCondition[i]->cutoffY;
                    if (math::abs(m_cell_coord[cell].x)  > options()->boundaryCondition[i]->cutoffX) value=0.;
                    if (math::abs(m_cell_coord[cell].y)  > options()->boundaryCondition[i]->cutoffY) value=0.;
                    // pinfo() << " cell coor " <<  m_cell_coord[cell] << "  " << value << " cut " << options()->boundaryCondition[i]->cutoffY;
                    m_internal_energy[cell] += value;
                }
            }
            if ( type == TypesMahyco::SuperGaussianEnergy ) {
                value = options()->boundaryCondition[i]->value();
                Cell cell = face.backCell();
                if ( cell.localId() != -1 ) {
                    // Fonction super Gaussienne limitée en taille (seule la dependance en Y est codé pour l'instant)
                    Real ay= m_cell_coord[cell].y * options()->boundaryCondition[i]->dependanceY;
                    value *= math::exp ( -math::pow ( ay,4. ) );
                    if (math::abs(m_cell_coord[cell].x)  > options()->boundaryCondition[i]->cutoffX) value=0.;
                    if (math::abs(m_cell_coord[cell].y)  > options()->boundaryCondition[i]->cutoffY) value=0.;
                    m_internal_energy[cell] += value;
                }
                cell = face.frontCell(); 
                if ( cell.localId() != -1 ) {
                    // Fonction super Gaussienne limitée en taille (seule la dependance en Y est codé pour l'instant)
                    Real ay= m_cell_coord[cell].y * options()->boundaryCondition[i]->dependanceY;
                    if (math::abs(m_cell_coord[cell].x)  > options()->boundaryCondition[i]->cutoffX) value=0.;
                    if (math::abs(m_cell_coord[cell].y)  > options()->boundaryCondition[i]->cutoffY) value=0.;
                    value *= math::exp ( -math::pow ( ay,4. ) );
                    m_internal_energy[cell] += value;
                }
            }
            
            if ( type == TypesMahyco::Pressure ) {
                Cell cell = face.backCell();
                if ( cell.localId() != -1 ) {
                    m_pressure[cell] = value;
                }
                cell = face.frontCell();
                if ( cell.localId() != -1 ) {
                    m_pressure[cell] = value;
                }
            }
            if ( type == TypesMahyco::GeometricPressure ||
                    type == TypesMahyco::LinearPressure ||
                    type == TypesMahyco::OnFilePressure ||
                    type == TypesMahyco::SuperGaussianPressure ||
                    type == TypesMahyco::ContactHerzPressure ) {

                if ( type != TypesMahyco::OnFilePressure ) {
                    // boucle sur les noeuds de la face
                    for ( Integer k = 0; k < nb_node; ++k ) {
                        value = options()->boundaryCondition[i]->value();
                        Node node = face.node ( k );
                        Real3 normale_sortante;
                        Real surface ( 0. );
                        if ( m_dimension == 2 ) {
                            Node node2;
                            if ( k<nb_node-1 ) {
                                node2 = face.node ( k+1 );
                            } else {
                                node2 = face.node ( 0 );
                            }
                            surface = ( m_node_coord[node2] - m_node_coord[node] ).normL2();
                        } else { // if (m_dimension == 3) {
                            Real3 face_vec1 = m_node_coord[face.node ( 2 )] - m_node_coord[face.node ( 0 )];
                            Real3 face_vec2 = m_node_coord[face.node ( 3 )] - m_node_coord[face.node ( 1 )];
                            surface = 0.5 * produit ( face_vec1.x, face_vec2.y, face_vec1.y, face_vec2.x ); // pour CDL Z
                            // surface = produit(face_vec1.y, face_vec2.z, face_vec1.z, face_vec2.y); // pour CDL X
                            // surface = - produit(face_vec2.x, face_vec1.z, face_vec2.z, face_vec1.x); // pour CDL Y
                        }
                        Cell cell, cellok;
                        cell = face.backCell();
                        if ( cell.localId() != -1 ) {
                            normale_sortante = ( m_face_coord[face]-m_cell_coord[cell] )
                                               / ( m_face_coord[face]-m_cell_coord[cell] ).normL2();
                            cellok = cell;
                        }
                        cell = face.frontCell();
                        if ( cell.localId() != -1 ) {
                            normale_sortante = ( m_face_coord[face]-m_cell_coord[cell] )
                                               / ( m_face_coord[face]-m_cell_coord[cell] ).normL2();
                            cellok = cell;
                        }

                        if ( type == TypesMahyco::LinearPressure ) {
                            // Fonction lineaire en X, Y ou Z et T
                            value += options()->boundaryCondition[i]->dependanceX * m_node_coord[node].x;
                            value += options()->boundaryCondition[i]->dependanceY * m_node_coord[node].y;
                            value += options()->boundaryCondition[i]->dependanceZ * m_node_coord[node].z;
                            value += options()->boundaryCondition[i]->dependanceT * m_global_time();
                        } else if ( type == TypesMahyco::SuperGaussianPressure ) {
                            Real power = options()->boundaryCondition[i]->power;
                            Real ay= m_node_coord[node].y * options()->boundaryCondition[i]->dependanceY;
                            value *= math::exp ( -math::pow ( ay,power ) );
                            Real ax= m_node_coord[node].x * options()->boundaryCondition[i]->dependanceX;
                            value *= math::exp ( -math::pow ( ax,power ) );
                            if (math::abs(m_node_coord[node].y)  > options()->boundaryCondition[i]->cutoffY) value=0.;
                        } else if ( type == TypesMahyco::ContactHerzPressure ) {
                            // Fonction Contact de Herz sur ZMAX
                            Real a2 = options()->boundaryCondition[i]->dependanceX;
                            Real b2 = options()->boundaryCondition[i]->dependanceY;
                            Real x2 = m_cell_coord[cellok].x * m_cell_coord[cellok].x;
                            Real y2 = m_cell_coord[cellok].y * m_cell_coord[cellok].y; // pour CDL Z
                            // Real r2 = m_cell_coord[cellok].y * m_cell_coord[cellok].y + m_cell_coord[cellok].z * m_cell_coord[cellok].z;// pour CDL X
                            // Real r2 = m_cell_coord[cellok].x * m_cell_coord[cellok].x + m_cell_coord[cellok].z * m_cell_coord[cellok].z;// pour CDL Y
                            if ( x2/a2 + y2/b2 < 1. ) {
                                value *= math::sqrt ( 1-x2/a2-y2/b2 );
                            } else {
                                value = 0.;
                            }
                        }
                        if ( m_global_time() < options()->boundaryCondition[i]->Tdebut ) {
                            value = 0.;
                        }
                        if ( m_global_time() > options()->boundaryCondition[i]->Tfin ) {
                            value = 0.;
                        }

                        // if (value != 0.) pinfo() << cellok.localId()  << " coord " << m_cell_coord[cellok] << " Pression " <<  r2 << " face " << normale_sortante << " dt " << dt << " et value " << value;
                        // pinfo() << face_group.size() << " face Id " << face.localId() ;
                        m_velocity[node] -= dt*one_over_nbnodeface*normale_sortante*value*surface/m_node_mass[node];
                        // pinfo() << " vitesse imposée" << m_velocity[node] << " m_node_mass " << m_node_mass[node] << " dt " << dt << " surface " << surface;
                    }
                } else {
                    Real alpha ( 0. );
                    value =0.;
                    // info() << " on cherche la CDL de pression " << m_table_temps[0] << " " << m_table_temps[1] << " " << m_table_temps[2];
                    Integer it ( -1 );
                    // valeur de la pression lue dans un fichier et interpoler
                    for ( Integer i = 0; i < m_taille_table; ++i ) {
                        // si je trouve un temps plus grand que le temps courant
                        if ( m_global_time() < m_table_temps[i] ) {
                            it = i;
                            break;
                        }
                        // sinon value reste 0.
                    }
                    if ( it>0 ) {
                        alpha = ( m_global_time() - m_table_temps[it-1] ) / ( m_table_temps[it] - m_table_temps[it-1] );
                        alpha = math::min ( alpha, 1. );
                        alpha = math::max ( alpha, 0. );
                        value = m_table_pression[it-1] + alpha * ( m_table_pression[it] - m_table_pression[it-1] );
                    }

                    // boucle sur les noeuds de la face
                    for ( Integer k = 0; k < nb_node; ++k ) {
                        Node node = face.node ( k );
                        Node node2;
                        if ( k<nb_node-1 ) {
                            node2 = face.node ( k+1 );
                        } else {
                            node2 = face.node ( 0 );
                        }
                        Real3 normale_sortante;
                        Real surface = ( m_node_coord[node2] - m_node_coord[node] ).normL2();
                        Cell cell;
                        cell = face.backCell();
                        if ( cell.localId() != -1 ) normale_sortante = ( m_face_coord[face]-m_cell_coord[cell] )
                                    / ( m_face_coord[face]-m_cell_coord[cell] ).normL2();
                        cell = face.frontCell();
                        if ( cell.localId() != -1 ) normale_sortante = ( m_face_coord[face]-m_cell_coord[cell] )
                                    / ( m_face_coord[face]-m_cell_coord[cell] ).normL2();
                        // if (cell.localId() == 0)  pinfo() << " vitesse avant" << m_velocity[node]<< " face " << normale_sortante << " dt " << dt << " et value " << value;
                        m_velocity[node] -= dt*one_over_nbnodeface*normale_sortante*value*surface/m_node_mass[node];
                        // if (cell.localId() == 0) pinfo() << " vitesse imposée" << m_velocity[node] << " m_node_mass " << m_node_mass[node] << " dt " << dt << " surface " << surface;
                    }
                }
            }
        }
    }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
InitGeometricValues()
{
    debug() << " Entree dans InitGeometricValues() ";
    ENUMERATE_NODE ( inode, allNodes() ) {
        m_node_coord_0[inode] = m_node_coord[inode];
    }
    Real one_over_nbnode = m_dimension == 2 ? .5  : .25 ;
    Real3 dirx = {1,0,0};
    Real3 diry = {0,1,0};
    Real3 dirz = {0,0,1};
    if ( m_dimension == 3 ) {
        ENUMERATE_FACE ( iFace, allFaces() ) {
            Face face = *iFace;
            Real3 face_vec1 = m_node_coord[face.node ( 2 )] - m_node_coord[face.node ( 0 )];
            Real3 face_vec2 = m_node_coord[face.node ( 3 )] - m_node_coord[face.node ( 1 )];
            m_face_normal[iFace].x = produit ( face_vec1.y, face_vec2.z, face_vec1.z, face_vec2.y );
            m_face_normal[iFace].y = - produit ( face_vec2.x, face_vec1.z, face_vec2.z, face_vec1.x );
            m_face_normal[iFace].z = produit ( face_vec1.x, face_vec2.y, face_vec1.y, face_vec2.x );
            m_face_normal[iFace] /= m_face_normal[iFace].normL2();
        }
    } else {
        ENUMERATE_FACE ( iFace, allFaces() ) {
            Face face = *iFace;
            m_face_normal[iFace].x = ( m_node_coord[face.node ( 1 )].y - m_node_coord[face.node ( 0 )].y );
            m_face_normal[iFace].y = ( m_node_coord[face.node ( 1 )].x - m_node_coord[face.node ( 0 )].x );
            m_face_normal[iFace] /= m_face_normal[iFace].normL2();
        }
    }
    ENUMERATE_FACE ( iFace, allFaces() ) {
        Face face = *iFace;
        m_face_coord[face] = 0.;
        for ( Integer inode = 0; inode < face.nbNode(); ++inode ) {
            m_face_coord[face] +=  one_over_nbnode * m_node_coord[face.node ( inode )];
        }
        // Calcul d'orientation des faces (initiale sur cartesien)
        if (std::fabs(math::dot(m_face_normal[face], dirx)) >= 1.0E-10) {                                   
            m_face_orientation[face] = 0;
        } else if (std::fabs(math::dot(m_face_normal[face], diry)) >= 1.0E-10) {                                   
            m_face_orientation[face] = 1;
        } else if (std::fabs(math::dot(m_face_normal[face], dirz)) >= 1.0E-10) {                                   
            m_face_orientation[face] = 2;
        }
    }
    ENUMERATE_CELL ( icell,allCells() ) {
        ENUMERATE_FACE ( iface, ( *icell ).faces() ) {
            const Face& face = *iface;
            Integer index = iface.index();
            m_outer_face_normal[icell][index] = ( m_face_coord[face]-m_cell_coord[icell] )
                                                / ( m_face_coord[face]-m_cell_coord[icell] ).normL2();
        } 
    }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------Vitesse Z fixée------------------------------------------*/

void MahycoModule::
computeGeometricValues()
{
    debug() << my_rank << " : " << " Entree dans computeGeometricValues() ";
    // Copie locale des coordonnées des sommets d'une maille
    Real3 coord[8];
    // Coordonnées des centres des faces
    Real3 face_coord[6];

    CellToAllEnvCellConverter all_env_cell_converter ( mm );
    Real racine = m_dimension == 2 ? .5  : 1./3. ;
    m_node_coord.synchronize();
    if ( m_dimension == 3 ) {
        ENUMERATE_CELL ( icell, allCells() ) {
            Cell cell = * icell;
            // Recopie les coordonnées locales (pour le cache)
            for ( NodeEnumerator inode ( cell.nodes() ); inode.index() < 8; ++inode ) {
                coord[inode.index()] = m_node_coord[inode];
            }
            // Calcul les coordonnées des centres des faces
            face_coord[0] = 0.25 * ( coord[0] + coord[3] + coord[2] + coord[1] );
            face_coord[1] = 0.25 * ( coord[0] + coord[4] + coord[7] + coord[3] );
            face_coord[2] = 0.25 * ( coord[0] + coord[1] + coord[5] + coord[4] );
            face_coord[3] = 0.25 * ( coord[4] + coord[5] + coord[6] + coord[7] );
            face_coord[4] = 0.25 * ( coord[1] + coord[2] + coord[6] + coord[5] );
            face_coord[5] = 0.25 * ( coord[2] + coord[3] + coord[7] + coord[6] );

            // Calcule les résultantes aux sommets
            computeCQs ( coord, face_coord, cell );
        }
    } else {
        Real3 npc[5];
        ENUMERATE_CELL ( icell, allCells() ) {
            Cell cell = * icell;
            // Recopie les coordonnées locales (pour le cache)
            for ( NodeEnumerator inode ( cell.nodes() ); inode.index() < cell.nbNode(); ++inode ) {
                coord[inode.index()] = m_node_coord[inode];
            }
            coord[4] = coord[0];
            for ( NodeEnumerator inode ( cell.nodes() ); inode.index() < cell.nbNode(); ++inode ) {
                npc[inode.index()+1].x = 0.5 * ( coord[inode.index()+1].y -  coord[inode.index()].y );
                npc[inode.index()+1].y = 0.5 * ( coord[inode.index()].x -  coord[inode.index()+1].x );
                // npc[inode.index()+1] = npc[inode.index()+1] / npc[inode.index()+1].normL2();
            }
            npc[0] = npc[4];
            for ( Integer ii = 0; ii < 4; ++ii ) {
                m_cell_cqs[icell][ii] = npc[ii+1] + npc[ii];
            }
        }
    }
    m_cell_cqs.synchronize();

    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        // Calcule le volume de la maille
        {
            Real volume = 0.;

            for ( Integer inode = 0; inode < cell.nbNode(); ++inode ) {
                volume += math::dot ( m_node_coord[cell.node ( inode )], m_cell_cqs[icell] [inode] );
                // pinfo() << cell.localId() << " coor " << m_node_coord[cell.node(inode)] << " et " << m_cell_cqs[icell] [inode];
                // m_node_volume[cell.node ( inode )] += volume;
            }
            volume /= m_dimension;

            m_cell_volume[cell] = volume;
            
            if ( volume < 0. ) {
                pinfo() << " ------------------------------------------------------------";
                pinfo() << " Volume Négatif";
                pinfo() << cell.uniqueId() << " : " << " calcul du volume=" << volume;
                pinfo() << "de coordonnées : "  << m_cell_coord[cell];
                pinfo() << " Energie Interne de la maille " << m_internal_energy[cell];
                pinfo() << " Pression de la maille " << m_pressure[cell];
                pinfo() << " Pseudo " << m_pseudo_viscosity[cell];
                pinfo() << " -------------------------------------";
                ENUMERATE_NODE(inode, cell.nodes()) {
                     Node node = *inode; 
                     pinfo() << node.uniqueId() << " Vitesse " << m_velocity[node];
                     ENUMERATE_CELL(icell, node.cells()) {
                          Cell cell2 = * icell;
                          pinfo() << " Environnement Cell :" << cell2.uniqueId() << " ESTMixte " << m_est_mixte[cell2] << " Energie Interne " << m_internal_energy[cell2]
                            << " Pression moyenne " << m_pressure[cell2] << " Pseudo " << m_pseudo_viscosity[cell2];
                          AllEnvCell all_env_cell = all_env_cell_converter[cell2];
                            if ( all_env_cell.nbEnvironment() !=1 ) {
                                ENUMERATE_CELL_ENVCELL ( ienvcell,all_env_cell ) {
                                    EnvCell ev = *ienvcell; 
                                    int index_env = ev.environmentId(); 
                                    pinfo() << " ENV INDEX  " << index_env << " Energie Interne ENV " << m_internal_energy[ev] << " Fraction de volume ENV " << m_fracvol[ev] << " Pression ENV " << m_pressure[ev];
                                }
                            }
                     }
                }
                exit(1);
            }
        }
        // Calcule la longueur caractéristique de la maille.
        {
            if ( options()->longueurCaracteristique() == "faces-opposees" ) {
                // Recopie les coordonnées locales (pour le cache)
                for ( NodeEnumerator inode ( cell.nodes() ); inode.index() < 8; ++inode ) {
                    coord[inode.index()] = m_node_coord[inode];
                }
                // Calcul les coordonnées des centres des faces
                face_coord[0] = 0.25 * ( coord[0] + coord[3] + coord[2] + coord[1] );
                face_coord[1] = 0.25 * ( coord[0] + coord[4] + coord[7] + coord[3] );
                face_coord[2] = 0.25 * ( coord[0] + coord[1] + coord[5] + coord[4] );
                face_coord[3] = 0.25 * ( coord[4] + coord[5] + coord[6] + coord[7] );
                face_coord[4] = 0.25 * ( coord[1] + coord[2] + coord[6] + coord[5] );
                face_coord[5] = 0.25 * ( coord[2] + coord[3] + coord[7] + coord[6] );

                Real3 median1 = face_coord[0] - face_coord[3];
                Real3 median2 = face_coord[2] - face_coord[5];
                Real3 median3 = face_coord[1] - face_coord[4];
                Real d1 = median1.normL2();
                Real d2 = median2.normL2();
                Real d3 = median3.normL2();

                Real dx_numerator = d1 * d2 * d3;
                Real dx_denominator = d1 * d2 + d1 * d3 + d2 * d3;
                m_caracteristic_length[icell] = dx_numerator / dx_denominator;
            } else if ( options()->longueurCaracteristique() == "racine-cubique-volume" ) {
                m_caracteristic_length[icell] = std::pow ( m_cell_volume[icell], racine );
            } else if ( options()->longueurCaracteristique() == "monodimX" ) {
                Real x1=m_node_coord[cell.node (0)].x;
                Real x2=m_node_coord[cell.node (2)].x;
                m_caracteristic_length[icell] = math::abs(x2-x1);
            } else if ( options()->longueurCaracteristique() == "monodimY" ) {
                Real y1=m_node_coord[cell.node (0)].y;
                Real y2=m_node_coord[cell.node (2)].y;
                m_caracteristic_length[icell] = math::abs(y2-y1);
            } else if ( options()->longueurCaracteristique() == "monodimZ" ) {
                Real z1=m_node_coord[cell.node (0)].z;
                Real z2=m_node_coord[cell.node (2)].z;
                m_caracteristic_length[icell] = math::abs(z2-z1);
            } else {
                info() << " pas de longeur caractéritique definie dans le .arc " << options()->longueurCaracteristique();
                subDomain()->timeLoopMng()->stopComputeLoop ( true );
            }
        }

    }

    // maille mixte
    // moyenne sur la maille
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        if ( all_env_cell.nbEnvironment() !=1 ) {
            ENUMERATE_CELL_ENVCELL ( ienvcell,all_env_cell ) {
                EnvCell ev = *ienvcell;
                Cell cell = ev.globalCell();
                m_cell_volume[ev] = m_fracvol[ev] * m_cell_volume[cell];
            }
        }
    }
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
    debug() << " Rentrée dans updatedensity";
    if ( options()->sansLagrange ) {
        return;
    }
    debug() << my_rank << " : " << " Entree dans updateDensity() ";
    ENUMERATE_ENV ( ienv,mm ) {
        IMeshEnvironment* env = *ienv;
        Real density_thresold = options()->environment[env->id()].eosModel()->getdensityDamageThresold ( env );
        ENUMERATE_ENVCELL ( ienvcell,env ) {
            EnvCell ev = *ienvcell;
            Cell cell = ev.globalCell();
            // pinfo() << my_rank << " : " << cell.uniqueId() << " " << m_cell_volume[ev];
            Real new_density = m_cell_mass[ev] / m_cell_volume[ev];
            // nouvelle density
            m_density[ev] = new_density;
            // Options de blocage de la densité à une valeur seuil 
            if (m_density[ev]/m_density_0[cell] < density_thresold) {
              // pinfo()  << ev.globalCell().localId() << " blocage densité " << m_density[ev]  << " et " <<  m_density_0[cell];
              m_density[ev] = density_thresold * m_density_0[cell];
              // pinfo()  << ev.globalCell().localId() << " donne blocage densité " << m_density[ev]  << " et " <<  m_density_0[cell];
            }
            
            // volume specifique de l'environnement au temps n+1/2
            m_tau_density[ev] =
                0.5 * ( 1.0 / m_density_n[ev] + 1.0 / m_density[ev] );
        }
    }

    
    ENUMERATE_CELL ( icell,allCells() ) {
        Cell cell = * icell;
        Real new_density = m_cell_mass[icell] / m_cell_volume[icell];
        // nouvelle density
        m_density[icell] = new_density;
        // volume specifique moyen au temps n+1/2
        m_tau_density[icell] =
            0.5 * ( 1.0 / m_density_n[icell] + 1.0 / m_density[icell] );
        // divergence de la vitesse mode A1
        m_div_u[icell] =
            1.0 / m_global_deltat()  * ( 1.0 / m_density[icell] - 1.0 / m_density_n[icell] )
            / m_tau_density[icell];
    }

    m_density.synchronize();
    m_tau_density.synchronize();
    m_div_u.synchronize();
}
/**
 *******************************************************************************
 * \file updateElasticityAndPlasticity
 *  * \brief Calcul de la contribution elasto-plastique
 *
 * \param
 *
 * \return
 *******************************************************************************
 */
void MahycoModule::updateElasticityAndPlasticity()
{
    debug() << "rentrée updateElasticityAndPlasticity ";
    if ( options()->sansLagrange ) {
        return;
    }
    // Calcul du gradient de vitesse
    options()->environment[0].elastoModel()->ComputeVelocityGradient();
    // Calcul du tenseur de déformation et de rotation
    options()->environment[0].elastoModel()->ComputeDeformationAndRotation();

    ENUMERATE_ENV ( ienv,mm ) {
        IMeshEnvironment* env = *ienv;
        // Calcul du travail elasto-plastique pour energie interne avec S à n
        options()->environment[env->id()].elastoModel()->ComputeElastoEnergie ( env,m_global_deltat.value() );
        //  Calcul du deviateur elasto-plastique
        options()->environment[env->id()].elastoModel()->ComputeElasticity ( env,m_global_deltat.value(),m_dimension );
        //  Calcul de la plasticité
        options()->environment[env->id()].elastoModel()->ComputePlasticity ( env,m_global_deltat.value(),m_dimension );
        
        // TODO : appliquer l'enommagement au tenseur des contraintes
        
        // Calcul du travail elasto-plastique pour energie interne avec S à n+1
        options()->environment[env->id()].elastoModel()->ComputeElastoEnergie ( env,m_global_deltat.value() );
        
    }
}
/**
 *******************************************************************************
 * \file DepotEnergy(IMeshEnvironment* env) 
 *  * \brief Calcul du dépot d'énergie
 *
 * \param
 *
 * \return
 *******************************************************************************
 */
void MahycoModule::DepotEnergy(IMeshEnvironment* env) 
{
    // energyDepot EnDepot = options()->environment[env->id()]->energyDepot[0];
    Real value = options()->environment[env->id()]->energyDepot[0]->valeurSourceEnergie();
    TypesMahyco::eEnergyDepot type = options()->environment[env->id()]->energyDepot[0]->type();
    Real3 origine = options()->environment[env->id()]->energyDepot[0]->origine();
    int i=0;
    if (type == TypesMahyco::DepotConstant) {
        ENUMERATE_ENVCELL ( ienvcell,env ) {
          EnvCell ev = *ienvcell;
          Cell cell = ev.globalCell();
          value = options()->environment[env->id()]->energyDepot[0]->valeurSourceEnergie();
          if (m_global_time() < options()->environment[env->id()]->energyDepot[0]->Tdebut ) value = 0.;
          if (m_global_time() > options()->environment[env->id()]->energyDepot[0]->Tfin ) value = 0.;
          if (math::abs(m_cell_coord[cell].x - origine.x)  > options()->environment[env->id()]->energyDepot[0]->cutoffX) value=0.;
          if (math::abs(m_cell_coord[cell].y - origine.y)  > options()->environment[env->id()]->energyDepot[0]->cutoffY) value=0.;
          m_internal_energy[ev] += value * m_global_deltat();
        }
    } else if (type == TypesMahyco::DepotLineaire) {
        ENUMERATE_ENVCELL ( ienvcell,env ) {
          EnvCell ev = *ienvcell;
          Cell cell = ev.globalCell();
          value = options()->environment[env->id()]->energyDepot[0]->valeurSourceEnergie();
          if (m_global_time() < options()->environment[env->id()]->energyDepot[0]->Tdebut ) value = 0.;
          if (m_global_time() > options()->environment[env->id()]->energyDepot[0]->Tfin ) value = 0.;
          // E + DE/DX * X 
          value += options()->environment[env->id()]->energyDepot[0]->dependanceX * m_cell_coord[cell].x;
          value += options()->environment[env->id()]->energyDepot[0]->dependanceY * m_cell_coord[cell].y;
          value += options()->environment[env->id()]->energyDepot[0]->dependanceZ * m_cell_coord[cell].z;
          // E + DE/DT * T 
          value += options()->environment[env->id()]->energyDepot[0]->dependanceT * m_global_time();
          if (math::abs(m_cell_coord[cell].x - origine.x)  > options()->environment[env->id()]->energyDepot[0]->cutoffX) value=0.;
          if (math::abs(m_cell_coord[cell].y - origine.y)  > options()->environment[env->id()]->energyDepot[0]->cutoffY) value=0.;
          m_internal_energy[ev] += value * m_global_deltat();
        }
    }  else if (type == TypesMahyco::DepotSuperGaussian) {
        Real power = options()->environment[env->id()]->energyDepot[0]->power;
        ENUMERATE_ENVCELL ( ienvcell,env ) {
          EnvCell ev = *ienvcell;
          Cell cell = ev.globalCell();
          value = options()->environment[env->id()]->energyDepot[0]->valeurSourceEnergie();
          if (m_global_time() < options()->environment[env->id()]->energyDepot[0]->Tdebut ) value = 0.;
          if (m_global_time() > options()->environment[env->id()]->energyDepot[0]->Tfin ) value = 0.;
          Real ay= m_cell_coord[cell].y * options()->environment[env->id()]->energyDepot[0]->dependanceY;
          value *= math::exp ( -math::pow ( ay, power) );
          Real ax= m_cell_coord[cell].x * options()->environment[env->id()]->energyDepot[0]->dependanceX;
          value *= math::exp ( -math::pow ( ax, power) );
          if (math::abs(m_cell_coord[cell].x - origine.x)  > options()->environment[env->id()]->energyDepot[0]->cutoffX) value=0.;
          if (math::abs(m_cell_coord[cell].y - origine.y)  > options()->environment[env->id()]->energyDepot[0]->cutoffY) value=0.;
          m_internal_energy[ev] += value * m_global_deltat();

        }
    }
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
    if ( options()->sansLagrange ) {
        return;
    }
    debug() << " Rentrée dans updateEnergyAndPressure";
    ENUMERATE_ENV ( ienv,mm ) {
        IMeshEnvironment* env = *ienv;
        if (options()->environment[env->id()]->energyDepot.size() != 0.) DepotEnergy(env);
    }
    if ( options()->pressionExplicite ) {
        updateEnergyAndPressureExplicite();
    } else {
        if ( options()->withNewton ) {
            updateEnergyAndPressurebyNewton();
        } else {
            updateEnergyAndPressureforGP();
        }
    }

    // maille mixte
    // moyenne sur la maille
    CellToAllEnvCellConverter all_env_cell_converter ( mm );
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        if ( all_env_cell.nbEnvironment() !=1 ) {
            m_internal_energy[cell] = 0.;
            ENUMERATE_CELL_ENVCELL ( ienvcell,all_env_cell ) {
                EnvCell ev = *ienvcell;
                m_internal_energy[cell] += m_mass_fraction[ev] * m_internal_energy[ev];
            }
        }
    }
    if ( ! options()->withProjection ) {
        // Calcul de la Pression si on ne fait pas de projection
        // Si on n'a pas la pression en explicite, on refait une derniere iteration 
        // pour calculer la pression à n+1
        if ( ! options()->pressionExplicite )
            for ( Integer i=0,n=options()->environment().size(); i<n; ++i ) {
                IMeshEnvironment* ienv = mm->environments() [i];
                // Calcul de la pression et de la vitesse du son
                options()->environment[i].eosModel()->applyEOS ( ienv );
                // Test endommagement --> Pression devient nulle ?
                options()->environment[i].eosModel()->Endommagement ( ienv );
            }
        // Calcul de la Pression moyenne    
        computePressionMoyenne();
    }
}
/*
 *******************************************************************************
*/
void MahycoModule::updateEnergyAndPressureExplicite()
{

    if ( options()->sansLagrange ) {
        return;
    }
    ENUMERATE_ENV ( ienv,mm ) {
        IMeshEnvironment* env = *ienv;
        ENUMERATE_ENVCELL ( ienvcell,env ) {
            EnvCell ev = *ienvcell;
            m_internal_energy[ev] -=
                ( m_pressure_n[ev] + m_pseudo_viscosity[ev] ) * ( 1.0 / m_density[ev] - 1.0 / m_density_n[ev] );
        }
        // puis on récupere la pression pour toutes les mailles
        options()->environment[env->id()].eosModel()->applyEOS ( env );
        // Test endommagement --> Pression devient nulle ?
        options()->environment[env->id()].eosModel()->Endommagement ( env );
    }
}
/*
 *******************************************************************************
*/
void MahycoModule::updateEnergyAndPressurebyNewton()
{

    if ( options()->sansLagrange ) {
        return;
    }
    debug() << " Entree dans updateEnergyAndPressure()";
    bool csts = options()->schemaCsts();
    bool pseudo_centree = options()->pseudoCentree();
    // Calcul de l'énergie interne
    if ( !csts ) {
        ENUMERATE_ENV ( ienv,mm ) {
            IMeshEnvironment* env = *ienv;
            Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst ( env );
            Real tension_limit = options()->environment[env->id()].eosModel()->getTensionLimitCst ( env );
            ENUMERATE_ENVCELL ( ienvcell,env ) {
                EnvCell ev = *ienvcell;
                Real pseudo ( 0. );
                if ( pseudo_centree &&
                        ( ( m_pseudo_viscosity_n[ev] + m_pseudo_viscosity[ev] ) * ( 1.0 / m_density[ev] - 1.0 / m_density_n[ev] ) < 0. ) ) {
                    pseudo = 0.5 * ( m_pseudo_viscosity[ev] + m_pseudo_viscosity_n[ev] );
                }
                if ( !pseudo_centree &&
                        ( m_pseudo_viscosity[ev] * ( 1.0 / m_density[ev] - 1.0 / m_density_n[ev] ) < 0. ) ) {
                    pseudo = m_pseudo_viscosity[ev];
                }

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
                double enew=0, e=m_internal_energy[ev], p=pn, c, dpde;
                int i = 0;

                // TODO : sauvegarde des indicateurs de phases à N

                while ( i<itermax && abs ( fvnr ( e, p, dpde, en, qnn1, pn, rn1, rn ) ) >=epsilon ) {
                    // TODO : reprise des indicateurs de phases à N
                    m_internal_energy[ev] = e;
                    options()->environment[env->id()].eosModel()->applyOneCellEOS ( env, ev );
                    p = m_pressure[ev];
                    c = m_sound_speed[ev];
                    dpde = m_dpde[ev];
                    enew = e - fvnr ( e, p, dpde, en, qnn1, pn, rn1, rn ) / fvnrderiv ( e, dpde, rn1, rn );
                    e = enew;
                    i++;
                }
                m_internal_energy[ev] = e;
                m_sound_speed[ev] = c;
                m_pressure[ev] = p;
            }
        }
    } else {
        ENUMERATE_ENV ( ienv,mm ) {
            IMeshEnvironment* env = *ienv;
            Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst ( env );
            Real tension_limit = options()->environment[env->id()].eosModel()->getTensionLimitCst ( env );
            ENUMERATE_ENVCELL ( ienvcell,env ) {
                EnvCell ev = *ienvcell;
                Cell cell=ev.globalCell();
                Real cqs_v_nplus1 ( 0. );
                Real cqs_v_n ( 0. );
                Real cqs_delta_v ( 0. );
                Real cqs_v_old_n ( 0. );
                for ( Integer inode = 0; inode < cell.nbNode(); ++inode ) {
                    cqs_v_nplus1 += math::dot ( m_velocity[cell.node ( inode )], m_cell_cqs[cell] [inode] )
                                    * m_global_deltat();
                    cqs_v_n += math::dot ( m_velocity[cell.node ( inode )], m_cell_cqs_n[cell] [inode] )
                               * m_global_deltat();
                    cqs_delta_v +=  math::dot ( m_velocity[cell.node ( inode )] - m_velocity_n[cell.node ( inode )], m_cell_cqs_n[cell] [inode] )
                                    * ( m_global_deltat() - m_global_old_deltat() );
                    cqs_v_old_n += math::dot ( m_velocity_n[cell.node ( inode )], m_cell_cqs_n[cell] [inode] )
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


                // TODO : sauvegarde des indicateurs de phases à N

                while ( i<itermax && abs ( f ( e, p, dpde, en, qn, pn, cn1, cn, m, qn1, cdn, cdon, qnm1 ) ) >=epsilon ) {
                    // TODO : reprise des indicateurs de phases à N
                    m_internal_energy[ev] = e;
                    options()->environment[env->id()].eosModel()->applyOneCellEOS ( env, ev );
                    p = m_pressure[ev];
                    c = m_sound_speed[ev];
                    dpde = m_dpde[ev];
                    enew = e - f ( e, p, dpde, en, qn, pn, cn1, cn, m, qn1, cdn, cdon, qnm1 ) / fderiv ( e, p, dpde, cn1, m );
                    e = enew;
                    i++;
                }
                m_internal_energy[ev] = e;
                m_sound_speed[ev] = c;
                m_pressure[ev] = p;
            }
        }
    }
}

/*
 *******************************************************************************
 */
void MahycoModule::
updateEnergyAndPressureforGP()
{
    if ( options()->sansLagrange ) {
        return;
    }
    debug() << " Entree dans updateEnergyAndPressure()";
    bool csts = options()->schemaCsts();
    bool pseudo_centree = options()->pseudoCentree();
    // Calcul de l'énergie interne
    if ( !csts ) {
        ENUMERATE_ENV ( ienv,mm ) {
            IMeshEnvironment* env = *ienv;
            Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst ( env );
            ENUMERATE_ENVCELL ( ienvcell,env ) {
                EnvCell ev = *ienvcell;
                Real pseudo ( 0. );
                Cell cell = ev.globalCell();
                if ( pseudo_centree &&
                        ( ( m_pseudo_viscosity_n[ev] + m_pseudo_viscosity[ev] ) * ( 1.0 / m_density[ev] - 1.0 / m_density_n[ev] ) < 0. ) ) {
                    pseudo = 0.5 * ( m_pseudo_viscosity[ev] + m_pseudo_viscosity_n[ev] );
                }
                if ( !pseudo_centree &&
                        ( m_pseudo_viscosity[ev] * ( 1.0 / m_density[ev] - 1.0 / m_density_n[ev] ) < 0. ) ) {
                    pseudo = m_pseudo_viscosity[ev];
                }

                Real denom_accrois_nrj ( 1 + 0.5 * ( adiabatic_cst - 1.0 ) *
                                         m_density[ev] *
                                         ( 1.0 / m_density[ev] -
                                           1.0 / m_density_n[ev] ) );
                Real numer_accrois_nrj ( m_internal_energy_n[ev] -
                                         ( 0.5 * m_pressure[ev] + pseudo ) *
                                         ( 1.0 / m_density[ev] - 1.0 / m_density_n[ev] ) );
                m_internal_energy[ev] = numer_accrois_nrj / denom_accrois_nrj;
            }
        }
    } else {
        ENUMERATE_ENV ( ienv,mm ) {
            IMeshEnvironment* env = *ienv;
            Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst ( env );
            ENUMERATE_ENVCELL ( ienvcell,env ) {
                EnvCell ev = *ienvcell;
                Cell cell = ev.globalCell();
                Real cqs_v_nplus1 ( 0. );
                Real cqs_v_n ( 0. );
                Real cqs_delta_v ( 0. );
                Real cqs_v_old_n ( 0. );
                for ( Integer inode = 0; inode < cell.nbNode(); ++inode ) {
                    cqs_v_nplus1 += math::dot ( m_velocity[cell.node ( inode )], m_cell_cqs[cell] [inode] )
                                    * m_global_deltat();
                    cqs_v_n += math::dot ( m_velocity[cell.node ( inode )], m_cell_cqs_n[cell] [inode] )
                               * m_global_deltat();
                    cqs_delta_v +=  math::dot ( m_velocity[cell.node ( inode )] - m_velocity_n[cell.node ( inode )], m_cell_cqs_n[cell] [inode] )
                                    * ( m_global_deltat() - m_global_old_deltat() );
                    cqs_v_old_n += math::dot ( m_velocity_n[cell.node ( inode )], m_cell_cqs_n[cell] [inode] )
                                   * m_global_old_deltat();
                }
                Real denom_accrois_nrj ( 1 + 0.5 * ( adiabatic_cst - 1.0 )
                                         * m_density[ev] * cqs_v_nplus1 / m_cell_mass[ev] ) ;
                Real numer_accrois_nrj ( m_internal_energy[ev]
                                         - ( 0.5 * ( m_pressure_n[ev] + m_pseudo_viscosity_n[ev] )
                                             * cqs_v_n / m_cell_mass[ev] )
                                         - ( 0.5 * m_pseudo_viscosity[ev]
                                             * cqs_v_nplus1 / m_cell_mass[ev] )
                                         + ( 0.25 * ( m_pressure_n[ev] + m_pseudo_viscosity_n[ev] )
                                             * cqs_delta_v / m_cell_mass[ev] )
                                         - ( 0.5 * ( m_pseudo_viscosity_n[ev] - m_pseudo_viscosity_nmoins1[ev] )
                                             * cqs_v_old_n / m_cell_mass[ev] )
                                       );

                m_internal_energy[ev] = numer_accrois_nrj / denom_accrois_nrj;
            }
        }
    }
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
    if ( options()->sansLagrange ) {
        return;
    }
    Real option = 0.;
    if ( options()->calculVitsonMax ) option = 1.;
    debug() << " Entree dans computePressionMoyenne() ";
    // maille mixte
    // moyenne sur la maille
    CellToAllEnvCellConverter all_env_cell_converter ( mm );
    ENUMERATE_CELL ( icell, allCells() ) {
        Cell cell = * icell;
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        if ( all_env_cell.nbEnvironment() !=1 ) {
            m_pressure[cell] = 0.;
            // m_sound_speed[icell] = 1.e-20;
            m_sound_speed[icell] = 0.;
            ENUMERATE_CELL_ENVCELL ( ienvcell,all_env_cell ) {
                EnvCell ev = *ienvcell;
                m_pressure[cell] += m_fracvol[ev] * m_pressure[ev];
                m_sound_speed[cell] = option * std::max ( m_sound_speed[ev], m_sound_speed[cell] );
                m_sound_speed[cell] += (1.-option) * m_fracvol[ev] * m_sound_speed[ev];
            }
        }
    }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
computeDeltaT()
{
    debug() << " Entree dans compute DT avec " << m_global_old_deltat()
            << " et " << options()->deltatInit()
            << " et " << m_global_deltat()
            << " et " << options()->deltatMax();

    m_global_old_deltat = m_global_deltat;
    m_old_deltat = m_global_old_deltat();

    Real new_dt = FloatInfo < Real >::maxValue();
    if ( options()->sansLagrange || options()->deltatConstant) {
        // on garde le meme pas de temps
        new_dt = options()->deltatInit();
    } else {
        CellToAllEnvCellConverter all_env_cell_converter ( mm );

        // Calcul du pas de temps pour le respect du critère de CFL
        Real minimum_aux = FloatInfo < Real >::maxValue();
        Int64 cell_id ( -1 ), nbenvcell ( -1 );
        Real cc ( 0. ), ll ( 0. );
        Real3 coord;
        Cell cell_min;
        Real vit_cell_min(0.);
        Real new_dt_proc;

        ENUMERATE_CELL ( icell, allCells() ) {
            Cell cell = * icell;
            Real cell_dx = m_caracteristic_length[icell];
            Real sound_speed = m_sound_speed[icell];
            Real vmax ( 0. );
            if ( /* options()->withProjection */ options()->remap()->isEuler() == "Full")
                for ( NodeEnumerator inode ( cell.nodes() ); inode.index() < cell.nbNode(); ++inode ) {
                    vmax = math::max ( m_velocity[inode].normL2(), vmax );
                }
            Real dx_sound = cell_dx / ( sound_speed + vmax );
            m_deltat_cell[icell] = options()->cfl() * dx_sound;
            minimum_aux = math::min ( minimum_aux, dx_sound );
            if ( minimum_aux == dx_sound ) {
                cell_id = cell.uniqueId();
                cc = m_sound_speed[icell];
                ll = m_caracteristic_length[icell];
                coord = m_cell_coord[icell];
                AllEnvCell all_env_cell = all_env_cell_converter[cell];
                nbenvcell = all_env_cell.nbEnvironment();
                cell_min = cell;
                vit_cell_min = vmax;
            }
        }
        
        new_dt = options()->cfl() * minimum_aux;
        new_dt_proc = new_dt;
        new_dt = parallelMng()->reduce ( Parallel::ReduceMin, new_dt );
        
        if (new_dt < .5*m_global_old_deltat() ) {
            pinfo() << " Chute du pas de temps " << new_dt;
            pinfo() << " nouveau pas de temps " << new_dt_proc << " par " << cell_id << " position : " << coord;
            pinfo() << " (avec " << nbenvcell << " envs) avec vitson = " << cc << " et longueur  =  " << ll << " et min " << minimum_aux;
            pinfo() << " Vitesse matérielle " << vit_cell_min;
            pinfo() << " Pseudo " << m_pseudo_viscosity[cell_min];
            AllEnvCell all_env_cell = all_env_cell_converter[cell_min];
             ENUMERATE_CELL_ENVCELL ( ienvcell,all_env_cell ) {
                EnvCell ev = *ienvcell;
                int index_env = ev.environmentId(); 
                pinfo() << " Vitson : env  " << index_env << " = " << m_sound_speed[ev] << " et pression= " << m_pressure[ev];
                pinfo() << " pour une fraction de " << m_fracvol[ev] << " et une Pseudo " << m_pseudo_viscosity[ev]; 
                pinfo() << " et une densité de " << m_density[ev] << " et avant densité de "  << m_density_n[ev];
             }
        }
        // respect de taux de croissance max
        new_dt = math::min ( new_dt, 1.05 * m_global_old_deltat() );
        // respect de la valeur max imposée par le fichier de données .plt
        debug() << " nouveau pas de temps " << new_dt << " par " << cell_id <<
                " (avec " << nbenvcell << " envs) avec " << cc << " " << ll << " et min " << minimum_aux;
        new_dt = math::min ( new_dt, options()->deltatMax() );
        // respect du pas de temps minimum
        if ( new_dt < options()->deltatMin() ) {
            pinfo() << " ------------------------------------------" ;
            pinfo() << " pas de temps minimum " << new_dt ;
            pinfo() << " nouveau pas de temps " << new_dt_proc << " par " << cell_id << " position : " << coord;
            pinfo() << " (avec " << nbenvcell << " envs) avec vitson = " << cc << " et longueur  =  " << ll << " et min " << minimum_aux;
            pinfo() << " Vitesse matérielle " << vit_cell_min;
            pinfo() << " Pseudo " << m_pseudo_viscosity[cell_min];
            AllEnvCell all_env_cell = all_env_cell_converter[cell_min];
            ENUMERATE_CELL_ENVCELL ( ienvcell,all_env_cell ) {
                EnvCell ev = *ienvcell;
                int index_env = ev.environmentId(); 
                pinfo() << " Vitson : env  " << index_env << " = " <<  m_sound_speed[ev] << " et pression= " << m_pressure[ev];
                pinfo() << " pour une fraction de " << m_fracvol[ev] << " et une Pseudo " << m_pseudo_viscosity[ev];
                pinfo() << " et une densité de " << m_density[ev] << " et avant densité de "  << m_density_n[ev];
            }
            pinfo() << " ------------------------------------------" ;
            exit ( 1 );
        }

        
        debug() << " nouveau pas de temps2 " << new_dt;

    }

    // Mise à jour du pas de temps
    m_global_deltat = new_dt;

    // Le dernier calcul se fait exactement au temps stopTime()
    Real stop_time  = options()->finalTime();
    bool not_yet_finish = ( m_global_time() < stop_time );
    bool finish = ( m_global_time() > stop_time );
    bool too_much = ( ( m_global_time()+new_dt ) > stop_time );


    debug() << " nouveau pas de temps " << new_dt;
    if ( not_yet_finish && too_much ) {
        new_dt = stop_time - m_global_time();
        subDomain()->timeLoopMng()->stopComputeLoop ( true );
    }
    if ( finish ) {
        subDomain()->timeLoopMng()->stopComputeLoop ( true );
    }

    debug() << " time " << m_global_time() << " et fin à " << stop_time;
    debug() << " not_yet_finish " << not_yet_finish;
    debug() << " finish " << finish;
    debug() << " too_much " << too_much;

}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
inline void MahycoModule::
computeCQs ( Real3 node_coord[8], Real3 face_coord[6], const Cell & cell )
{
    const Real3 c0 = face_coord[0];
    const Real3 c1 = face_coord[1];
    const Real3 c2 = face_coord[2];
    const Real3 c3 = face_coord[3];
    const Real3 c4 = face_coord[4];
    const Real3 c5 = face_coord[5];

    // Calcul des normales face 1 :
    const Real3 n1a04 = 0.5 * math::vecMul ( node_coord[0] - c0, node_coord[3] - c0 );
    const Real3 n1a03 = 0.5 * math::vecMul ( node_coord[3] - c0, node_coord[2] - c0 );
    const Real3 n1a02 = 0.5 * math::vecMul ( node_coord[2] - c0, node_coord[1] - c0 );
    const Real3 n1a01 = 0.5 * math::vecMul ( node_coord[1] - c0, node_coord[0] - c0 );

    // Calcul des normales face 2 :
    const Real3 n2a05 = 0.5 * math::vecMul ( node_coord[0] - c1, node_coord[4] - c1 );
    const Real3 n2a12 = 0.5 * math::vecMul ( node_coord[4] - c1, node_coord[7] - c1 );
    const Real3 n2a08 = 0.5 * math::vecMul ( node_coord[7] - c1, node_coord[3] - c1 );
    const Real3 n2a04 = 0.5 * math::vecMul ( node_coord[3] - c1, node_coord[0] - c1 );

    // Calcul des normales face 3 :
    const Real3 n3a01 = 0.5 * math::vecMul ( node_coord[0] - c2, node_coord[1] - c2 );
    const Real3 n3a06 = 0.5 * math::vecMul ( node_coord[1] - c2, node_coord[5] - c2 );
    const Real3 n3a09 = 0.5 * math::vecMul ( node_coord[5] - c2, node_coord[4] - c2 );
    const Real3 n3a05 = 0.5 * math::vecMul ( node_coord[4] - c2, node_coord[0] - c2 );

    // Calcul des normales face 4 :
    const Real3 n4a09 = 0.5 * math::vecMul ( node_coord[4] - c3, node_coord[5] - c3 );
    const Real3 n4a10 = 0.5 * math::vecMul ( node_coord[5] - c3, node_coord[6] - c3 );
    const Real3 n4a11 = 0.5 * math::vecMul ( node_coord[6] - c3, node_coord[7] - c3 );
    const Real3 n4a12 = 0.5 * math::vecMul ( node_coord[7] - c3, node_coord[4] - c3 );

    // Calcul des normales face 5 :
    const Real3 n5a02 = 0.5 * math::vecMul ( node_coord[1] - c4, node_coord[2] - c4 );
    const Real3 n5a07 = 0.5 * math::vecMul ( node_coord[2] - c4, node_coord[6] - c4 );
    const Real3 n5a10 = 0.5 * math::vecMul ( node_coord[6] - c4, node_coord[5] - c4 );
    const Real3 n5a06 = 0.5 * math::vecMul ( node_coord[5] - c4, node_coord[1] - c4 );

    // Calcul des normales face 6 :
    const Real3 n6a03 = 0.5 * math::vecMul ( node_coord[2] - c5, node_coord[3] - c5 );
    const Real3 n6a08 = 0.5 * math::vecMul ( node_coord[3] - c5, node_coord[7] - c5 );
    const Real3 n6a11 = 0.5 * math::vecMul ( node_coord[7] - c5, node_coord[6] - c5 );
    const Real3 n6a07 = 0.5 * math::vecMul ( node_coord[6] - c5, node_coord[2] - c5 );

    // Calcul des résultantes aux sommets :
    m_cell_cqs[cell] [0] = ( 5. * ( n1a01 + n1a04 + n2a04 + n2a05 + n3a05 + n3a01 ) +
                             ( n1a02 + n1a03 + n2a08 + n2a12 + n3a06 + n3a09 ) ) * ( 1. / 12. );
    m_cell_cqs[cell] [1] = ( 5. * ( n1a01 + n1a02 + n3a01 + n3a06 + n5a06 + n5a02 ) +
                             ( n1a04 + n1a03 + n3a09 + n3a05 + n5a10 + n5a07 ) ) * ( 1. / 12. );
    m_cell_cqs[cell] [2] = ( 5. * ( n1a02 + n1a03 + n5a07 + n5a02 + n6a07 + n6a03 ) +
                             ( n1a01 + n1a04 + n5a06 + n5a10 + n6a11 + n6a08 ) ) * ( 1. / 12. );
    m_cell_cqs[cell] [3] = ( 5. * ( n1a03 + n1a04 + n2a08 + n2a04 + n6a08 + n6a03 ) +
                             ( n1a01 + n1a02 + n2a05 + n2a12 + n6a07 + n6a11 ) ) * ( 1. / 12. );
    m_cell_cqs[cell] [4] = ( 5. * ( n2a05 + n2a12 + n3a05 + n3a09 + n4a09 + n4a12 ) +
                             ( n2a08 + n2a04 + n3a01 + n3a06 + n4a10 + n4a11 ) ) * ( 1. / 12. );
    m_cell_cqs[cell] [5] = ( 5. * ( n3a06 + n3a09 + n4a09 + n4a10 + n5a10 + n5a06 ) +
                             ( n3a01 + n3a05 + n4a12 + n4a11 + n5a07 + n5a02 ) ) * ( 1. / 12. );
    m_cell_cqs[cell] [6] = ( 5. * ( n4a11 + n4a10 + n5a10 + n5a07 + n6a07 + n6a11 ) +
                             ( n4a12 + n4a09 + n5a06 + n5a02 + n6a03 + n6a08 ) ) * ( 1. / 12. );
    m_cell_cqs[cell] [7] = ( 5. * ( n2a08 + n2a12 + n4a12 + n4a11 + n6a11 + n6a08 ) +
                             ( n2a04 + n2a05 + n4a09 + n4a10 + n6a07 + n6a03 ) ) * ( 1. / 12. );
}

/*---------------------------------------------------------------------------*/
Real MahycoModule::produit ( Real A, Real B, Real C, Real D )
{
    return ( A*B-C*D );
}
/*---------------------------------------------------------------------------*/
inline Real MahycoModule::norme ( Real E, Real F, Real G )
{
    return ( sqrt ( ( E*E )+ ( F*F )+ ( G*G ) ) );
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_MAHYCO ( MahycoModule );
