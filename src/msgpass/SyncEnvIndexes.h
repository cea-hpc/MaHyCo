#ifndef MSG_PASS_SYNC_ENV_INDEXES_H
#define MSG_PASS_SYNC_ENV_INDEXES_H

#include "accenv/AcceleratorUtils.h"
#include "accenv/MultiEnvUtils.h"
#include "msgpass/SyncItems.h"
#include <arcane/materials/IMeshMaterialMng.h>
#include <arcane/materials/IMeshMaterialVariableSynchronizer.h>
#include <arcane/utils/MultiArray2.h>

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/* Encapsule la liste des EnvVarIndex(es) à envoyer/recevoir                 */
/*---------------------------------------------------------------------------*/
class SyncEnvIndexes {
 public:
  SyncEnvIndexes(MatVarSpace mvs, IMeshMaterialMng* mm,
    Int32ConstArrayView neigh_ranks, 
    SyncItems<Cell>* sync_cells,
    AccMemAdviser* acc_mem_adv,
    bool evlist_per_env_needed);
  virtual ~SyncEnvIndexes() {}

  //! Recalcule les listes de EnvIndex "owned" ("shared") et "ghost"
  void updateEnvIndexes();

  MatVarSpace space() const { return m_mvs; }

  auto nbOwnedMviPn() const {
    return m_nb_owned_mvi_pn.constView();
  }

  auto nbGhostMviPn() const {
    return m_nb_ghost_mvi_pn.constView();
  }

  // Listes par voisin

  auto ownedMviPn() const {
    return ConstMultiArray2View<MatVarIndex>(m_buf_owned_mvi.constView(),
        m_indexes_owned_mvi_pn.constView(), m_nb_owned_mvi_pn.constView());
  }

  auto ghostMviPn() const {
    return ConstMultiArray2View<MatVarIndex>(m_buf_ghost_mvi.constView(),
        m_indexes_ghost_mvi_pn.constView(), m_nb_ghost_mvi_pn.constView());
  }
  // Listes par environnement

  auto allEviPenv() const {
    if (!m_evlist_per_env_needed) throw FatalErrorException(A_FUNCINFO, "Liste des EnvCell par environnement doit être activée");
    return ConstMultiArray2View<EnvVarIndex>(m_buf_all_evi_4env.constView(),
        m_indexes_all_evi_penv.constView(), m_nb_all_evi_penv.constView());
  }

  auto ownEviPenv() const {
    if (!m_evlist_per_env_needed) throw FatalErrorException(A_FUNCINFO, "Liste des EnvCell par environnement doit être activée");
    return ConstMultiArray2View<EnvVarIndex>(m_buf_owned_evi_4env.constView(),
        m_indexes_owned_evi_penv.constView(), m_nb_owned_evi_penv.constView());
  }

  auto privateEviPenv() const {
    if (!m_evlist_per_env_needed) throw FatalErrorException(A_FUNCINFO, "Liste des EnvCell par environnement doit être activée");
    return ConstMultiArray2View<EnvVarIndex>(m_buf_private_evi_4env.constView(),
        m_indexes_private_evi_penv.constView(), m_nb_private_evi_penv.constView());
  }

  auto sharedEviPenv() const {
    if (!m_evlist_per_env_needed) throw FatalErrorException(A_FUNCINFO, "Liste des EnvCell par environnement doit être activée");
    return ConstMultiArray2View<EnvVarIndex>(m_buf_shared_evi_4env.constView(),
        m_indexes_shared_evi_penv.constView(), m_nb_shared_evi_penv.constView());
  }
 
 protected:

  //! Calcule les listes de EnvIndex(es) par environnement à partir d'un groupe de mailles
  void _eviListFromGroup(CellGroup cell_group,
      UniqueArray<EnvVarIndex>& buf_evi, 
      IntegerUniqueArray&      indexes_evi_penv,
      IntegerUniqueArray&      nb_evi_penv
      );

 protected:

  IMeshMaterialMng*                  m_mesh_material_mng=nullptr;
  AccMemAdviser*                     m_acc_mem_adv=nullptr;
  IMeshMaterialVariableSynchronizer* m_mmvs=nullptr;
  SyncItems<Cell>*                   m_sync_cells=nullptr;
  MatVarSpace                        m_mvs;
  Integer                            m_nb_nei;
  bool                               m_evlist_per_env_needed=false;

  // "shared" ou "owned" : les MatVarIndex(es) intérieurs au sous-domaine et qui doivent être envoyés
  // "ghost" : les MatVarIndex(es) fantômes pour lesquels on va recevoir des informations
  // _pn : _per_neigh, info par voisin
  
  // Nb d'MatVarIndex(es) par voisin et listes des MatVarIndex(es) mêmes
  // Owned ou Shared
  UniqueArray<MatVarIndex> m_buf_owned_mvi;
  IntegerUniqueArray       m_indexes_owned_mvi_pn;
  IntegerUniqueArray       m_nb_owned_mvi_pn;

  // Ghost
  UniqueArray<MatVarIndex> m_buf_ghost_mvi;
  IntegerUniqueArray       m_indexes_ghost_mvi_pn;
  IntegerUniqueArray       m_nb_ghost_mvi_pn;

  // Listes des EnvVarIndex(es) par environnement des groupes de mailles ...
  // "all"
  UniqueArray<EnvVarIndex> m_buf_all_evi_4env    ;
  IntegerUniqueArray       m_indexes_all_evi_penv;
  IntegerUniqueArray       m_nb_all_evi_penv     ;

  // "owned"
  UniqueArray<EnvVarIndex> m_buf_owned_evi_4env    ;
  IntegerUniqueArray       m_indexes_owned_evi_penv;
  IntegerUniqueArray       m_nb_owned_evi_penv     ;

  // "private"
  UniqueArray<EnvVarIndex> m_buf_private_evi_4env    ;
  IntegerUniqueArray       m_indexes_private_evi_penv;
  IntegerUniqueArray       m_nb_private_evi_penv     ;

  // "shared"
  UniqueArray<EnvVarIndex> m_buf_shared_evi_4env    ;
  IntegerUniqueArray       m_indexes_shared_evi_penv;
  IntegerUniqueArray       m_nb_shared_evi_penv;
};

#endif
