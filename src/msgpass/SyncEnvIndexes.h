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
    AccMemAdviser* acc_mem_adv);
  virtual ~SyncEnvIndexes() {}

  //! Recalcule les listes de EnvIndex "owned" ("shared") et "ghost"
  void updateEnvIndexes();

  MatVarSpace space() const { return m_mvs; }

  auto nbOwnedEviPn() const {
    return m_nb_owned_evi_pn.constView();
  }

  auto nbGhostEviPn() const {
    return m_nb_ghost_evi_pn.constView();
  }

  // Listes par voisin

  auto ownedEviPn() const {
    return ConstMultiArray2View<EnvVarIndex>(m_buf_owned_evi.constView(),
        m_indexes_owned_evi_pn.constView(), m_nb_owned_evi_pn.constView());
  }

  auto ghostEviPn() const {
    return ConstMultiArray2View<EnvVarIndex>(m_buf_ghost_evi.constView(),
        m_indexes_ghost_evi_pn.constView(), m_nb_ghost_evi_pn.constView());
  }

  // Listes par environnement

  auto allEviPenv() const {
    return ConstMultiArray2View<EnvVarIndex>(m_buf_all_evi_4env.constView(),
        m_indexes_all_evi_penv.constView(), m_nb_all_evi_penv.constView());
  }

  auto ownEviPenv() const {
    return ConstMultiArray2View<EnvVarIndex>(m_buf_owned_evi_4env.constView(),
        m_indexes_owned_evi_penv.constView(), m_nb_owned_evi_penv.constView());
  }

  auto privateEviPenv() const {
    return ConstMultiArray2View<EnvVarIndex>(m_buf_private_evi_4env.constView(),
        m_indexes_private_evi_penv.constView(), m_nb_private_evi_penv.constView());
  }

  auto sharedEviPenv() const {
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

  // "shared" ou "owned" : les EnvVarIndex(es) intérieurs au sous-domaine et qui doivent être envoyés
  // "ghost" : les EnvVarIndex(es) fantômes pour lesquels on va recevoir des informations
  // _pn : _per_neigh, info par voisin
  
  // Nb d'EnvVarIndex(es) par voisin et listes des EnvVarIndex(es) mêmes
  // Owned ou Shared
  UniqueArray<EnvVarIndex> m_buf_owned_evi;
  IntegerUniqueArray       m_indexes_owned_evi_pn;
  IntegerUniqueArray       m_nb_owned_evi_pn;

  // Ghost
  UniqueArray<EnvVarIndex> m_buf_ghost_evi;
  IntegerUniqueArray       m_indexes_ghost_evi_pn;
  IntegerUniqueArray       m_nb_ghost_evi_pn;

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
