#ifndef MSG_PASS_SYNC_ENV_INDEXES_H
#define MSG_PASS_SYNC_ENV_INDEXES_H

#include "accenv/AcceleratorUtils.h"
#include "accenv/MultiEnvUtils.h"
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
    Int32ConstArrayView neigh_ranks, AccMemAdviser* acc_mem_adv);
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

  auto ownedEviPn() const {
    return ConstMultiArray2View<EnvVarIndex>(m_buf_owned_evi.constView(),
        m_indexes_owned_evi_pn.constView(), m_nb_owned_evi_pn.constView());
  }

  auto ghostEviPn() const {
    return ConstMultiArray2View<EnvVarIndex>(m_buf_ghost_evi.constView(),
        m_indexes_ghost_evi_pn.constView(), m_nb_ghost_evi_pn.constView());
  }

 protected:

  IMeshMaterialMng*                  m_mesh_material_mng=nullptr;
  AccMemAdviser*                     m_acc_mem_adv=nullptr;
  IMeshMaterialVariableSynchronizer* m_mmvs=nullptr;
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
};

#endif
