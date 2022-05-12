#include "msgpass/SyncEnvIndexes.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
SyncEnvIndexes::SyncEnvIndexes(MatVarSpace mvs, IMeshMaterialMng* mm,
    Int32ConstArrayView neigh_ranks, AccMemAdviser* acc_mem_adv) :
  m_mesh_material_mng    (mm),
  m_acc_mem_adv          (acc_mem_adv),
  m_mvs                  (mvs),
  m_nb_nei               (neigh_ranks.size()),
  m_buf_owned_evi        (platform::getAcceleratorHostMemoryAllocator()),
  m_indexes_owned_evi_pn (platform::getAcceleratorHostMemoryAllocator()),
  m_nb_owned_evi_pn      (platform::getAcceleratorHostMemoryAllocator()),
  m_buf_ghost_evi        (platform::getAcceleratorHostMemoryAllocator()),
  m_indexes_ghost_evi_pn (platform::getAcceleratorHostMemoryAllocator()),
  m_nb_ghost_evi_pn      (platform::getAcceleratorHostMemoryAllocator())
{
  if (!m_mesh_material_mng) {
    throw NotSupportedException(A_FUNCINFO, "Unsupported m_mesh_material_mng nullptr");
  }

  if (m_mvs == MatVarSpace::MaterialAndEnvironment) {
    m_mmvs = m_mesh_material_mng->_allCellsMatEnvSynchronizer();
  } else if (m_mvs == MatVarSpace::Environment) {
    m_mmvs = m_mesh_material_mng->_allCellsEnvOnlySynchronizer();
  } else {
    throw NotSupportedException(A_FUNCINFO,
        String::format("Invalid space ={0}",(int)m_mvs));
  }
}

/*---------------------------------------------------------------------------*/
/* Recalcule les listes de EnvIndex "owned" ("shared") et "ghost"            */
/*---------------------------------------------------------------------------*/
void SyncEnvIndexes::updateEnvIndexes() {

  m_mmvs->checkRecompute();
  // TODO : Il faudrait pouvoir récupérer si ça a été recalculé et faire la suite
  // uniquement dans ce cas-là

  // "shared" ou "owned" : les EnvIndex(es) intérieurs au sous-domaine et qui doivent être envoyés
  // "ghost" : les EnvIndex(es) fantômes pour lesquels on va recevoir des informations
  m_indexes_owned_evi_pn.resize(m_nb_nei);
  m_indexes_ghost_evi_pn.resize(m_nb_nei);
  m_nb_owned_evi_pn.resize(m_nb_nei);
  m_nb_ghost_evi_pn.resize(m_nb_nei);

  Integer accu_nb_owned=0;
  Integer accu_nb_ghost=0;
  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    m_nb_owned_evi_pn[inei] = m_mmvs->sharedItems(inei).size();
    m_nb_ghost_evi_pn[inei] = m_mmvs->ghostItems(inei).size();

    m_indexes_owned_evi_pn[inei] = accu_nb_owned;
    m_indexes_ghost_evi_pn[inei] = accu_nb_ghost;
    
    accu_nb_owned += m_nb_owned_evi_pn[inei];
    accu_nb_ghost += m_nb_ghost_evi_pn[inei];
  }
  m_buf_owned_evi.resize(accu_nb_owned);
  m_buf_ghost_evi.resize(accu_nb_ghost);

  // On construit des multi-vues sur des zones allouées en mémoire managées
  MultiArray2View<EnvVarIndex> owned_evi_pn(m_buf_owned_evi.view(),
      m_indexes_owned_evi_pn.constView(), m_nb_owned_evi_pn.constView());

  MultiArray2View<EnvVarIndex> ghost_evi_pn(m_buf_ghost_evi.view(),
      m_indexes_ghost_evi_pn.constView(), m_nb_ghost_evi_pn.constView());

  // On convertit les MatVarIndex(es) en EnvVarIndex(es)
  auto mvi2evi = [](ConstArrayView<MatVarIndex> lmvis, ArrayView<EnvVarIndex> levis)
  {
    for(Integer i=0 ; i<lmvis.size() ; ++i) {
      const MatVarIndex& mvi = lmvis[i];
      levis[i] = EnvVarIndex(mvi.arrayIndex(), mvi.valueIndex());
    }
  };

  for(Integer inei=0 ; inei<m_nb_nei ; ++inei) {
    mvi2evi(m_mmvs->sharedItems(inei), owned_evi_pn[inei]);
    mvi2evi(m_mmvs->ghostItems(inei) , ghost_evi_pn[inei]);
  }

  // Puisque les adresses ont potentiellement changé : 
  // "Conseils" mémoire
  if (m_acc_mem_adv) {
    m_acc_mem_adv->setReadMostly(m_buf_owned_evi       .view());
    m_acc_mem_adv->setReadMostly(m_indexes_owned_evi_pn.view());
    m_acc_mem_adv->setReadMostly(m_nb_owned_evi_pn     .view());
    m_acc_mem_adv->setReadMostly(m_buf_ghost_evi       .view());
    m_acc_mem_adv->setReadMostly(m_indexes_ghost_evi_pn.view());
    m_acc_mem_adv->setReadMostly(m_nb_ghost_evi_pn     .view());
  }
}

