#ifndef ACC_ENV_MULTI_ENV_MNG_H
#define ACC_ENV_MULTI_ENV_MNG_H

#include <arcane/materials/MeshMaterialVariableRef.h>

#include "accenv/AcceleratorUtils.h"
#include "accenv/MultiEnvUtils.h"
#include "msgpass/VarSyncMng.h"

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/* Manager pour la descritpion du multi-environnement sur le maillage        */
/*---------------------------------------------------------------------------*/
class MultiEnvMng {
 public:
  MultiEnvMng(IMeshMaterialMng* mm, ax::Runner& runner, VarSyncMng* vsync_mng, AccMemAdviser* acc_mem_adv);

  virtual ~MultiEnvMng();

  //! Calcul des cell_id globaux : permet d'associer à chaque maille impure (mixte)
  // l'identifiant de la maille globale
  void computeMultiEnvGlobalCellId();

  void checkMultiEnvGlobalCellId();

  //! Préparer les données multi-environnement pour l'accélérateur
  // A appeler quand la carte des environnements change
  void updateMultiEnv(VarSyncMng* vsync_mng);

  //! Remplit les tableaux IsActiveCell et IsActiveNode à partir des groupes d'items actifs
  void setActiveItemsFromGroups(CellGroup active_cells, NodeGroup active_nodes);
  
  //! Vue sur le stockage pour utilisation en lecture sur GPU
  MultiEnvCellViewIn viewIn(ax::RunCommand& command) {
    return MultiEnvCellViewIn(command, m_max_nb_env, m_nb_env, m_l_env_arrays_idx, m_l_env_values_idx, m_env_id);
  }
  
  //! les queues pour traiter les environnements de façon asynchrone
  MultiAsyncRunQueue* multiEnvQueue() { return m_menv_queue; }

  //! Nb d'environnement par maille
  VariableCellInteger nbEnv() {
    return m_nb_env;
  }

  //! Si maille pure: identifiant de l'environnement de la maille, si maille mixte: opposé du nombre d'environnements+1, si maille vide: -1
  VariableCellInteger envId() {
    return m_env_id;
  }

  //! globalCell()[envcell] : local id de la maille globale
  MaterialVariableCellInteger globalCell() {
    return m_global_cell;
  }

  //! Indique si la maille est active (= non 100% vide)
  VariableCellBool isActiveCell() {
    return m_is_active_cell;
  }

  //! Indique si le noeud est actif (= dont au moins une de ces mailles est active)
  VariableNodeBool isActiveNode() {
    return m_is_active_node;
  }

 public:

  //! Remplissage
  void buildStorage(ax::Runner& runner, Materials::MaterialVariableCellInteger& v_global_cell);

  //! Verification
  void checkStorage(Materials::MaterialVariableCellInteger& v_global_cell);

 protected:

  IMeshMaterialMng* m_mesh_material_mng=nullptr;
  ax::Runner& m_runner;
  AccMemAdviser* m_acc_mem_adv = nullptr;

  Integer m_max_nb_env;

  /* Les variables Arcane qui décrivent la répartition des environnements sur le maillage */
  VariableCellInteger m_nb_env;  //! Nb d'env par maille
  UniqueArray<Int16> m_l_env_arrays_idx; //! liste des indexes des env par maille
  VariableCellArrayInteger m_l_env_values_idx;
  VariableCellInteger m_env_id;  //! Si maille pure: identifiant de l'environnement de la maille, si maille mixte: opposé du nombre d'environnements+1, si maille vide: -1
  MaterialVariableCellInteger m_global_cell;  //! m_global_cell[envcell] : local id de la maille globale
  VariableCellBool m_is_active_cell;  //! Indique si la maille est active (= non 100% vide)
  VariableNodeBool m_is_active_node;  //! Indique si le noeud est actif (= dont au moins une de ces mailles est active)

  // Les queues asynchrones d'exéution
  MultiAsyncRunQueue* m_menv_queue=nullptr; //!< les queues pour traiter les environnements de façon asynchrone
};

#endif

