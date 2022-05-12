#ifndef ACC_ENV_DEFAULT_SERVICE_H
#define ACC_ENV_DEFAULT_SERVICE_H

#include "accenv/IAccEnv.h"
#include "arcane/materials/MeshMaterialVariableRef.h"

// Ces fichiers doivent être inclus avant AccEnvDefault_axl.h
#include "accenv/AccEnvDefaultOptions.h"
#include "msgpass/VarSyncMngOptions.h"
#include "accenv/AccEnvDefault_axl.h"

using namespace Arcane;

class AccEnvDefaultService : public ArcaneAccEnvDefaultObject
{
 public:
  AccEnvDefaultService(const ServiceBuildInfo & sbi);

  virtual ~AccEnvDefaultService();

 public:
  void initAcc() override;

  ax::Runner& runner() override { return m_runner; }
  ax::RunQueue newQueue() override { return makeQueue(m_runner); }
  Ref<ax::RunQueue> refQueueAsync(eQueuePriority qp=QP_default) override;

  AccMemAdviser* accMemAdv() override { return m_acc_mem_adv; }

  UnstructuredMeshConnectivityView& connectivityView() override { return m_connectivity_view; }
  const UnstructuredMeshConnectivityView& connectivityView() const override { return m_connectivity_view; }
  
  void initMesh(IMesh* mesh) override;

  Span<const Int16> nodeIndexInCells() const override { return m_node_index_in_cells.constSpan(); }

  Integer maxNodeCell() const override { return 8; }

  void initMultiEnv(IMeshMaterialMng* mesh_material_mng) override;
  MultiAsyncRunQueue* multiEnvQueue() override { return m_menv_queue; }

  void computeMultiEnvGlobalCellId(IMeshMaterialMng* mesh_material_mng) override;
  void checkMultiEnvGlobalCellId(IMeshMaterialMng* mesh_material_mng) override;
  void updateMultiEnv(IMeshMaterialMng* mesh_material_mng) override;

  MultiEnvCellStorage* multiEnvCellStorage() override { return m_menv_cell; }

  VarSyncMng* vsyncMng() override { return m_vsync_mng; }

 protected:

  void _computeNodeIndexInCells();

 protected:
  ax::Runner m_runner;
  AccMemAdviser* m_acc_mem_adv=nullptr;
  UnstructuredMeshConnectivityView m_connectivity_view;

  UniqueArray<Int16> m_node_index_in_cells;

  //! Description/accès aux mailles multi-env
  MultiEnvCellStorage* m_menv_cell=nullptr;

  // Les queues asynchrones d'exéution
  MultiAsyncRunQueue* m_menv_queue=nullptr; //!< les queues pour traiter les environnements de façon asynchrone

  //! Pour "synchroniser" les items fantômes
  VarSyncMng* m_vsync_mng=nullptr;
};

#endif

