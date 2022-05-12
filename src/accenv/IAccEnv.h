#ifndef IACCENV_H
#define IACCENV_H

#include <arcane/ItemTypes.h>
#include "accenv/AcceleratorUtils.h"
#include "accenv/MultiEnvUtils.h"
#include "msgpass/VarSyncMng.h"
#include "arcane/UnstructuredMeshConnectivity.h"
#include "arcane/materials/IMeshMaterialMng.h"
#include "arcane/IMesh.h"

using namespace Arcane;

class IAccEnv {
 public:
  IAccEnv() {}
  virtual ~IAccEnv() {}

 public:
  virtual void initAcc() = 0;

  virtual ax::Runner& runner() = 0;
  virtual ax::RunQueue newQueue() = 0;
  virtual Ref<ax::RunQueue> refQueueAsync(eQueuePriority qp=QP_default) = 0;

  virtual AccMemAdviser* accMemAdv() = 0;

  virtual UnstructuredMeshConnectivityView& connectivityView() = 0;
  virtual const UnstructuredMeshConnectivityView& connectivityView() const = 0;

  virtual void initMesh(IMesh* mesh) = 0;

  virtual Span<const Int16> nodeIndexInCells() const = 0;

  virtual Integer maxNodeCell() const = 0;

  virtual void initMultiEnv(IMeshMaterialMng* mesh_material_mng) = 0;
  virtual MultiAsyncRunQueue* multiEnvQueue() = 0;

  virtual void computeMultiEnvGlobalCellId(IMeshMaterialMng* mesh_material_mng) = 0;
  virtual void checkMultiEnvGlobalCellId(IMeshMaterialMng* mesh_material_mng) = 0;
  virtual void updateMultiEnv(IMeshMaterialMng* mesh_material_mng) = 0;

  virtual MultiEnvCellStorage* multiEnvCellStorage() = 0;

  virtual VarSyncMng* vsyncMng() = 0;
};

#endif
