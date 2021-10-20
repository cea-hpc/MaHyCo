#ifndef IACCENV_H
#define IACCENV_H

#include <arcane/ItemTypes.h>
#include "accenv/AcceleratorUtils.h"
#include "arcane/UnstructuredMeshConnectivity.h"
#include "cartesian/interface/ICartesianMesh.h"

using namespace Arcane;

class IAccEnv {
 public:
  IAccEnv() {}
  virtual ~IAccEnv() {}

 public:
  virtual void initAcc() = 0;

  virtual ax::Runner& runner() = 0;
  virtual ax::RunQueue newQueue() = 0;

  virtual UnstructuredMeshConnectivityView& connectivityView() = 0;
  virtual const UnstructuredMeshConnectivityView& connectivityView() const = 0;

  virtual void initMesh(ICartesianMesh* cartesian_mesh) = 0;

  virtual Span<const Int16> nodeIndexInCells() const = 0;
  virtual Span<const Int16> nodeIndexInFaces() const = 0;

  virtual Integer maxNodeCell() const = 0;
  virtual Integer maxNodeFace() const = 0;
};

#endif
