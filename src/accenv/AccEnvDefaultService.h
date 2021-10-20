#ifndef ACC_ENV_DEFAULT_SERVICE_H
#define ACC_ENV_DEFAULT_SERVICE_H

#include "accenv/IAccEnv.h"
#include "accenv/AccEnvDefault_axl.h"

using namespace Arcane;

class AccEnvDefaultService : public ArcaneAccEnvDefaultObject
{
 public:
  AccEnvDefaultService(const ServiceBuildInfo & sbi)
    : ArcaneAccEnvDefaultObject(sbi) {}

  virtual ~AccEnvDefaultService() {}

 public:
  void initAcc() override;

  ax::Runner& runner() override { return m_runner; }
  ax::RunQueue newQueue() override { return makeQueue(m_runner); }

  UnstructuredMeshConnectivityView& connectivityView() override { return m_connectivity_view; }
  const UnstructuredMeshConnectivityView& connectivityView() const override { return m_connectivity_view; }
  void initMesh(IMesh* mesh) override;

 protected:
  ax::Runner m_runner;
  UnstructuredMeshConnectivityView m_connectivity_view;
};

#endif

