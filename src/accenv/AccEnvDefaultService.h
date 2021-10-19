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

 protected:
  ax::Runner m_runner;
};

#endif

