#include "AccEnvModule.h"

#include <arcane/materials/IMeshMaterialMng.h>

#include "accenv/SingletonIAccEnv.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Acc_env {

/*---------------------------------------------------------------------------*/
/* Constructeur */
/*---------------------------------------------------------------------------*/
AccEnvModule::
  AccEnvModule(const ModuleBuildInfo& mbi)
: ArcaneAccEnvObject(mbi)
{}

/*---------------------------------------------------------------------------*/
/* Destructeur */
/*---------------------------------------------------------------------------*/
AccEnvModule::~AccEnvModule() {
}

/*---------------------------------------------------------------------------*/
/* Initialise les environnements pour les accélérateurs */
/*---------------------------------------------------------------------------*/

void AccEnvModule::
accBuild() 
{
  prof_acc_begin("AccEnvModule::accBuild");

  m_acc_env = SingletonIAccEnv::accEnv(subDomain());
  m_acc_env->initAcc();

  prof_acc_end("AccEnvModule::accBuild");
}

/*---------------------------------------------------------------------------*/
/* Initialise des informations liées au maillage pour les accélérateurs */
/*---------------------------------------------------------------------------*/

void AccEnvModule::
initMesh() 
{
  prof_acc_begin("AccEnvModule::initMesh");

  m_acc_env->initMesh(defaultMesh());

  prof_acc_end("AccEnvModule::initMesh");
}

/*---------------------------------------------------------------------------*/
/* Initialise des informations liées au multi-environnement pour les accélérateurs */
/*---------------------------------------------------------------------------*/

void AccEnvModule::
initMultiEnv() 
{
  prof_acc_begin("AccEnvModule::initMultiEnv");

  IMeshMaterialMng* mesh_material_mng = IMeshMaterialMng::getReference(defaultMesh());
  m_acc_env->createMultiEnvMng(mesh_material_mng);

  prof_acc_end("AccEnvModule::initMultiEnv");
}

/*---------------------------------------------------------------------------*/
/* Amorce l'instrumentation pour le profiling */
/*---------------------------------------------------------------------------*/

void AccEnvModule::
startInstrument() 
{
  prof_acc_start_capture();
}

/*---------------------------------------------------------------------------*/
/* Termine l'instrumentation pour le profiling */
/*---------------------------------------------------------------------------*/

void AccEnvModule::
stopInstrument() 
{
  prof_acc_stop_capture();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_ACCENV(AccEnvModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}  // namespace Acc_env
