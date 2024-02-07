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
/* Initialise les environnements avec l'API pré-Accélérateur */
/*---------------------------------------------------------------------------*/

void AccEnvModule::
paccBuild() 
{
  prof_acc_begin("AccEnvModule::paccBuild");

  accBuild();
  m_acc_env->initPAcc();

  prof_acc_end("AccEnvModule::paccBuild");
}

/*---------------------------------------------------------------------------*/
/* Initialise des informations liées au maillage pour les accélérateurs */
/*---------------------------------------------------------------------------*/

void AccEnvModule::
initMesh() 
{
  prof_acc_begin("AccEnvModule::initMesh");

  m_acc_env->initMesh(defaultMesh());

  // active la fonctionnalité des allenvell pour les runcommand
  // également possible de le faire via la variable d'env ARCANE_ALLENVCELL_FOR_RUNCOMMAND
  IMeshMaterialMng* mesh_material_mng = IMeshMaterialMng::getReference(defaultMesh());
  mesh_material_mng->enableCellToAllEnvCellForRunCommand(true);

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

  // active la fonctionnalité des allenvcell pour les runcommand
  // également possible de le faire via la variable d'environnement ARCANE_ALLENVCELL_FOR_RUNCOMMAND
  mesh_material_mng->enableCellToAllEnvCellForRunCommand(true);

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
