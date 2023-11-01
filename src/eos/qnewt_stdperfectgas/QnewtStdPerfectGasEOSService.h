// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef QNEWT_STDPERFECTGASEOSSERVICE_H
#define QNEWT_STDPERFECTGASEOSSERVICE_H

#include "arcane/materials/CellToAllEnvCellConverter.h"
#include "arcane/materials/ComponentSimd.h"
#include "arcane/materials/EnvCellVector.h"
#include "arcane/materials/IMeshBlock.h"
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/IMeshMaterial.h"
#include "arcane/materials/IMeshMaterialMng.h"
#include "arcane/materials/MatCellVector.h"
#include "arcane/materials/MatConcurrency.h"
#include "arcane/materials/MaterialVariableBuildInfo.h"
#include "arcane/materials/MeshBlockBuildInfo.h"
#include "arcane/materials/MeshEnvironmentBuildInfo.h"
#include "arcane/materials/MeshEnvironmentVariableRef.h"
#include "arcane/materials/MeshMaterialIndirectModifier.h"
#include "arcane/materials/MeshMaterialModifier.h"
#include "arcane/materials/MeshMaterialVariableDependInfo.h"
#include "arcane/materials/MeshMaterialVariableRef.h"
#include "arcane/materials/MeshMaterialVariableSynchronizerList.h"
#include "eos/IEquationOfState.h"
#include "qnewt_stdperfectgas/QnewtStdPerfectGasEOS_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;


class IAccEnv;

namespace Qnewt_stdperfectgas {

/**
 * Représente le modèle d'équation d'état <em>Gaz Parfait</em>
 */
class QnewtStdPerfectGasEOSService 
: public ArcaneQnewtStdPerfectGasEOSObject
{
public:
  /** Constructeur de la classe */
  QnewtStdPerfectGasEOSService(const ServiceBuildInfo & sbi);
  
  /** Destructeur de la classe */
  virtual ~QnewtStdPerfectGasEOSService() {};

public:
  /** 
   *  Initialise l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et l'énergie interne. 
   */
  virtual void initEOS(IMeshEnvironment* env) override;
  /** 
   *  Applique l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et la pression. 
   */
  virtual void applyEOS(IMeshEnvironment* env) override;
   /** 
   *  Applique l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et la pression pour une cellule
   */
  virtual void applyOneCellEOS(IMeshEnvironment* env, EnvCell ev) override;
  /** 
   *  Renvoie la constante adiabatic de l'environnement. 
   */
  Real getAdiabaticCst(IMeshEnvironment* env) override;
  /** 
   *  Renvoie la constante tension limit de l'environnement. 
   */
  Real getTensionLimitCst(IMeshEnvironment* env) override;

private:
  // Pour l'utilisation des accélérateurs
  IAccEnv* m_acc_env=nullptr;
};

} // Qnewt_stdperfectgas

#endif
