// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef TABULEEEOSSERVICE_H
#define TABULEEEOSSERVICE_H

#include <map>
#include <iostream>
#include <fstream>
#include <string> 
#include <filesystem>

#include "eos/IEquationOfState.h"
#include "eos/tabulee/TabuleeEOS_axl.h"
#include "arcane/materials/IMeshMaterialMng.h"
#include "arcane/materials/IMeshMaterial.h"
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/IMeshBlock.h"
#include "arcane/materials/MeshMaterialModifier.h"
#include "arcane/materials/MeshMaterialVariableRef.h"
#include "arcane/materials/MeshEnvironmentVariableRef.h"
#include "arcane/materials/MaterialVariableBuildInfo.h"
#include "arcane/materials/MeshBlockBuildInfo.h"
#include "arcane/materials/MeshEnvironmentBuildInfo.h"
#include "arcane/materials/MeshMaterialVariableDependInfo.h"
#include "arcane/materials/CellToAllEnvCellConverter.h"
#include "arcane/materials/MatCellVector.h"
#include "arcane/materials/EnvCellVector.h"
#include "arcane/materials/MatConcurrency.h"
#include "arcane/materials/MeshMaterialIndirectModifier.h"
#include "arcane/materials/MeshMaterialVariableSynchronizerList.h"
#include "arcane/materials/ComponentSimd.h"
using namespace Arcane;
using namespace Arcane::Materials;


/**
 * Représente le modèle d'équation d'état tabulée
 */
class TabuleeEOSService 
: public ArcaneTabuleeEOSObject
{
public:
  /** Constructeur de la classe */
  TabuleeEOSService(const ServiceBuildInfo & sbi)
    : ArcaneTabuleeEOSObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~TabuleeEOSService() {};

public:
  /** 
   *  Initialise l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et l'énergie interne. 
   */
  virtual void initEOS(IMeshEnvironment* env);
   /** 
   *  Initialise de certains elements de l'équation d'état en reprise. 
   */
  virtual void ReinitEOS(IMeshEnvironment* env);
  /** 
   *  Applique l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et la pression. 
   */
  virtual void applyEOS(IMeshEnvironment* env);
   /** 
   *  Applique l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et la pression pour une cellule
   */
  virtual void applyOneCellEOS(IMeshEnvironment* env, EnvCell ev);
   /** 
   *  Applique un endommagement dans la maille (pression nulle)
   *  si la maille est endessous de la tension limite. 
   */
  virtual void Endommagement(IMeshEnvironment* env);
  /** 
   *  Renvoie la constante adiabatic de l'environnement. 
   */
  virtual Real getAdiabaticCst(IMeshEnvironment* env);
  /** 
   *  Renvoie la constante tension limit de l'environnement. 
   */
  virtual Real getTensionLimitCst(IMeshEnvironment* env);
  /** 
   *  Renvoie la chaleur spécifique de l'environnement. 
   */
  virtual Real getSpecificHeatCst(IMeshEnvironment* env);
  /** 
   *  Renvoie le rapport seuil de densité  
   */
  virtual Real getdensityDamageThresold(IMeshEnvironment* env);
private:
  // taille des tableaux = nombre de maille de l'environement
  Integer nbmail;
  double *m_dtime,*m_rho,*m_ene,*m_Pres,*m_Temp;
  double *m_Frac1,*m_Frac2,*m_Frac3,*m_Frac4,*m_Frac5,*m_Frac6;
  double *m_pente_dpde,*m_Cv, *m_cs2,*m_conv;
  
  std::map<std::string, std::vector<std::vector<double>>> data;
  double emin, emax, dmin, dmax, tmin, tmax;
  void initdata();
  
  bool calculPetT(double d0, double e0, double& t0, double& p0, std::string interp);
  bool calculPetE(double t0, double d0, double& p0, double& e0, std::string interp);
  
  std::vector<double> splitAndConvertToDoubleForT(const std::string& line, char delimiter) {
    std::vector<std::string> tokens = splitString(line, delimiter);
    std::vector<double> result;
	
	std::vector<std::string> tab = splitString(  tokens[1]);
	for (int i = 2 ; i< tab.size() ; ++i) 
	  if (i%2 == 0) {
            result.push_back(std::stod(tab[i]));
	}
    return result;
  };
  
  std::vector<double> splitAndConvertToDouble(const std::string& line, char delimiter) {
    std::vector<std::string> tokens = splitString(line, delimiter);
    std::vector<double> result;
	
	for (int i = 4 ; i< tokens.size() ; ++i) 
	  if (i%2 == 0) {
            result.push_back(std::stod(tokens[i]));
	}
    return result;
  };
  
  std::vector<std::string> splitString(const std::string& str, char delimiter = ' ') {
    std::vector<std::string> tokens;
    std::istringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
  };
};

#endif
