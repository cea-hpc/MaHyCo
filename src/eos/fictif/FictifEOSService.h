﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef FICTIFEOSSERVICE_H
#define FICTIFEOSSERVICE_H

#include "eos/IEquationOfState.h"
#include "eos/fictif/FictifEOS_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Représente le modèle d'équation d'état <em>Gaz Parfait</em>
 */
class FictifEOSService 
: public ArcaneFictifEOSObject
{
public:
  /** Constructeur de la classe */
  FictifEOSService(const ServiceBuildInfo & sbi)
    : ArcaneFictifEOSObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~FictifEOSService() {};

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
};

#endif
