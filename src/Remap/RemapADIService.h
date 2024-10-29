// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
#ifndef IREMAPADISERVICE_H
#define IREMAPADISERVICE_H

#include "TypesMahyco.h"

#include "Remap/IRemap.h"
#include "Remap/RemapADI_axl.h"

using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Représente le service de Remap version ADI
 */
class RemapADIService 
: public ArcaneRemapADIObject
{
public:
  /** Constructeur de la classe */
  RemapADIService(const ServiceBuildInfo & sbi)
    : ArcaneRemapADIObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~RemapADIService() {};
  
  struct interval {
    double inf, sup;
  };
  ICartesianMesh* m_cartesian_mesh;
  Materials::IMeshMaterialMng* mm;

public:
   /**
   * main du remap
   **/
  virtual void appliRemap(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env);
   /**
   * resize les variables du remap
   **/
  virtual  void resizeRemapVariables(Integer nb_vars_to_project, Integer nb_env);
  /**
   * synchronisation des valeurs aux cellules 
   **/
  virtual void synchronizeUremap();
  
  
  virtual Integer getOrdreProjection();
  virtual bool hasProjectionPenteBorne();
  virtual bool hasProjectionSimplePente();    
  virtual bool hasConservationEnergieTotale();   
  virtual String isEuler();
  
    /**
   * fonction final de la projection
   **/
  virtual void remapVariables(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env) ;
   
private:
    
  /**
   * Calcule les résultantes aux noeuds d'une maille hexaédrique.
   * La méthode utilisée est celle du découpage en quatre triangles.
   * Méthode appelée par le point d'entrée \c computeGeometricValues()
   */
  inline void computeCQs(Real3 node_coord[8],Real3 face_coord[6],const Cell& cell);
  /** 
   * calcul du nouveau volume dit Euler si différent du maillage initial    
   **/ 
   void computeVolumeEuler(Integer idir);
  /**
   * calcul des gradients aux faces
   **/
   void computeGradPhiFace(Integer idir, Integer nb_vars_to_project, Integer nb_env);
  /**
   * calcul des gradients aux faces ou flux aux faces 
   **/
   void computeGradPhiCell(Integer idir, Integer nb_vars_to_project, Integer nb_env);
  /**
   * calcul des flux aux faces des cellules 
   **/
   void computeUpwindFaceQuantitiesForProjection(Integer idir, Integer nb_vars_to_project, Integer nb_env);
  /**
   * calcul des valeurs aux cellules 
   **/
   void computeUremap(Integer idir, Integer nb_vars_to_project, Integer nb_env, Integer withDualProjection);

  /**
   * fonction pour la phase de projection duales
   **/
   void computeDualUremap(Integer idir, Integer nb_env);
  /**
   * synchronisation des valeurs aux noeuds 
   **/
   void synchronizeDualUremap();
  /**
   * calcul le gradient aux mailles et reconstruction limitée d'une quantité Phi 
   **/
   void computeAndLimitGradPhi(Integer projectionLimiterId, Face frontFace, Face backFace, 
                                      Cell cell, Cell frontcell, Cell backcell, int nb_vars) ;
  /**
   * calcul le gradient aux mailles et reconstruction limitée d'une quantité dual Phi Dual
   **/
   void computeAndLimitGradPhiDual(Integer projectionLimiterId, 
                                          Node inode, Node frontnode, Node backnode,
                           Real3 grad_front, Real3 grad_back, Real h0, Real hplus, Real hmoins);
  /**
   * calcul le gradient aux mailles d'une quantité dual Phi Dual
   **/
   void computeDualGradPhi(Node inode, Node frontfrontnode, Node frontnode, 
                                   Node backnode, Node backbacknode, Integer idir);

    /**
   * Fonction limiteur de gradient (maille regulier)
   **/
    Real fluxLimiter(Integer projectionLimiterId, double r);
  /**
   * Fonction limiteur de gradient (maille irregulier)
   **/
    Real fluxLimiterG(Integer projectionLimiterId, double gradplus,
                           double gradmoins, double y0, double yplus,
                           double ymoins, double h0, double hplus,
                           double hmoins);
  /**
   * Calcul des seuils de monotonie des reconstructions simple pente 
   **/
    Real computeY0(Integer projectionLimiterId, double y0, double yplus,
                        double ymoins, double h0, double hplus, double hmoins,
                        Integer type);
   /**
   * Calcul des abscisses des pointes d'appui de la reconstuctions pente-borne
   **/
    Real computexgxd(double y0, double yplus, double ymoins, double h0,
                          double y0plus, double y0moins, Integer type);
   /**
   * Calcul des oordonnes des pointes d'appui de la reconstuctions pente-borne
   **/
  
    Real computeygyd(double y0, double yplus, double ymoins, double h0,
                          double y0plus, double y0moins, double grady,
                          Integer type);
  
  /**
   * Calcul des flux pente-borne (integration de  la reconstruction en trois morceaux) 
   * cas des mailles mixtes ou le flux de volume determine les flux de masse et d'energie
   **/
   void computeFluxPP(Cell cell, Cell frontcell, Cell backcell, 
                                     Real face_normal_velocity, 
                                     Real deltat_n, Integer type, Real flux_threshold, 
                                     Integer projectionPenteBorneComplet, 
                                     Real dual_normal_velocity,
                                     Integer calcul_flux_dual,
                                     RealArrayView Flux, RealArrayView Flux_dual,
                                     int nbmat, int nb_vars);
  /**
   * Calcul des flux pente-borne (integration de  la reconstruction en trois morceaux) 
   * cas des mailles pures ou le flux de masse determine le flux d'energie
   **/
   void computeFluxPPPure(Cell cell, Cell frontcell, Cell backcell, 
                                     Real face_normal_velocity, 
                                     Real deltat_n, Integer type, Real flux_threshold, 
                                     Integer projectionPenteBorneComplet, 
                                     Real dual_normal_velocity,
                                     Integer calcul_flux_dual,
                                     RealArrayView Flux, RealArrayView Flux_dual,
                                     int nbmat, int nb_vars);

  /**
   * calcul des flux d'une variable Phi final 
   **/
   Real computeRemapFlux(Integer projectionOrder, Integer projectionAvecPlateauPente,
        Real face_normal_velocity, Real3 face_normal,
        Real face_length, Real phi_face,
        Real3 outer_face_normal, Real3 exy, Real deltat_n);
  
  Real INTY(double X, double x0, double y0, double x1, double y1);
  
  /**
   * Fonctions pour l'ordre 3
   **/
  Real2 define_interval(double a, double b);
  Real2 intersection(Real2 I1, Real2 I2);
  double evaluate_grad(double hm, double h0, double hp, double ym, double y0,
                       double yp);
  double evaluate_ystar(double hmm, double hm, double hp, double hpp,
                        double ymm, double ym, double yp, double ypp,
                        double gradm, double gradp);
  double evaluate_fm(double x, double dx, double up, double du, double u6);
  
  double evaluate_fp(double x, double dx, double um, double du, double u6);
  double ComputeFluxOrdre3(double ymmm, double ymm, double ym, double yp,
                           double ypp, double yppp, double hmmm, double hmm,
                           double hm, double hp, double hpp, double hppp,
                           double v_dt);
  
  
  Real m_arithmetic_thresold = 1.e-300;
  
};

#endif
