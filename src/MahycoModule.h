// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef MAHYCOMODULE_H
#define MAHYCOMODULE_H

#include "TypesMahyco.h"

#include <arcane/geometry/IGeometryMng.h>
// Ajout au PIF
#include "arcane/utils/List.h"
#include "arcane/utils/OStringStream.h"
#include "arcane/utils/ValueChecker.h"
#include "arcane/utils/SimdOperation.h"

#include "arcane/IUnitTest.h"
#include "arcane/ITimeLoopMng.h"
#include "arcane/ITimeLoopService.h"
#include "arcane/ITimeLoop.h"
#include "arcane/TimeLoopEntryPointInfo.h"
#include "arcane/IMesh.h"
#include "arcane/IItemFamily.h"
#include "arcane/ItemPrinter.h"
#include "arcane/IParallelMng.h"
#include "arcane/IMeshModifier.h"
#include "arcane/IMeshUtilities.h"
#include "arcane/IMeshPartitioner.h"
#include "arcane/VariableDependInfo.h"

#include "arcane/Concurrency.h"
#include "arcane/VariableView.h"

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

#include "arcane/cea/ICartesianMesh.h"
#include "arcane/cea/CellDirectionMng.h"
#include "arcane/cea/FaceDirectionMng.h"
#include "arcane/cea/NodeDirectionMng.h"
#include "arcane/cea/CartesianConnectivity.h"
// fin ajout au PIF

#include "Mahyco_axl.h"

using namespace Arcane;

/**
 * Représente un module d'hydrodynamique lagrangienne très simplifié :
 *   - le seul type de maille supporté est l'hexaèdre,
 *   - pas de pseudo viscosité supportée,
 *   - le seul type de calcul de longueur caractéristique supporté est celui utilisant les médianes,
 *   - le seul type de condition aux limites supporté est d'imposer une composante de la vitesse sur une surface,
 *   - la masse nodale est supposée constante et n'est pas recalculée à chaque itération,
 *   - aucun test de cohérence des valeurs (pression positive, volume positif, ...)  n'est effectué.
 *  
 * La liste des opérations effectuées par le module est la suivante :
 *   - calcul des forces de pression,
 *   - calcul de l'impulsion,
 *   - prise en compte des conditions aux limites,
 *   - déplacement des noeuds,
 *   - calcul des nouvelles valeurs géométriques : volume des mailles, longueur caractéristique des mailles,
 *     resultantes aux sommets de chaque maille,
 *   - calcul de la densité,
 *   - calcul de la pression et de l'énergie par l'équation d'état. Ce calcul est effectué par un service
 *     ARCANE. Deux implémentations sont disponibles pour le service : gaz parfait, et "stiffened" gaz.
 *   - calcul du nouveau pas de temps.
 * 
 */
class MahycoModule
: public ArcaneMahycoObject
{
 public:
  /** Constructeur de la classe */
  MahycoModule(const ModuleBuildInfo& mbi)
    : ArcaneMahycoObject(mbi) {}
  /** Destructeur de la classe */
  ~MahycoModule() {}
  
  struct interval {
    double inf, sup;
  };
 public:
  /** 
   *  Initialise le module. 
   *  L'initialisation comporte deux parties distinctes:
   *  - la première partie où il faut indiquer la taille des variables
   *    tableaux. Dans notre cas, il s'agit de \c m_cell_cqs et
   *    \c m_viscosity_force, qui sont toutes deux des variables
   *    aux mailles possédant une valeur pour chaque noeud de chaque
   *    maille. Comme on ne supporte que les héxaèdres, il y a 8 valeurs 
   *    par maille,
   *  - la deuxième partie qui consiste à initialiser les variables avec
   *    leur valeur de départ. Pour les variables \c Pressure, \c Density et
   *    \c AdiabaticCst, c'est ARCANE qui les initialisent directement
   *    à partir du jeu de donnée. La variable \c NodeCoord est aussi
   *    initialisée par l'architecture lors de la lecture du maillage. Les
   *    autres variables sont calculées comme suit :
   *    - le pas de temps initial est donné par le jeu de donnée,
   *    - les valeurs géométriques (longueur caractéristique, volume et
   *      résultantes aux sommets) sont calculées à partir des coordonnées des
   *      noeuds,
   *    - la masse des mailles est calculée à partir de sa densité et de
   *      son volume,
   *    - la masse des mailles et la masse nodale. La masse d'une maille
   *      est calculée à partir de sa densité et de son volume,
   *    - la masse nodale est calculée en ajoutant les contributions de
   *      chaque maille connecté à un noeud donné. Chaque maille
   *      contribue pour 1/8ème de sa masse à la masse nodale de chacun de ses
   *      sommets,
   *    - l'énergie interne et la vitesse du son sont calculées en fonction 
   *      de l' équation d'état.
   *    - appelle hydroStartInitEnvAndMat() pour la création des env et des mat
   *      en fonction du maillage
   */
  virtual void hydroStartInit();
  virtual void hydroStartInitEnvAndMat();
  virtual void PrepareFaceGroup();
  
  /** 
   *  Initialise le module en cas de reprise. 
   *  Afin d'éviter de sauvegarder le volume des mailles, cette méthode
   *  recalcule le volume en fonction des coordonnées.
   */
  virtual void hydroContinueInit();
  
  /** 
   *  Sauvegarde des variables à l'instant n 
   */
  virtual void saveValuesAtN();

   /** 
   * Calcule la masse des mailles
   */
  virtual void computeCellMass();
   /** 
   * Calcule la masse des noeuds
   */
  virtual void computeNodeMass();

  /** 
   * Calcule la pseudo viscosité au maille
   * au temps courant \f$t^{n}\f$. Pour chaque maille, 
   * il s'y a une partie linéaire et quadratique vis-vis du saut de vitesse .
   * Calcule la pseudo viscosité au maille au temps courant \q$t^{n}\f$.
   */
  virtual void computeArtificialViscosity();
  
		
  /**
   * Calcule la force (\c m_force) qui s'applique aux noeuds en
   * ajoutant l'éventuelle contribution de la pseudo-viscosité. Calcule 
   * ensuite la nouvelle vitesse (\c m_velocity) aux noeuds a n+1/2
   * en fonction de celle à n-1/2
   */
  virtual void updateVelocity();
  /**
   * Calcul de la vitesse de n a n-1/2
   * Pour CSTS
   */
  virtual void updateVelocityBackward();
  /**
   * Calcul de la vitesse de n+1/2 a n+1
   * Pour CSTS
   */
  virtual void updateVelocityForward();
  /* 
   * pour les cas d'avection pure
   */
  virtual void updateVelocityWithoutLagrange();

		
  /**
   * Applique les conditions aux limites.
   * Les conditions aux limites dépendent des options du
   * jeu de données. Dans cette implémentation, une condition aux limites
   * possède les propriétés suivantes :
   * - un type: trois types sont supportés: contraindre la composante
   * \f$x\f$ du vecteur vitesse, contraindre la composante \f$y\f$ du vecteur
   * vitesse ou contraindre la composante \f$z\f$ du vecteur vitesse,
   * - une valeur: il s'agit d'un réel indiquant la valeur de la
   * contrainte,
   * - une surface: il s'agit de la surface sur laquelle s'applique la
   * contrainte.
   * 
   * Appliquer les conditions aux limites consiste donc à fixer une
   * composante d'un vecteur vitesse pour chaque noeud de chaque face de
   * chaque surface sur laquelle on impose une condition aux limites.
   */		
  virtual void applyBoundaryCondition();
		
  /**
   * Modifie les coordonnées (\c m_node_coord)
   * des noeuds d'après la valeur du vecteur vitesse et du pas de temps.
   */
  virtual void updatePosition();
		
  /**
   * Ce point d'entrée regroupe l'ensemble des calculs géométriques
   * utiles pour le schéma. Dans notre cas, il s'agit pour chaque maille :
   * - de calculer sa longueur caractéristique,
   * - de calculer les résultantes à ses sommets,
   * - de calculer son volume.
   
   * Pour optimiser le calcul (utilisation du cache), à chaque itération 
   * sur une maille, sont stockées localement les coordonnées de ses noeuds 
   * et celles du centre de ses faces.
   */
  virtual void computeGeometricValues();
   /**
      Ce point d'entrée regroupe l'ensemble des calculs géométriques
   * utiles pour le schéma de projection (normal aux faces et centres des faces)
   */
  virtual void InitGeometricValues();
  
  /**
   * Calcule la nouvelle valeur de la densité des
   * mailles, en considérant que la masse d'une maille est constante au
   * cours du temps. Dans ce cas, la nouvelle densité est égale à la masse
   * divisée par le nouveau volume.
   */
  virtual void updateDensity();
		
  /**
   * Ce point d'entrée calcule l'énergie interne, la pression et la vitesse
   * du son dans la maille en faisant appel au service d'équation d'état.
   */
  virtual void updateEnergyAndPressure();

    /**
   * Ce point d'entrée calcule la pression moyenne dans la maille.
   */
  virtual void computePressionMoyenne();
		
  /**
   * Détermine la valeur du pas de temps pour l'itération suivante. 
   * Le pas de temps est contraint par :
   * - la valeur de la CFL,
   * - les valeurs \c deltatMin() et \c deltatMax() du jeu de données,
   * - la valeur du temps final. Lors de la dernière itération, le pas
   *   de temps doit être tel qu'on s'arrête exactement au temps spécifié
   *   dans le jeu de données (\c finalTime()).
   */
  virtual void computeDeltaT();
  
  /**
   * Calcul de quantites aux faces pour la projection :
   *    DxLagrange, du milieu, de la longueur des faces et de leur vitesse normale
   */
  virtual void computeFaceQuantitesForRemap();
  
  /**
   * Remplissage des variables de la projection et de la projection duale
   *         varlp->ULagrange (variables aux mailles)
   *                           de 0 à nbmat-1 : volume partiel,
   *                           de nbmat à 2*nbmat-1 : masse partielle
   *                           de 2*nbmat à 3*nbmat-1 : energie partielle
   *                           de 3*nbmat à 3*nbmat+1 : quantite de mouvement
   *                           3*nbmat+2 : enegie cinetique
   *                           3*nbmat+3 : pseudo-viscosite * volume
   *
   *         varlp->UDualLagrange (variables aux noeuds)
   *                           0 : masse
   *                           1 à 2 : quantite de mouvement
   *                           3 : energie cinetique
   *
   *  Pour l'option projection avec limiteurs pente-borne
   *
   *         varlp->Phi (variables aux mailles)
   *                           de 0 à nbmat-1 : fraction volumique
   *                           de nbmat à 2*nbmat-1 : densite partielle
   *                           de 2*nbmat à 3*nbmat-1 : energie specifique
   *partielle de 3*nbmat à 3*nbmat+1 : vitesse 3*nbmat+2 : enegie cinetique
   *specifique 3*nbmat+3 : pseudo-viscosite
   *
   *         varlp->DualPhi (variables aux noeuds)
   *                           0 : densite moyenne
   *                           1 à 2 : vitesse
   *                           3 : energie cinetique specifique
   */

  virtual void computeVariablesForRemap();
  virtual void remap();
  virtual void remapVariables();
  virtual void computeGradPhiFace(Integer idir, String name);
  virtual void computeGradPhiCell(Integer idir);
  virtual void computeUremap(Integer idir);
  virtual void synchronizeUremap();
  virtual void computeDualUremap(Integer idir, String name);
  virtual void synchronizeDualUremap();
  virtual void computeUpwindFaceQuantitiesForProjection(Integer idir, String name);
  virtual void computeAndLimitGradPhi(Integer projectionLimiterId, Face frontFace, Face backFace, Cell cell, Cell frontcell, Cell backcell);
  
  virtual void computeAndLimitGradPhiDual(Integer projectionLimiterId, 
                                          Node inode, Node frontnode, Node backnode,
                           Real3 grad_front, Real3 grad_back, Real h0, Real hplus, Real hmoins);
  
  virtual void computeDualGradPhi(Node inode, Node frontfrontnode, Node frontnode, 
                                   Node backnode, Node backbacknode, Integer idir);

  /** Retourne le numéro de version du module */
  virtual VersionInfo versionInfo() const { return VersionInfo(1,0,0); }
  
 private:
  
  /**
   * Calcule les résultantes aux noeuds d'une maille hexaédrique.
   * La méthode utilisée est celle du découpage en quatre triangles.
   * Méthode appelée par le point d'entrée \c computeGeometricValues()
   */
  inline void computeCQs(Real3 node_coord[8],Real3 face_coord[6],const Cell& cell);
  inline Real produit(Real A, Real B, Real C, Real D);
  inline Real norme(Real E, Real F, Real G);
  inline Real INTY(double X, double x0, double y0, double x1, double y1);
  inline Real fluxLimiter(Integer projectionLimiterId, double r);
  inline Real fluxLimiterG(Integer projectionLimiterId, double gradplus,
                           double gradmoins, double y0, double yplus,
                           double ymoins, double h0, double hplus,
                           double hmoins);
  inline Real computeY0(Integer projectionLimiterId, double y0, double yplus,
                        double ymoins, double h0, double hplus, double hmoins,
                        Integer type);
  inline Real computexgxd(double y0, double yplus, double ymoins, double h0,
                          double y0plus, double y0moins, Integer type);
  
  inline Real computeygyd(double y0, double yplus, double ymoins, double h0,
                          double y0plus, double y0moins, double grady,
                          Integer type);
  void computeFluxPP(Cell cell, Cell frontcell, Cell backcell, 
                                     Real face_normal_velocity, 
                                     Real deltat_n, Integer type, Real flux_threshold, 
                                     Integer projectionPenteBorneComplet, 
                                     Real dual_normal_velocity,
                                     Integer calcul_flux_dual,
                                     RealArrayView Flux, RealArrayView Flux_dual
                                    );
  void computeFluxPPPure(Cell cell, Cell frontcell, Cell backcell, 
                                     Real face_normal_velocity, 
                                     Real deltat_n, Integer type, Real flux_threshold, 
                                     Integer projectionPenteBorneComplet, 
                                     Real dual_normal_velocity,
                                     Integer calcul_flux_dual,
                                     RealArrayView Flux, RealArrayView Flux_dual
                                    );
  
  interval define_interval(double a, double b);
  interval intersection(interval I1, interval I2);
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
  Real computeRemapFlux(Integer projectionOrder, Integer projectionAvecPlateauPente,
        Real face_normal_velocity, Real3 face_normal,
        Real face_length, Real phi_face,
        Real3 outer_face_normal, Real3 exy, Real deltat_n);
  
  // fonctions d'initialisation des variables
  void hydroStartInitCasTest();
  void hydroStartInitVar();
  void initMatSOD();
  void initVarSOD();
  void initMatBiSOD();
  void initVarBiSOD();
  void initMatRider(Real3 Xb);
  void initMatRiderMono(Real3 Xb);
  void initVarRider(Real3 Xb);
  void initVarRiderMono(Real3 Xb);
  void initVarBiSedov();
  void initMatBiSedov();
  void initVarSedov();
  void initMatSedov();
  
  ICartesianMesh* m_cartesian_mesh;
  /* variables membre */
  Materials::IMeshMaterialMng* mm;
  Real m_deltat_n;
  Real m_deltat_nplus1;
  Integer m_nb_vars_to_project;
  Integer m_nb_env;
  Integer sens_projection;
  
  Integer my_rank;
  enum limiteur {
    minmod, 
    superBee,
    vanLeer ,
    minmodG ,
    superBeeG,
    vanLeerG,
    arithmeticG,
    ultrabeeG,
  }; 
  enum Test {
    // cas test
     UnitTestCase = 0,
     Sedov = 1,
     TriplePoint = 2,
     SodCaseX = 3,
     SodCaseY = 4,
     SodCaseZ = 5,
     NohTestCase = 6,
     AdvectionX = 7,
     AdvectionY = 8,
     AdvectionVitX = 9,
     AdvectionVitY = 10,
     BiSedov = 11,
     BiTriplePoint = 12,
     BiSodCaseX = 13,
     BiSodCaseY = 14,
     BiSodCaseZ = 15,
     BiShockBubble = 16,
     BiNohTestCase = 17,
     BiAdvectionX = 18,
     BiAdvectionY = 19,
     BiAdvectionVitX = 20,
     BiAdvectionVitY = 21,
     BiImplosion = 22,
     MonoRiderTx = 23,
     MonoRiderTy = 24,
     MonoRiderT45 = 25,
     MonoRiderRotation = 26,
     MonoRiderVortex = 27,
     MonoRiderDeformation = 28,
     MonoRiderVortexTimeReverse = 29,
     MonoRiderDeformationTimeReverse = 30,
     RiderTx = 31,
     RiderTy = 32,
     RiderT45 = 33,
     RiderRotation = 34,
     RiderVortex = 35,
     RiderDeformation = 36,
     RiderVortexTimeReverse = 37,
     RiderDeformationTimeReverse = 38,
  };
  const double Pi = 3.14159265359;

};

#endif
