// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef MAHYCOMODULE_H
#define MAHYCOMODULE_H

#include "TypesMahyco.h"

#include <arcane/geometry/IGeometryMng.h>
// Ajout au PIF
#include "arcane/utils/List.h"
#include "arcane/utils/OStringStream.h"
#include "arcane/utils/ValueChecker.h"

#include "arcane/utils/Simd.h"
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

#include "cartesian/interface/ICartesianMesh.h"
#include "cartesian/interface/CellDirectionMng.h"
#include "cartesian/interface/FaceDirectionMng.h"
#include "cartesian/interface/NodeDirectionMng.h"
#include "cartesian/interface/CartesianConnectivity.h"
// fin ajout au PIF

// Ajout pour accélérateur
#include "accenv/IAccEnv.h"
#include "accenv/AcceleratorUtils.h"
//

#include "Mahyco_axl.h"

// Pour les définitions, il faut finir par GCC car Clang et ICC définissent
// la macro __GNU__
// Pour CLANG, il n'y a pas encore d'équivalent au pragma ivdep.
// Celui qui s'en approche le plus est:
//   #pragma clang loop vectorize(enable)
// mais il ne force pas la vectorisation.
#ifdef __clang__
#  define PRAGMA_IVDEP_VALUE "clang loop vectorize(enable)"
#else
#  ifdef __INTEL_COMPILER
#    define PRAGMA_IVDEP_VALUE "ivdep"
#  else
#    ifdef __GNUC__
// S'assure qu'on compile avec la vectorisation même en '-O2'
#      pragma GCC optimize ("-ftree-vectorize")
#      define PRAGMA_IVDEP_VALUE "GCC ivdep"
#    endif
#  endif
#endif

//#undef PRAGMA_IVDEP_VALUE

#ifdef PRAGMA_IVDEP_VALUE
#define PRAGMA_IVDEP _Pragma(PRAGMA_IVDEP_VALUE)
#else
#define PRAGMA_IVDEP
#define PRAGMA_IVDEP_VALUE ""
#endif

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
  MahycoModule(const ModuleBuildInfo& mbi);
  /** Destructeur de la classe */
  ~MahycoModule();
  
  struct interval {
    double inf, sup;
  };

  // les paramètres pour appliquer les conditions aux limites sur des noeuds de bord
  struct BoundaryCondition
  {
    NodeGroup nodes; //!< le groupe de noeuds sur lequel s'applique la CL
    NodeVectorView boundary_nodes; //!< vue relative à ce groupe de noeuds
    Real value; //!< la valeur appliquée à la composante de vitesse
    TypesMahyco::eBoundaryCondition type; //!< le type de CL
  };

 public:

  //! Initialise l'environnement pour les accélérateurs
  void accBuild() override;

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
   *    c'est ARCANE qui les initialisent directement
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
   * Calcul générique de m_force et de v_velocity_out
   * Remarque : cette méthode ne peut pas être private ou protected
   * car elle déporte du calcul sur accélérateur
   */
  void updateForceAndVelocity(Real dt,
      const MaterialVariableCellReal& v_pressure,
      const MaterialVariableCellReal& v_pseudo_viscosity,
      const VariableCellArrayReal3& v_cell_cqs,
      const VariableNodeReal3& v_velocity_in,
      VariableNodeReal3& v_velocity_out);

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
  /*
   * Cacul de l'energie et de la pression par une méthode de Newton
   **/
  virtual void updateEnergyAndPressurebyNewton();  
  /*
   * Cacul de l'energie et de la pression par une méthode directe pour les gaz parfait
   **/
  virtual void updateEnergyAndPressureforGP();
    /**
   * Ce point d'entrée calcule la pression moyenne dans la maille.
   */
  virtual void computePressionMoyenne();
	
  /**
   * Calcul d'un pas de temps à partir des grandeurs hydrodynamiques
   * et affectation d'informations sur la maille qui fait le pas de temps
   *
   * \return Valeur du pas de temps hydro
   */
  template<typename DtCellInfoType>
  Real computeHydroDeltaT(DtCellInfoType &dt_cell_info);

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
  
  /**
   * Spécialisation de computeVariablesForRemap dans le cas penteborne=0
   * Est publique car fait appel à l'accélérateur
   **/
  void computeVariablesForRemap_PBorn0();
  
  /**
   * point d'entree pour la phase de projection
   **/
  virtual void remap();
 
  /** Retourne le numéro de version du module */
  virtual VersionInfo versionInfo() const { return VersionInfo(1,0,0); }
  
 private:
 
  /**
   * Calcule les résultantes aux noeuds d'une maille hexaédrique.
   * La méthode utilisée est celle du découpage en quatre triangles.
   * Méthode appelée par le point d'entrée \c computeGeometricValues()
   */
  ARCCORE_HOST_DEVICE inline void computeCQs(Real3 node_coord[8], Real3 face_coord[6], Span<Real3> out_cqs);
  
  // inline void computeCQsSimd(SimdReal3 node_coord[8],SimdReal3 face_coord[6],SimdReal3 cqs[8]);

  /**
   * A appeler après hydroStartInitEnvAndMat pour préparer
   * traitement des environnements sur accélérateur
   */
  void _initEnvForAcc();

  /** Les listes de faces XMIN, XMAX, YMIN ... doivent être construites au
   *  préalable par un appel à PrepareFaceGroup()
   */
  void _initBoundaryConditionsForAcc();

  /** Construit le maillage cartésien et les managers par direction
   */
  CartesianInterface::ICartesianMesh* _initCartMesh();

  /**
   * Fonctions diverses
   **/
  Real produit(Real A, Real B, Real C, Real D);
  inline Real norme(Real E, Real F, Real G);
  /** */
  /* la fonction dont on cherche un zero */
  /** */
  inline double fvnr(double e, double p, double dpde,
         double en, double qnn1, double pn, double rn1, double rn) 
         {return e-en+0.5*(p+pn+2.*qnn1)*(1./rn1-1./rn);};

  /** */
  /* la derivee de f */
  /** */
  inline double fvnrderiv(double e, double dpde, double rn1, double rn)
    {return 1.+0.5*dpde*(1./rn1-1./rn);};
    
   /** */
   /* la fonction dont on cherche le zero */
   /** */
   inline double f(double e, double p,double dpde,
		  double en,double qn, double pn, double cn1,
		  double cn, double m, double qn1, double cdn, double cdon, double qnm1) 
     {return e-en+0.5*(p+qn1)*cn1/m +0.5*(pn+qn)*cn/m-0.25*(pn+qn)*cdn/m + 0.5*(qn1-qnm1)*cdon/m;};
   /** */
   /* la derivee de f */
   /** */
   inline double fderiv(double e, double p, double dpde, double cn1, double m)
     {return 1.+0.5*dpde*cn1/m;}; 

  
  /* variables membre */
  // CartesianInterface:: = Arcane:: ou Cartesian::
  CartesianInterface::ICartesianMesh* m_cartesian_mesh;
  Materials::IMeshMaterialMng* mm;
  Integer m_nb_vars_to_project;
  Integer m_nb_env;
  Integer my_rank;
  Integer m_dimension;
 
  // Pour l'utilisation des accélérateurs
  IAccEnv* m_acc_env=nullptr;

  UniqueArray<BoundaryCondition> m_boundary_conditions;

  // Va contenir eosModel()->getAdiabaticCst(env), accessible à la fois sur CPU et GPU
  NumArray<Real,1> m_adiabatic_cst_env;
};

#endif
