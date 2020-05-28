#include <iostream>
#include <iomanip>
#include <type_traits>
#include <limits>
#include <utility>
#include <cmath>
#include <cfenv>
//#pragma STDC FENV_ACCESS ON

// Kokkos headers
#include <Kokkos_Core.hpp>
#include <Kokkos_hwloc.hpp>

// Project headers
#include "mesh/CartesianMesh2DGenerator.h"
#include "mesh/CartesianMesh2D.h"
#include "mesh/PvdFileWriter2D.h"
#include "utils/Utils.h"
#include "utils/Timer.h"
#include "types/Types.h"
#include "types/MathFunctions.h"
#include "types/ArrayOperations.h"

using namespace nablalib;

class EucclhydRemap
{
public:
  static const int dim=2;
  static const int nbmatmax=3;
  static const int nbequamax = 3*nbmatmax+2+1; // (volumes, masses, energies internes) * nbmatmax + vitesses + energie cinétique  
  struct Options
  {
    // Should be const but usefull to set them from main args
    RealArray1D<dim> ex = {{1.0, 0.0}};
    RealArray1D<dim> ey = {{0.0, 1.0}};
    RealArray1D<dim> zeroVect = {{0.0, 0.0}};
    RealArray2D<dim,dim> zeroMat = {{{0.0, 0.0}, {0.0, 0.0}}};
    RealArray1D<nbmatmax> zeroVectmat = {{0.0, 0.0, 0.0}};
    RealArray1D<nbequamax> Uzero = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    int UnitTestCase = 0;
    int SedovTestCase = 1;
    int TriplePoint = 2;
    int SodCase = 4;
    int NohTestCase = 5;
    int BiUnitTestCase = 10;
    int BiSedovTestCase = 11;
    int BiTriplePoint = 12;
    int BiShockBubble = 13;
    int BiSodCase = 14;
    int BiNohTestCase = 15;
    int eosPerfectGas = 100;
    int symmetry = 200;
    int imposedVelocity = 201;
    int freeSurface = 202;
    int minmod = 300;
    int superBee = 301;
    int vanLeer = 302;
    int minmodG = 1300;
    int superBeeG = 1301;
    int vanLeerG = 1302;
    int arithmeticG = 1303;
    int testCase = SedovTestCase;
    double final_time = 1.0;
    double output_time = final_time;
    double gamma = 1.4;
    RealArray1D<nbmatmax> gammap = {{1.4, 1.4, 1.4}};
    double cfl = 0.45 ;
    double X_LENGTH = 1.2;
    double Y_LENGTH = X_LENGTH;
    int X_EDGE_ELEMS = 30;
    int Y_EDGE_ELEMS = 30;
    double X_EDGE_LENGTH = X_LENGTH / X_EDGE_ELEMS;
    double Y_EDGE_LENGTH = Y_LENGTH / Y_EDGE_ELEMS;
    int max_time_iterations = 500000000;
    double u0 = 0.0;
    double p0 = 1.0;
    double rho0 = 1.0;
    double threshold = 1.0E-16;
    double deltat_init = 0.;
    double deltat_min = 1.0E-10;
    int eos = eosPerfectGas;
    int spaceOrder = 2;
    int projectionOrder = 2;
    int projectionLimiterId = superBee;
    int projectionLimiterIdPure = arithmeticG;
    int projectionAvecPlateauPente = 0;
    int projectionConservative = 0;
    int projectionLimiteurMixte = 0;

    int leftFluxBC=0;
    RealArray1D<nbequamax> leftFluxBCValue = Uzero;
    int rightFluxBC=0;
    RealArray1D<nbequamax> rightFluxBCValue = Uzero;
    int bottomFluxBC=0;
    RealArray1D<nbequamax> bottomFluxBCValue = Uzero;
    int topFluxBC=0;
    RealArray1D<nbequamax> topFluxBCValue = Uzero;
    int FluxBC = leftFluxBC+rightFluxBC+bottomFluxBC+topFluxBC;
    
    int leftBC = 0;
    RealArray1D<dim> leftBCValue = ey;
	  	  
    int rightBC = 0;
    RealArray1D<dim> rightBCValue = ey;
	  
    int topBC = 0;
    RealArray1D<dim> topBCValue = ex;
	  
    int bottomBC = 0;
    RealArray1D<dim> bottomBCValue = ex;

    int DragModel;
    int Kliatchko = 20;
    int Classique = 21;
    int KliatchkoDragModel = 20;
    
    double Reynolds_min = 1.e-4;
    double Reynolds_max = 1.e3;
    double Drag = 10.;
  };
  Options* options;

private:
  
  CartesianMesh2D* mesh;
  PvdFileWriter2D writer;
  PvdFileWriter2D writerpart;
  int nbPartMax;
  int nbPart=0;
  int nbNodes, nbCells, nbFaces, nbCellsOfNode, nbFacesOfCell, nbNodesOfCell, nbInnerNodes, nbNodesOfFace, nbInnerVerticalFaces, nbInnerHorizontalFaces, nbNeighbourCells, nbCellsOfFace, nbLeftNodes, nbRightNodes, nbTopNodes, nbBottomNodes, nbTopLeftNode, nbTopRightNode, nbBottomLeftNode, nbBottomRightNode;

  // constantes
  double Pi = 3.14159265359;
  double viscosity = 1.715e-5;
  double Cp = 2.84;
  double Pr = 0.7;
  // Global Variables
  int n, nbCalls;
  bool x_then_y_n, x_then_y_nplus1;
  double t_n, t_nplus1, deltat_n, deltat_nplus1, lastDump;
  int imat;
  double ETOTALE_L, ETOTALE_T, ETOTALE_0;
  double MASSET_L, MASSET_T, MASSET_0;

        // cells a debuguer
        int dbgcell1=-389;
        int dbgcell2=-126;
        int dbgcell3=-156;
        int face_debug1=-2033;
        int face_debug2=-410;
        int test_debug=1;

	// Connectivity Variables
	Kokkos::View<RealArray1D<dim>*> X;
	Kokkos::View<RealArray1D<dim>*> Xc;
	Kokkos::View<RealArray1D<dim>*> Xf;
	Kokkos::View<RealArray1D<dim>*> XfLagrange;
	Kokkos::View<RealArray1D<dim>*> XLagrange;
	Kokkos::View<RealArray1D<dim>*> XcLagrange;
	Kokkos::View<double*> Xc_x;
        Kokkos::View<double*> Xc_y;
	Kokkos::View<RealArray1D<dim>**> lpc_n;
	Kokkos::View<RealArray1D<dim>**> nplus;
	Kokkos::View<RealArray1D<dim>**> nminus;
	Kokkos::View<double**> lplus;
	Kokkos::View<double**> lminus;
	Kokkos::View<double*> p;
        Kokkos::View<RealArray1D<nbmatmax>*> pp;
	Kokkos::View<double*> m;
        Kokkos::View<RealArray1D<nbmatmax>*> mp;
	Kokkos::View<double*> v;
	Kokkos::View<double*> vLagrange;
	Kokkos::View<double*> LfLagrange;
	Kokkos::View<double*> HvLagrange;
	Kokkos::View<RealArray1D<nbmatmax>*> vpLagrange;
	Kokkos::View<double*> perim;
	Kokkos::View<double*> vitson;
	Kokkos::View<RealArray1D<nbmatmax>*> vitsonp;
	Kokkos::View<double*> deltaxLagrange;
	Kokkos::View<double*> faceLength;
	Kokkos::View<double*> faceLengthLagrange;
	Kokkos::View<double*> faceNormalVelocity;
	Kokkos::View<RealArray1D<dim>*> faceNormal;
	Kokkos::View<RealArray1D<dim>**> outerFaceNormal;
	Kokkos::View<double*> deltatc;
	Kokkos::View<RealArray1D<dim>*> Vnode_n;
	Kokkos::View<RealArray1D<dim>*> Vnode_nplus1;
	Kokkos::View<RealArray1D<dim>*> Vnode_n0;  
	Kokkos::View<double*> ETOT_0;
	Kokkos::View<double*> ETOT_T;
	Kokkos::View<double*> ETOT_L; 
	Kokkos::View<double*> MTOT_0;
	Kokkos::View<double*> MTOT_T;
	Kokkos::View<double*> MTOT_L;
	Kokkos::View<double*> rho_n;
	Kokkos::View<double*> rho_nplus1;
	Kokkos::View<double*> rho_n0;
	Kokkos::View<RealArray1D<nbmatmax>*> rhop_n;
	Kokkos::View<RealArray1D<nbmatmax>*> rhop_nplus1;
	Kokkos::View<RealArray1D<nbmatmax>*> rhop_n0;
	Kokkos::View<RealArray1D<dim>*> V_n;
	Kokkos::View<RealArray1D<dim>*> V_nplus1;
	Kokkos::View<RealArray1D<dim>*> V_n0;
	Kokkos::View<double*> Vxc;
	Kokkos::View<double*> Vyc;
	Kokkos::View<double*> eps_n;
	Kokkos::View<double*> eps_nplus1;
	Kokkos::View<double*> eps_n0;  
	Kokkos::View<double*> delta_ec;  
	Kokkos::View<RealArray1D<nbmatmax>*> epsp_n;
	Kokkos::View<RealArray1D<nbmatmax>*> epsp_nplus1;
	Kokkos::View<RealArray1D<nbmatmax>*> epsp_n0;
	Kokkos::View<RealArray1D<nbequamax>*> ULagrange;
	Kokkos::View<RealArray1D<nbequamax>*> Uremap1;
	Kokkos::View<RealArray1D<nbequamax>*> Uremap2;
	Kokkos::View<RealArray1D<nbequamax>*> gradPhiFace1;
	Kokkos::View<RealArray1D<nbequamax>*> gradPhiFace2;
	Kokkos::View<RealArray1D<nbequamax>*> gradPhi1;
	Kokkos::View<RealArray1D<nbequamax>*> gradPhi2;
	Kokkos::View<RealArray1D<nbequamax>*> phiFace1;
	Kokkos::View<RealArray1D<nbequamax>*> phiFace2;
	Kokkos::View<RealArray1D<nbequamax>*> deltaPhiFaceAv;
	Kokkos::View<RealArray1D<nbequamax>*> deltaPhiFaceAr;
	Kokkos::View<RealArray1D<nbequamax>*> Phi;
	Kokkos::View<double**> p_extrap;
        Kokkos::View<RealArray1D<nbmatmax>**> pp_extrap;
	Kokkos::View<RealArray1D<dim>**> V_extrap;
	Kokkos::View<RealArray1D<dim>*> gradp;
        Kokkos::View<RealArray1D<dim>*> gradp1;
        Kokkos::View<RealArray1D<dim>*> gradp2;
        Kokkos::View<RealArray1D<dim>*> gradp3;
        Kokkos::View<RealArray1D<dim>*> gradf1;
        Kokkos::View<RealArray1D<dim>*> gradf2;
        Kokkos::View<RealArray1D<dim>*> gradf3;
	Kokkos::View<RealArray2D<dim,dim>*> gradV;
	Kokkos::View<RealArray1D<dim>**> F_n;
	Kokkos::View<RealArray1D<dim>**> F_nplus1;
	Kokkos::View<RealArray1D<dim>**> F_n0;
	Kokkos::View<RealArray1D<dim>**> F1_n;
	Kokkos::View<RealArray1D<dim>**> F1_nplus1;
	Kokkos::View<RealArray1D<dim>**> F2_n;
	Kokkos::View<RealArray1D<dim>**> F2_nplus1;
	Kokkos::View<RealArray1D<dim>**> F3_n;
	Kokkos::View<RealArray1D<dim>**> F3_nplus1;
	Kokkos::View<RealArray1D<dim>*> G;
	Kokkos::View<RealArray2D<dim,dim>**> M;
	Kokkos::View<RealArray2D<dim,dim>**> M1;
	Kokkos::View<RealArray2D<dim,dim>**> M2;
	Kokkos::View<RealArray2D<dim,dim>**> M3;
	Kokkos::View<RealArray2D<dim,dim>*> Mnode;
        Kokkos::View<RealArray1D<nbmatmax>*> fracmass;
        Kokkos::View<RealArray1D<nbmatmax>*> fracvol;
        Kokkos::View<RealArray1D<nbmatmax>*> fracvolnode;
	Kokkos::View<double*> fracvol1;
	Kokkos::View<double*> fracvol2;
	Kokkos::View<double*> fracvol3;
        
	Kokkos::View<double*> vpart;
        Kokkos::View<double*> wpart;
	Kokkos::View<double*> mpart;
	Kokkos::View<double*> rpart;
	Kokkos::View<double*> rhopart;
        Kokkos::View<double*> Cdpart;
        Kokkos::View<double*> Mcpart;
        Kokkos::View<double*> Repart;
        Kokkos::View<double*> Temppart;
        Kokkos::View<RealArray1D<dim>*> Xpart_n0;
        Kokkos::View<RealArray1D<dim>*> Xpart_n;
        Kokkos::View<RealArray1D<dim>*> Xpart_nplus1;
        Kokkos::View<RealArray1D<dim>*> ForceGradp;
	Kokkos::View<RealArray1D<dim>*> Vpart_n0;
	Kokkos::View<RealArray1D<dim>*> Vpart_n;
	Kokkos::View<RealArray1D<dim>*> Vpart_nplus1;
        Kokkos::View<int*> mixte; 
        Kokkos::View<int*> pure; 
        Kokkos::View<int*> ICellp;
        Kokkos::View<int*> IMatp;
	Kokkos::View<double*> fracpart;
        Kokkos::View<vector<int>*> listpart;
  
	utils::Timer global_timer;
	utils::Timer cpu_timer;
	utils::Timer io_timer;
        const size_t maxHardThread = Kokkos::DefaultExecutionSpace::max_hardware_threads();

public:
	EucclhydRemap(Options* aOptions, CartesianMesh2D* aCartesianMesh2D, string output)
	: options(aOptions)
	, mesh(aCartesianMesh2D)
	, writer("EucclhydRemap", output)
	, writerpart("Particules", output)
	, nbNodes(mesh->getNbNodes())
	, nbPartMax(1)
	, nbCells(mesh->getNbCells())
	, nbFaces(mesh->getNbFaces())
	, nbCellsOfNode(CartesianMesh2D::MaxNbCellsOfNode)
	, nbFacesOfCell(CartesianMesh2D::MaxNbFacesOfCell)
	, nbNodesOfCell(CartesianMesh2D::MaxNbNodesOfCell)
	, nbInnerNodes(mesh->getNbInnerNodes())
	, nbNodesOfFace(CartesianMesh2D::MaxNbNodesOfFace)
	, nbInnerVerticalFaces(mesh->getNbInnerVerticalFaces())
	, nbInnerHorizontalFaces(mesh->getNbInnerHorizontalFaces())
	, nbNeighbourCells(CartesianMesh2D::MaxNbNeighbourCells)
	, nbCellsOfFace(CartesianMesh2D::MaxNbCellsOfFace)
	, nbLeftNodes(mesh->getNbLeftNodes())
	, nbRightNodes(mesh->getNbRightNodes())
	, nbTopNodes(mesh->getNbTopNodes())
	, nbBottomNodes(mesh->getNbBottomNodes())
	, nbTopLeftNode(mesh->getNbTopLeftNode())
	, nbTopRightNode(mesh->getNbTopRightNode())
	, nbBottomLeftNode(mesh->getNbBottomLeftNode())
	, nbBottomRightNode(mesh->getNbBottomRightNode())
	, x_then_y_n(true)
	, x_then_y_nplus1(true)
	, t_n(0.0)
	, t_nplus1(0.0)
	, deltat_n(options->deltat_init)
	, deltat_nplus1(options->deltat_init)
	, nbCalls(0)
	, lastDump(0.0)
	, vpart("VolumePart", nbPartMax)
	, wpart("WeightPart", nbPartMax)
	, mpart("MassePart", nbPartMax)
	, rpart("RayonPart", nbPartMax)
	, rhopart("RhoPart", nbPartMax)
	, Repart("RePart", nbPartMax)
	, Cdpart("Cdpart", nbPartMax)
	, Mcpart("Cdpart", nbPartMax)
	, ForceGradp("ForceGradp", nbCells)
	, listpart("listepart", nbCells)
	, fracpart("fracPart", nbCells)
	, Xpart_n0("Xpart_n0", nbPartMax)
	, Xpart_n("Xpart_n", nbPartMax)
	, Xpart_nplus1("Xpart_nplus1", nbPartMax)	  
	, Vpart_n0("Vpart_n0", nbPartMax)
	, Vpart_n("Vpart_n", nbPartMax)
	, Vpart_nplus1("Vpart_nplus1", nbPartMax)
	, ICellp("Icellp", nbPartMax)
	, IMatp("Imatp", nbPartMax)
	, X("X", nbNodes)
	, Xc("Xc", nbCells)
	, Xf("Xf", nbFaces)
	, XfLagrange("Xf", nbFaces)
	, LfLagrange("LfLagrange", nbFaces)
	, XLagrange("XLagrange", nbNodes)
	, XcLagrange("XcLagrange", nbCells)
	, HvLagrange("HvLagrange", nbCells)
	, Xc_x("Xc_x", nbCells)
	, Xc_y("Xc_y", nbCells)
	, lpc_n("lpc_n", nbNodes, nbCellsOfNode)
	, nplus("nplus", nbNodes, nbCellsOfNode)
	, nminus("nminus", nbNodes, nbCellsOfNode)
	, lplus("lplus", nbNodes, nbCellsOfNode)
	, lminus("lminus", nbNodes, nbCellsOfNode)
	, p("p", nbCells)
	  //, p1("p1", nbCells) à réactiver pour les sorties
	  //, p2("p2", nbCells)
	, pp("pp", nbCells)	
	, m("m", nbCells)	
	, mp("mp", nbCells)	
	  //, m1("m1", nbCells)	à réactiver pour les sorties
	  //, m2("m2", nbCells)
	, v("v", nbCells)
	, fracmass("fracmass", nbCells)
	  //, c1("c1", nbCells) à réactiver pour les sorties
	  //, c2("c2", nbCells)
	  //, c3("c3", nbCells)
	, mixte("mixte", nbCells)
	, pure("pure", nbCells)
	, fracvol("fracvol", nbCells)
	, fracvolnode("fracvolnode", nbNodes)
	, fracvol1("fracvol1", nbCells)
	, fracvol2("fracvol2", nbCells)
	, fracvol3("fracvol3", nbCells)
	, vLagrange("vLagrange", nbCells)
	, vpLagrange("vpLagrange", nbCells)
	, perim("perim", nbCells)
	, vitson("vitson", nbCells)
	, vitsonp("vitsonp", nbCells)
	, deltaxLagrange("deltaxLagrange", nbFaces)
	, faceLength("faceLength", nbFaces)
	, faceLengthLagrange("faceLengthLagrange", nbFaces)
	, faceNormalVelocity("faceNormalVelocity", nbFaces)
	, faceNormal("faceNormal", nbFaces)
	, outerFaceNormal("outerFaceNormal", nbCells, nbFacesOfCell)
	, deltatc("deltatc", nbCells)
	, Vnode_n("Vnode_n", nbNodes)
	, Vnode_nplus1("Vnode_nplus1", nbNodes)
	, Vnode_n0("Vnode_n0", nbNodes)
	, ETOT_0("ETOT_0", nbCells)
	, ETOT_T("ETOT_T", nbCells)
	, ETOT_L("ETOT_L", nbCells)
	, MTOT_0("MTOT_0", nbCells)
	, MTOT_T("MTOT_T", nbCells)
	, MTOT_L("MTOT_L", nbCells)
	, rho_n("rho_n", nbCells)
	, rho_nplus1("rho_nplus1", nbCells)
	, rho_n0("rho_n0", nbCells)
	, rhop_n("rhop_n", nbCells)
	, rhop_nplus1("rhop_nplus1", nbCells)
	, rhop_n0("rhop_n0", nbCells)
	, V_n("V_n", nbCells)
	, V_nplus1("V_nplus1", nbCells)
	, V_n0("V_n0", nbCells)
	, Vxc("Vxc", nbCells)
	, Vyc("Vyc", nbCells)
	, eps_n("eps_n", nbCells)
	, eps_nplus1("eps_nplus1", nbCells)
	, eps_n0("eps_n0", nbCells)
	, epsp_n("epsp_n", nbCells)
	, epsp_nplus1("epsp_nplus1", nbCells)
	, epsp_n0("epsp_n0", nbCells)
	, delta_ec("delta_ec", nbCells)
	, ULagrange("ULagrange", nbCells)
	, Uremap1("Uremap1", nbCells)
	, Uremap2("Uremap2", nbCells)
	, gradPhiFace1("gradPhiFace1", nbFaces)
	, gradPhiFace2("gradPhiFace2", nbFaces)
	, gradPhi1("gradPhi1", nbCells)
	, gradPhi2("gradPhi2", nbCells)
	, phiFace1("phiFace1", nbFaces)
	, phiFace2("phiFace2", nbFaces)
	, deltaPhiFaceAv("deltaPhiFaceAv", nbCells)
	, deltaPhiFaceAr("deltaPhiFaceAr", nbCells)
	, Phi("Phi", nbCells)
	, p_extrap("p_extrap", nbCells, nbNodesOfCell)
	  //, p1_extrap("p1_extrap", nbCells, nbNodesOfCell)
	  //, p2_extrap("p2_extrap", nbCells, nbNodesOfCell)
	, pp_extrap("pp_extrap", nbCells, nbNodesOfCell)
	, V_extrap("V_extrap", nbCells, nbNodesOfCell)
	, gradp("gradp", nbCells)
	, gradp1("gradp1", nbCells)
	, gradp2("gradp2", nbCells)
	, gradp3("gradp3", nbCells)
	, gradf1("gradf1", nbCells)
	, gradf2("gradf2", nbCells)
	, gradf3("gradf3", nbCells)
	, gradV("gradV", nbCells)
	, F_n("F_n", nbNodes, nbCellsOfNode)
	, F_nplus1("F_nplus1", nbNodes, nbCellsOfNode)
	, F_n0("F_n0", nbNodes, nbCellsOfNode)
	, F1_n("F1_n", nbNodes, nbCellsOfNode)
	, F1_nplus1("F1_nplus1", nbNodes, nbCellsOfNode)
	, F2_n("F2_n", nbNodes, nbCellsOfNode)
	, F2_nplus1("F2_nplus1", nbNodes, nbCellsOfNode)
	, F3_n("F3_n", nbNodes, nbCellsOfNode)
	, F3_nplus1("F3_nplus1", nbNodes, nbCellsOfNode)
	, G("G", nbNodes)
	, M("M", nbNodes, nbCellsOfNode)
	, M1("M1", nbNodes, nbCellsOfNode)
	, M2("M2", nbNodes, nbCellsOfNode)
	, M3("M3", nbNodes, nbCellsOfNode)
	, Mnode("Mnode", nbNodes)
	{
		// Copy node coordinates
		const auto& gNodes = mesh->getGeometry()->getNodes();
		Kokkos::parallel_for(nbNodes, KOKKOS_LAMBDA(const int& rNodes)
		{
			X(rNodes) = gNodes[rNodes];
		});
	}

private:
  #include "Init.cc"
  #include "PhaseLagrange.cc"
  #include "PhaseRemap1.cc"
  #include "PhaseRemap2.cc"
  #include "PhaseRemapFinal.cc"
  #include "Utiles.cc"
  #include "UtilesRemap.cc"
  #include "ConditionsLimites.cc"
  #include "SchemaParticules.cc"
  
  /**
   * Job dumpVariables called @2.0 in executeTimeLoopN method.
   * In variables: Xc_x, Xc_y, eps_n, m, p, rho_n, t_n, v
   * Out variables: 
   */
  KOKKOS_INLINE_FUNCTION
  void dumpVariables() noexcept
  {
	  
    std::cout << " Deltat = " << deltat_n << std::endl;
    std::cout << " ---------------------------" << " Energie totale(t=0) = " << ETOTALE_0 << " Energie totale(Lagrange) = " << ETOTALE_L << " Energie totale(time) = " << ETOTALE_T << std::endl;
    std::cout << " ---------------------------" << " Masse totale(t=0) = " << MASSET_0 << " Masse totale(Lagrange) = " << MASSET_L << " Masse totale(time) = " << MASSET_T << std::endl;
    nbCalls++;
    if (!writer.isDisabled() && (t_n >= lastDump + options->output_time || t_n == 0. ))
      {
	cpu_timer.stop();
	io_timer.start();
	std::map<string, double*> cellVariables;
	std::map<string, double*> nodeVariables;
	std::map<string, double*> partVariables;
	cellVariables.insert(pair<string,double*>("Pressure", p.data()));
	cellVariables.insert(pair<string,double*>("Density", rho_n.data()));
	cellVariables.insert(pair<string,double*>("F1", fracvol1.data()));
	cellVariables.insert(pair<string,double*>("F2", fracvol2.data()));
	cellVariables.insert(pair<string,double*>("F3", fracvol3.data()));
	cellVariables.insert(pair<string,double*>("VelocityX", Vxc.data()));
	cellVariables.insert(pair<string,double*>("VelocityY", Vyc.data()));
	cellVariables.insert(pair<string,double*>("Energy", eps_n.data()));
	partVariables.insert(pair<string,double*>("VolumePart", vpart.data()));
	partVariables.insert(pair<string,double*>("VxPart", Vpart_n[0].data()));
	partVariables.insert(pair<string,double*>("VyPart", Vpart_n[1].data()));
	auto quads = mesh->getGeometry()->getQuads();
	writer.writeFile(nbCalls, t_n, nbNodes, X.data(), nbCells, quads.data(), cellVariables, nodeVariables);
	writerpart.writeFile(nbCalls, t_n, nbPart, Xpart_n.data(), 0, quads.data(), cellVariables, partVariables);
	lastDump = t_n;
	std::cout << " time = " << t_n << " sortie demandée " << std::endl;
	io_timer.stop();
	cpu_timer.start();
      }
  }

  /**
   * Job executeTimeLoopN called @4.0 in simulate method.
   * In variables: F_n, F_nplus1, G, M, Mnode, ULagrange, Uremap1, Uremap2, V_extrap, V_n, Vnode_n, Vnode_nplus1, X, XLagrange, Xc, XcLagrange, Xc_x, Xc_y, Xf, bottomBC, bottomBCValue, c, cfl, deltat_n, deltat_nplus1, deltatc, deltaxLagrange, eos, eosPerfectGas, eps_n, faceLength, faceNormal, faceNormalVelocity, gamma, gradPhi1, gradPhi2, gradPhiFace1, gradPhiFace2, gradV, gradp, leftBC, leftBCValue, lminus, lpc_n, lplus, m, nminus, nplus, outerFaceNormal, p, p_extrap, perim, phiFace1, phiFace2, projectionLimiterId, projectionOrder, rho_n, rightBC, rightBCValue, spaceOrder, t_n, topBC, topBCValue, v, vLagrange, x_then_y_n
   * Out variables: F_nplus1, G, M, Mnode, ULagrange, Uremap1, Uremap2, V_extrap, V_nplus1, Vnode_nplus1, XLagrange, XcLagrange, c, deltat_nplus1, deltatc, deltaxLagrange, eps_nplus1, faceNormalVelocity, gradPhi1, gradPhi2, gradPhiFace1, gradPhiFace2, gradV, gradp, m, p, p_extrap, phiFace1, phiFace2, rho_nplus1, t_nplus1, vLagrange, x_then_y_nplus1
   */
  KOKKOS_INLINE_FUNCTION
  void executeTimeLoopN() noexcept
  {
    n = 0;
    bool continueLoop = true;
    do
      {
	global_timer.start();
	cpu_timer.start();
	n++;
	if (n!=1)
	  std::cout << "[" << __CYAN__ << __BOLD__ << setw(3) << n << __RESET__ "] time = " << __BOLD__
		    << setiosflags(std::ios::scientific) << setprecision(8) << setw(16) << t_n << __RESET__;
	
	if (nbPart !=0) switchalpharho_rho();		
	computeEOS(); // @1.0
	if (nbPart !=0) switchrho_alpharho();
	computeGradients(); // @1.0
	computeMass(); // @1.0		
	computeDissipationMatrix(); // @2.0
	computedeltatc(); // @2.0		
	dumpVariables(); // @2.0
	extrapolateValue(); // @2.0
	computeG(); // @3.0
	computeNodeDissipationMatrixAndG(); // @3.0
	computedeltat(); // @3.0
	computeBoundaryNodeVelocities(); // @4.0
	computeNodeVelocity(); // @4.0
	updateTime(); // @4.0
	computeFaceVelocity(); // @5.0
	computeLagrangePosition(); // @5.0
	computeSubCellForce(); // @5.0
	computeLagrangeVolumeAndCenterOfGravity(); // @6.0
	computeFacedeltaxLagrange(); // @7.0			
	updateCellCenteredLagrangeVariables(); // @7.0
	computeGradPhiFace1(); // @8.0
	computeGradPhi1(); // @9.0
	computeUpwindFaceQuantitiesForProjection1(); // @10.0
	computeUremap1(); // @11.0			
	computeGradPhiFace2(); // @12.0
	computeGradPhi2(); // @13.0
	computeUpwindFaceQuantitiesForProjection2(); // @14.0
	computeUremap2(); // @15.0
	remapCellcenteredVariable(); // @16.0
	if (nbPart !=0) {
	  updateParticlePosition();
	  updateParticleCoefficients();
	  updateParticleVelocity();
	  updateParticleRetroaction();
	}
		
	// Evaluate loop condition with variables at time n
	continueLoop = (n + 1 < options->max_time_iterations && t_nplus1 < options->final_time);
		
	if (continueLoop)
	  {
	    // Switch variables to prepare next iteration
	    std::swap(x_then_y_nplus1, x_then_y_n);
	    std::swap(t_nplus1, t_n);
	    std::swap(deltat_nplus1, deltat_n);
	    std::swap(Vnode_nplus1, Vnode_n);
	    std::swap(rho_nplus1, rho_n);
	    std::swap(rhop_nplus1, rhop_n);
	    std::swap(V_nplus1, V_n);
	    std::swap(eps_nplus1, eps_n);
	    std::swap(epsp_nplus1, epsp_n);
	    std::swap(F_nplus1, F_n);
	    std::swap(F1_nplus1, F1_n);
	    std::swap(F2_nplus1, F2_n);
	    std::swap(F3_nplus1, F3_n);
            std::swap(Vpart_nplus1, Vpart_n);
            std::swap(Xpart_nplus1, Xpart_n);
	  }
		
	cpu_timer.stop();
	global_timer.stop();
		
	// Timers display
	if (!writer.isDisabled())
	  std::cout << " {CPU: " << __BLUE__ << cpu_timer.print(true) << __RESET__ ", IO: " << __BLUE__ << io_timer.print(true) << __RESET__ "} ";
	else
	  std::cout << " {CPU: " << __BLUE__ << cpu_timer.print(true) << __RESET__ ", IO: " << __RED__ << "none" << __RESET__ << "} ";
			
	// Progress
	std::cout << utils::progress_bar(n, options->max_time_iterations, t_n, options->final_time, 30);
	std::cout << __BOLD__ << __CYAN__ << utils::Timer::print(
								 utils::eta(n, options->max_time_iterations, t_n, options->final_time, deltat_n, global_timer), true)
		  << __RESET__ << "\r";
	std::cout.flush();
		
	cpu_timer.reset();
	io_timer.reset();
      } while (continueLoop);
  }
  /**
   * Job computedeltat called @3.0 in executeTimeLoopN method.
   * In variables: cfl, deltat_n, deltatc
   * Out variables: deltat_nplus1
   */
  KOKKOS_INLINE_FUNCTION
  void computedeltat() noexcept
  {
    double reduction10(numeric_limits<double>::max());
    {	
      Kokkos::Min<double> reducer(reduction10);
      Kokkos::parallel_reduce("reduction10", nbCells, KOKKOS_LAMBDA(const int& cCells, double& x)
			      {
				reducer.join(x, deltatc(cCells));
			      }, reducer);
    }
    deltat_nplus1 = MathFunctions::min(options->cfl * reduction10, deltat_n * 1.05);
    if (deltat_nplus1 < options->deltat_min)
      {
	std::cerr << "Fin de la simulation par pas de temps minimum " << deltat_nplus1 << " < " << options->deltat_min << std::endl;		    
	Kokkos::finalize();
	exit(1);
      }		    
  }
  /**
   * Job updateTime called @4.0 in executeTimeLoopN method.
   * In variables: deltat_nplus1, t_n
   * Out variables: t_nplus1
   */
  KOKKOS_INLINE_FUNCTION
  void updateTime() noexcept
  {
    t_nplus1 = t_n + deltat_nplus1;
  }



public:
	void simulate()
	{
		std::cout << "\n" << __BLUE_BKG__ << __YELLOW__ << __BOLD__ <<"\tStarting EucclhydRemap ..." << __RESET__ << "\n\n";
		
		std::cout << "[" << __GREEN__ << "MESH" << __RESET__ << "]      X=" << __BOLD__ << options->X_EDGE_ELEMS << __RESET__ << ", Y=" << __BOLD__ << options->Y_EDGE_ELEMS
			<< __RESET__ << ", X length=" << __BOLD__ << options->X_EDGE_LENGTH << __RESET__ << ", Y length=" << __BOLD__ << options->Y_EDGE_LENGTH << __RESET__ << std::endl;
		
		if (Kokkos::hwloc::available()) {
			std::cout << "[" << __GREEN__ << "TOPOLOGY" << __RESET__ << "]  NUMA=" << __BOLD__ << Kokkos::hwloc::get_available_numa_count()
				<< __RESET__ << ", Cores/NUMA=" << __BOLD__ << Kokkos::hwloc::get_available_cores_per_numa()
				<< __RESET__ << ", Threads/Core=" << __BOLD__ << Kokkos::hwloc::get_available_threads_per_core() << __RESET__ << std::endl;
		} else {
			std::cout << "[" << __GREEN__ << "TOPOLOGY" << __RESET__ << "]  HWLOC unavailable cannot get topological informations" << std::endl;
		}
		
		// std::cout << "[" << __GREEN__ << "KOKKOS" << __RESET__ << "]    " << __BOLD__ << (is_same<MyLayout,Kokkos::LayoutLeft>::value?"Left":"Right")" << __RESET__ << " layout" << std::endl;
		
		if (!writer.isDisabled())
			std::cout << "[" << __GREEN__ << "OUTPUT" << __RESET__ << "]    VTK files stored in " << __BOLD__ << writer.outputDirectory() << __RESET__ << " directory" << std::endl;
		else
			std::cout << "[" << __GREEN__ << "OUTPUT" << __RESET__ << "]    " << __BOLD__ << "Disabled" << __RESET__ << std::endl;

		computeCornerNormal(); // @1.0
		initMeshGeometryForCells(); // @1.0
		initVpAndFpc(); // @1.0
		initBoundaryConditions();
		initCellInternalEnergy(); // @2.0
		initCellVelocity(); // @2.0		
		initDensity(); // @2.0	
		initMeshGeometryForFaces(); // @2.0
		if (nbPart !=0) {
		  initPart();
		  updateParticleCoefficients();
		  switchrho_alpharho(); // on travaille avec alpharho sauf pour l'EOS
		}  
		setUpTimeLoopN(); // @3.0
		executeTimeLoopN(); // @4.0
		std::cout << __YELLOW__ << "\n\tDone ! Took " << __MAGENTA__ << __BOLD__ << global_timer.print() << __RESET__ << std::endl;
	}
};

int main(int argc, char* argv[]) 
{
	Kokkos::initialize(argc, argv);
	auto o = new EucclhydRemap::Options();
	string output;
	if (argc == 12) {
		o->X_EDGE_ELEMS = std::atoi(argv[1]);
		o->Y_EDGE_ELEMS = std::atoi(argv[2]);
		o->X_EDGE_LENGTH = std::atof(argv[3]);
		o->Y_EDGE_LENGTH = std::atof(argv[4]);
		o->testCase = std::atof(argv[5]);
		o->final_time = std::atof(argv[6]);
		o->output_time = std::atof(argv[7]);
		o->projectionLimiterId = std::atof(argv[8]);
		o->projectionAvecPlateauPente = std::atof(argv[9]);
		o->projectionLimiteurMixte = std::atof(argv[10]);
		o->projectionLimiterIdPure = std::atof(argv[11]);
	} else if (argc == 13) {
		o->X_EDGE_ELEMS = std::atoi(argv[1]);
		o->Y_EDGE_ELEMS = std::atoi(argv[2]);
		o->X_EDGE_LENGTH = std::atof(argv[3]);
		o->Y_EDGE_LENGTH = std::atof(argv[4]);
		o->testCase = std::atof(argv[5]);
		o->final_time = std::atof(argv[6]);
		o->output_time = std::atof(argv[7]);
		o->projectionLimiterId = std::atof(argv[8]);
		o->projectionAvecPlateauPente = std::atof(argv[9]);
		o->projectionLimiteurMixte = std::atof(argv[10]);		
		o->projectionLimiterIdPure = std::atof(argv[11]);
		output = argv[11];
	} else if (argc == 2) {
		o->testCase = std::atof(argv[1]);
	} else if (argc != 1) {
		std::cerr << "[ERROR] Wrong number of arguments. Expecting 4 or 5 args: X Y Xlength Ylength (output)." << std::endl;
		std::cerr << "(X=100, Y=10, Xlength=0.01, Ylength=0.01 output=current directory with no args)" << std::endl;
		exit(1);
	}

	auto nm = CartesianMesh2DGenerator::generate(o->X_EDGE_ELEMS, o->Y_EDGE_ELEMS, o->X_EDGE_LENGTH, o->Y_EDGE_LENGTH);
	auto c = new EucclhydRemap(o, nm, output);
	c->simulate();
	delete c;
	delete nm;
	delete o;
	Kokkos::finalize();
	return 0;
}
