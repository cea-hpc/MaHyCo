#ifndef EUCCLHYDREMAP_H
#define EUCCLHYDREMAP_H

/*---------------------------------------*/
/*---------------------------------------*/

#include <stddef.h>                       // for size_t
#include <Kokkos_Core.hpp>                // for KOKKOS_LAMBDA
#include <OpenMP/Kokkos_OpenMP_Exec.hpp>  // for OpenMP::impl_is_initialized
#include <algorithm>                      // for copy
#include <array>                          // for array
#include <string>                         // for allocator, string
#include <vector>                         // for vector
#include "EucclhydRemap.h"
#include "mesh/CartesianMesh2D.h"  // for CartesianMesh2D, CartesianM...
#include "mesh/MeshGeometry.h"     // for MeshGeometry
#include "mesh/PvdFileWriter2D.h"  // for PvdFileWriter2D
#include "types/Types.h"           // for RealArray1D, RealArray2D
#include "utils/Timer.h"           // for Timer

/*---------------------------------------*/
/*---------------------------------------*/
using namespace nablalib;

class EucclhydRemap {
 public:
  static const int dim = 2;
  static const int nbmatmax = 3;
  static const int nbequamax =
      3 * nbmatmax + 2 + 1;  // (volumes, masses, energies internes) * nbmatmax
                             // + vitesses + energie cinétique
  struct Options {
    // Should be const but usefull to set them from main args
    RealArray1D<dim> ex = {{1.0, 0.0}};
    RealArray1D<dim> ey = {{0.0, 1.0}};
    RealArray1D<dim> zeroVect = {{0.0, 0.0}};
    RealArray2D<dim, dim> zeroMat = {{{0.0, 0.0}, {0.0, 0.0}}};
    RealArray1D<nbmatmax> zeroVectmat = {{0.0, 0.0, 0.0}};
    RealArray1D<nbequamax> Uzero = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
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
    double cfl = 0.45;
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

    int leftFluxBC = 0;
    RealArray1D<nbequamax> leftFluxBCValue = Uzero;
    int rightFluxBC = 0;
    RealArray1D<nbequamax> rightFluxBCValue = Uzero;
    int bottomFluxBC = 0;
    RealArray1D<nbequamax> bottomFluxBCValue = Uzero;
    int topFluxBC = 0;
    RealArray1D<nbequamax> topFluxBCValue = Uzero;
    int FluxBC = leftFluxBC + rightFluxBC + bottomFluxBC + topFluxBC;

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
  int nbPart = 0;
  int nbNodes, nbCells, nbFaces, nbCellsOfNode, nbFacesOfCell, nbNodesOfCell,
      nbInnerNodes, nbNodesOfFace, nbInnerVerticalFaces, nbInnerHorizontalFaces,
      nbNeighbourCells, nbCellsOfFace, nbLeftNodes, nbRightNodes, nbTopNodes,
      nbBottomNodes, nbTopLeftNode, nbTopRightNode, nbBottomLeftNode,
      nbBottomRightNode;

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
  int dbgcell1 = -389;
  int dbgcell2 = -126;
  int dbgcell3 = -156;
  int face_debug1 = -2033;
  int face_debug2 = -410;
  int test_debug = 1;

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
  Kokkos::View<RealArray2D<dim, dim>*> gradV;
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
  Kokkos::View<RealArray2D<dim, dim>**> M;
  Kokkos::View<RealArray2D<dim, dim>**> M1;
  Kokkos::View<RealArray2D<dim, dim>**> M2;
  Kokkos::View<RealArray2D<dim, dim>**> M3;
  Kokkos::View<RealArray2D<dim, dim>*> Mnode;
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
  // const size_t maxHardThread =
  // Kokkos::DefaultExecutionSpace::max_hardware_threads();

 public:
  EucclhydRemap(Options* aOptions, CartesianMesh2D* aCartesianMesh2D,
                string output)
      : options(aOptions),
        mesh(aCartesianMesh2D),
        writer("EucclhydRemap", output),
        writerpart("Particules", output),
        nbNodes(mesh->getNbNodes()),
        nbPartMax(1),
        nbCells(mesh->getNbCells()),
        nbFaces(mesh->getNbFaces()),
        nbCellsOfNode(CartesianMesh2D::MaxNbCellsOfNode),
        nbFacesOfCell(CartesianMesh2D::MaxNbFacesOfCell),
        nbNodesOfCell(CartesianMesh2D::MaxNbNodesOfCell),
        nbInnerNodes(mesh->getNbInnerNodes()),
        nbNodesOfFace(CartesianMesh2D::MaxNbNodesOfFace),
        nbInnerVerticalFaces(mesh->getNbInnerVerticalFaces()),
        nbInnerHorizontalFaces(mesh->getNbInnerHorizontalFaces()),
        nbNeighbourCells(CartesianMesh2D::MaxNbNeighbourCells),
        nbCellsOfFace(CartesianMesh2D::MaxNbCellsOfFace),
        nbLeftNodes(mesh->getNbLeftNodes()),
        nbRightNodes(mesh->getNbRightNodes()),
        nbTopNodes(mesh->getNbTopNodes()),
        nbBottomNodes(mesh->getNbBottomNodes()),
        nbTopLeftNode(mesh->getNbTopLeftNode()),
        nbTopRightNode(mesh->getNbTopRightNode()),
        nbBottomLeftNode(mesh->getNbBottomLeftNode()),
        nbBottomRightNode(mesh->getNbBottomRightNode()),
        x_then_y_n(true),
        x_then_y_nplus1(true),
        t_n(0.0),
        t_nplus1(0.0),
        deltat_n(options->deltat_init),
        deltat_nplus1(options->deltat_init),
        nbCalls(0),
        lastDump(0.0),
        vpart("VolumePart", nbPartMax),
        wpart("WeightPart", nbPartMax),
        mpart("MassePart", nbPartMax),
        rpart("RayonPart", nbPartMax),
        rhopart("RhoPart", nbPartMax),
        Repart("RePart", nbPartMax),
        Cdpart("Cdpart", nbPartMax),
        Mcpart("Cdpart", nbPartMax),
        ForceGradp("ForceGradp", nbCells),
        listpart("listepart", nbCells),
        fracpart("fracPart", nbCells),
        Xpart_n0("Xpart_n0", nbPartMax),
        Xpart_n("Xpart_n", nbPartMax),
        Xpart_nplus1("Xpart_nplus1", nbPartMax),
        Vpart_n0("Vpart_n0", nbPartMax),
        Vpart_n("Vpart_n", nbPartMax),
        Vpart_nplus1("Vpart_nplus1", nbPartMax),
        ICellp("Icellp", nbPartMax),
        IMatp("Imatp", nbPartMax),
        X("X", nbNodes),
        Xc("Xc", nbCells),
        Xf("Xf", nbFaces),
        XfLagrange("Xf", nbFaces),
        LfLagrange("LfLagrange", nbFaces),
        XLagrange("XLagrange", nbNodes),
        XcLagrange("XcLagrange", nbCells),
        HvLagrange("HvLagrange", nbCells),
        Xc_x("Xc_x", nbCells),
        Xc_y("Xc_y", nbCells),
        lpc_n("lpc_n", nbNodes, nbCellsOfNode),
        nplus("nplus", nbNodes, nbCellsOfNode),
        nminus("nminus", nbNodes, nbCellsOfNode),
        lplus("lplus", nbNodes, nbCellsOfNode),
        lminus("lminus", nbNodes, nbCellsOfNode),
        p("p", nbCells),
        pp("pp", nbCells),
        m("m", nbCells),
        mp("mp", nbCells),
        v("v", nbCells),
        fracmass("fracmass", nbCells),
        mixte("mixte", nbCells),
        pure("pure", nbCells),
        fracvol("fracvol", nbCells),
        fracvolnode("fracvolnode", nbNodes),
        fracvol1("fracvol1", nbCells),
        fracvol2("fracvol2", nbCells),
        fracvol3("fracvol3", nbCells),
        vLagrange("vLagrange", nbCells),
        vpLagrange("vpLagrange", nbCells),
        perim("perim", nbCells),
        vitson("vitson", nbCells),
        vitsonp("vitsonp", nbCells),
        deltaxLagrange("deltaxLagrange", nbFaces),
        faceLength("faceLength", nbFaces),
        faceLengthLagrange("faceLengthLagrange", nbFaces),
        faceNormalVelocity("faceNormalVelocity", nbFaces),
        faceNormal("faceNormal", nbFaces),
        outerFaceNormal("outerFaceNormal", nbCells, nbFacesOfCell),
        deltatc("deltatc", nbCells),
        Vnode_n("Vnode_n", nbNodes),
        Vnode_nplus1("Vnode_nplus1", nbNodes),
        Vnode_n0("Vnode_n0", nbNodes),
        ETOT_0("ETOT_0", nbCells),
        ETOT_T("ETOT_T", nbCells),
        ETOT_L("ETOT_L", nbCells),
        MTOT_0("MTOT_0", nbCells),
        MTOT_T("MTOT_T", nbCells),
        MTOT_L("MTOT_L", nbCells),
        rho_n("rho_n", nbCells),
        rho_nplus1("rho_nplus1", nbCells),
        rho_n0("rho_n0", nbCells),
        rhop_n("rhop_n", nbCells),
        rhop_nplus1("rhop_nplus1", nbCells),
        rhop_n0("rhop_n0", nbCells),
        V_n("V_n", nbCells),
        V_nplus1("V_nplus1", nbCells),
        V_n0("V_n0", nbCells),
        Vxc("Vxc", nbCells),
        Vyc("Vyc", nbCells),
        eps_n("eps_n", nbCells),
        eps_nplus1("eps_nplus1", nbCells),
        eps_n0("eps_n0", nbCells),
        epsp_n("epsp_n", nbCells),
        epsp_nplus1("epsp_nplus1", nbCells),
        epsp_n0("epsp_n0", nbCells),
        delta_ec("delta_ec", nbCells),
        ULagrange("ULagrange", nbCells),
        Uremap1("Uremap1", nbCells),
        Uremap2("Uremap2", nbCells),
        gradPhiFace1("gradPhiFace1", nbFaces),
        gradPhiFace2("gradPhiFace2", nbFaces),
        gradPhi1("gradPhi1", nbCells),
        gradPhi2("gradPhi2", nbCells),
        phiFace1("phiFace1", nbFaces),
        phiFace2("phiFace2", nbFaces),
        deltaPhiFaceAv("deltaPhiFaceAv", nbCells),
        deltaPhiFaceAr("deltaPhiFaceAr", nbCells),
        Phi("Phi", nbCells),
        p_extrap("p_extrap", nbCells, nbNodesOfCell),
        pp_extrap("pp_extrap", nbCells, nbNodesOfCell),
        V_extrap("V_extrap", nbCells, nbNodesOfCell),
        gradp("gradp", nbCells),
        gradp1("gradp1", nbCells),
        gradp2("gradp2", nbCells),
        gradp3("gradp3", nbCells),
        gradf1("gradf1", nbCells),
        gradf2("gradf2", nbCells),
        gradf3("gradf3", nbCells),
        gradV("gradV", nbCells),
        F_n("F_n", nbNodes, nbCellsOfNode),
        F_nplus1("F_nplus1", nbNodes, nbCellsOfNode),
        F_n0("F_n0", nbNodes, nbCellsOfNode),
        F1_n("F1_n", nbNodes, nbCellsOfNode),
        F1_nplus1("F1_nplus1", nbNodes, nbCellsOfNode),
        F2_n("F2_n", nbNodes, nbCellsOfNode),
        F2_nplus1("F2_nplus1", nbNodes, nbCellsOfNode),
        F3_n("F3_n", nbNodes, nbCellsOfNode),
        F3_nplus1("F3_nplus1", nbNodes, nbCellsOfNode),
        G("G", nbNodes),
        M("M", nbNodes, nbCellsOfNode),
        M1("M1", nbNodes, nbCellsOfNode),
        M2("M2", nbNodes, nbCellsOfNode),
        M3("M3", nbNodes, nbCellsOfNode),
        Mnode("Mnode", nbNodes) {
    // Copy node coordinates
    const auto& gNodes = mesh->getGeometry()->getNodes();
    Kokkos::parallel_for(nbNodes, KOKKOS_LAMBDA(const int& rNodes) {
      X(rNodes) = gNodes[rNodes];
    });
  }

 private:
  void computeBoundaryNodeVelocities() noexcept;
  RealArray1D<dim> nodeVelocityBoundaryCondition(int BC,
                                                 RealArray1D<dim> BCValue,
                                                 RealArray2D<dim, dim> Mp,
                                                 RealArray1D<dim> Gp);
  RealArray1D<dim> nodeVelocityBoundaryConditionCorner(
      int BC1, RealArray1D<dim> BCValue1, int BC2, RealArray1D<dim> BCValue2,
      RealArray2D<dim, dim> Mp, RealArray1D<dim> Gp);
  RealArray1D<nbequamax> computeBoundaryFluxes(int proj, int cCells,
                                               RealArray1D<dim> exy);

  void initBoundaryConditions() noexcept;
  void initMeshGeometryForCells() noexcept;
  void initVpAndFpc() noexcept;
  void initCellInternalEnergy() noexcept;
  void initCellVelocity() noexcept;
  void initDensity() noexcept;
  void initMeshGeometryForFaces() noexcept;
  void initPart() noexcept;
  void setUpTimeLoopN() noexcept;

  void computeCornerNormal() noexcept;
  void computeEOS() noexcept;
  void computeGradients() noexcept;
  void computeMass() noexcept;
  void computeDissipationMatrix() noexcept;
  void computedeltatc() noexcept;
  void extrapolateValue() noexcept;
  void computeG() noexcept;
  void computeNodeDissipationMatrixAndG() noexcept;
  void computeNodeVelocity() noexcept;
  void computeFaceVelocity() noexcept;
  void computeLagrangePosition() noexcept;
  void computeSubCellForce() noexcept;
  void computeLagrangeVolumeAndCenterOfGravity() noexcept;
  void computeFacedeltaxLagrange() noexcept;
  void updateCellCenteredLagrangeVariables() noexcept;

  void computeGradPhiFace1() noexcept;
  void computeGradPhi1() noexcept;
  void computeUpwindFaceQuantitiesForProjection1() noexcept;
  void computeUremap1() noexcept;

  void computeGradPhiFace2() noexcept;
  void computeGradPhi2() noexcept;
  void computeUpwindFaceQuantitiesForProjection2() noexcept;
  void computeUremap2() noexcept;

  void remapCellcenteredVariable() noexcept;

  void updateParticlePosition() noexcept;
  void updateParticleCoefficients() noexcept;
  void updateParticleVelocity() noexcept;
  void updateParticleRetroaction() noexcept;
  void switchalpharho_rho() noexcept;
  void switchrho_alpharho() noexcept;

  RealArray2D<2, 2> inverse(RealArray2D<2, 2> a);
  double divideNoExcept(double a, double b);
  template <size_t N, size_t M>
  RealArray2D<N, M> tensProduct(RealArray1D<N> a, RealArray1D<M> b);
  double crossProduct2d(RealArray1D<2> a, RealArray1D<2> b);

  double fluxLimiter(int projectionLimiterId, double r);
  double fluxLimiterPP(int projectionLimiterId, double gradplus,
                       double gradmoins, double y0, double yplus, double ymoins,
                       double h0, double hplus, double hmoins);
  double computeY0(int projectionLimiterId, double y0, double yplus,
                   double ymoins, double h0, double hplus, double hmoins,
                   int type);
  double computexgxd(double y0, double yplus, double ymoins, double h0,
                     double y0plus, double y0moins, int type);
  double computeygyd(double y0, double yplus, double ymoins, double h0,
                     double y0plus, double y0moins, double grady, int type);
  double INTY(double X, double x0, double y0, double x1, double y1);
  double INT2Y(double X, double x0, double y0, double x1, double y1);
  template <size_t d>
  RealArray1D<d> computeAndLimitGradPhi(
      int projectionLimiterId, RealArray1D<d> gradphiplus,
      RealArray1D<d> gradphimoins, RealArray1D<d> phi, RealArray1D<d> phiplus,
      RealArray1D<d> phimoins, double h0, double hplus, double hmoins);
  template <size_t d>
  RealArray1D<d> computeIntersectionPP(
      RealArray1D<d> gradphi, RealArray1D<d> phi, RealArray1D<d> phiplus,
      RealArray1D<d> phimoins, double h0, double hplus, double hmoins,
      double face_normal_velocity, double deltat_n, int type, int cell,
      double flux_threhold);
  template <size_t d>
  RealArray1D<d> computeIntersectionPPPure(
      RealArray1D<d> gradphi, RealArray1D<d> phi, RealArray1D<d> phiplus,
      RealArray1D<d> phimoins, double h0, double hplus, double hmoins,
      double face_normal_velocity, double deltat_n, int type, int cell,
      double flux_threhold);
  template <size_t d>
  RealArray1D<d> computeUpwindFaceQuantities(
      RealArray1D<dim> face_normal, double face_normal_velocity, double delta_x,
      RealArray1D<dim> x_f, RealArray1D<d> phi_cb, RealArray1D<d> grad_phi_cb,
      RealArray1D<dim> x_cb, RealArray1D<d> phi_cf, RealArray1D<d> grad_phi_cf,
      RealArray1D<dim> x_cf);
  RealArray1D<dim> xThenYToDirection(bool x_then_y_);
  template <size_t d>
  RealArray1D<d> computeRemapFlux(int projectionAvecPlateauPente,
                                  double face_normal_velocity,
                                  RealArray1D<dim> face_normal,
                                  double face_length, RealArray1D<d> phi_face,
                                  RealArray1D<dim> outer_face_normal,
                                  RealArray1D<dim> exy, double deltat_n);

  /**
   * Job dumpVariables called @2.0 in executeTimeLoopN method.
   * In variables: Xc_x, Xc_y, eps_n, m, p, rho_n, t_n, v
   * Out variables:
   */
  void dumpVariables() noexcept;

  /**
   * Job executeTimeLoopN called @4.0 in simulate method.
   * In variables: F_n, F_nplus1, G, M, Mnode, ULagrange, Uremap1, Uremap2,
   * V_extrap, V_n, Vnode_n, Vnode_nplus1, X, XLagrange, Xc, XcLagrange, Xc_x,
   * Xc_y, Xf, bottomBC, bottomBCValue, c, cfl, deltat_n, deltat_nplus1,
   * deltatc, deltaxLagrange, eos, eosPerfectGas, eps_n, faceLength, faceNormal,
   * faceNormalVelocity, gamma, gradPhi1, gradPhi2, gradPhiFace1, gradPhiFace2,
   * gradV, gradp, leftBC, leftBCValue, lminus, lpc_n, lplus, m, nminus, nplus,
   * outerFaceNormal, p, p_extrap, perim, phiFace1, phiFace2,
   * projectionLimiterId, projectionOrder, rho_n, rightBC, rightBCValue,
   * spaceOrder, t_n, topBC, topBCValue, v, vLagrange, x_then_y_n Out variables:
   * F_nplus1, G, M, Mnode, ULagrange, Uremap1, Uremap2, V_extrap, V_nplus1,
   * Vnode_nplus1, XLagrange, XcLagrange, c, deltat_nplus1, deltatc,
   * deltaxLagrange, eps_nplus1, faceNormalVelocity, gradPhi1, gradPhi2,
   * gradPhiFace1, gradPhiFace2, gradV, gradp, m, p, p_extrap, phiFace1,
   * phiFace2, rho_nplus1, t_nplus1, vLagrange, x_then_y_nplus1
   */
  void executeTimeLoopN() noexcept;

  /**
   * Job computedeltat called @3.0 in executeTimeLoopN method.
   * In variables: cfl, deltat_n, deltatc
   * Out variables: deltat_nplus1
   */
  void computedeltat() noexcept;

  /**
   * Job updateTime called @4.0 in executeTimeLoopN method.
   * In variables: deltat_nplus1, t_n
   * Out variables: t_nplus1
   */
  void updateTime() noexcept;

 public:
  void simulate();
};

#include "Utiles-Impl.h"
#include "UtilesRemap-Impl.h"

#endif  // EUCCLHYDREMAP_H