#include <math.h>                   // for floor, sqrt
#include <Kokkos_Core.hpp>          // for deep_copy
#include <algorithm>                // for copy
#include <array>                    // for array
#include <iostream>                 // for operator<<, basic_ostream::ope...
#include <vector>                   // for allocator, vector
#include "EucclhydRemap.h"          // for EucclhydRemap, EucclhydRemap::...
#include "mesh/CartesianMesh2D.h"   // for CartesianMesh2D
#include "types/ArrayOperations.h"  // for minus, multiply, plus, divide
#include "types/MathFunctions.h"    // for max, norm, dot
#include "types/MultiArray.h"       // for operator<<
#include "utils/Utils.h"            // for indexOf

void EucclhydRemap::initBoundaryConditions() noexcept {
  if (options->testCase == options->SodCase ||
      options->testCase == options->BiSodCase) {
    // maillage 200 5 0.005 0.02
    options->leftBC = options->symmetry;
    options->leftBCValue = options->ey;

    options->rightBC = options->symmetry;
    options->rightBCValue = options->ey;

    options->topBC = options->symmetry;
    options->topBCValue = options->ex;

    options->bottomBC = options->symmetry;
    options->bottomBCValue = options->ex;
  } else if (options->testCase == options->BiShockBubble) {
    // maillage 520 64 0.00125 0.00125
    options->leftBC = options->symmetry;
    options->leftBCValue = options->ey;

    options->rightBC = options->imposedVelocity;
    options->rightBCValue = {{-124.824, 0.0}};
    options->rightFluxBC = 1;
    options->rightFluxBCValue = {
        {1.0, 0.0, 0.0, 1.0, 0.0, 0.0, -124.824, 0.0, 250000, 0.0, 0.0, 0.0}};

    options->topBC = options->symmetry;
    options->topBCValue = options->ex;

    options->bottomBC = options->symmetry;
    options->bottomBCValue = options->ex;

  } else if (options->testCase == options->SedovTestCase ||
             options->testCase == options->BiSedovTestCase) {
    // const ℕ leftBC = symmetry; const ℝ[2] leftBCValue = ey;
    // const ℕ rightBC = imposedVelocity; const ℝ[2] rightBCValue = zeroVect;
    // const ℕ topBC = imposedVelocity; const ℝ[2] topBCValue = zeroVect;
    // const ℕ bottomBC = symmetry; const ℝ[2] bottomBCValue = ex;

    options->leftBC = options->symmetry;
    options->leftBCValue = options->ey;

    options->rightBC = options->imposedVelocity;
    options->rightBCValue = options->zeroVect;

    options->topBC = options->imposedVelocity;
    options->topBCValue = options->zeroVect;

    options->bottomBC = options->symmetry;
    options->bottomBCValue = options->ex;
  } else if (options->testCase == options->TriplePoint ||
             options->testCase == options->BiTriplePoint) {
    // maillage 140 60 0.0005 0.0005
    options->leftBC = options->symmetry;
    options->leftBCValue = options->ey;

    options->rightBC = options->symmetry;
    options->rightBCValue = options->ey;

    options->topBC = options->symmetry;
    options->topBCValue = options->ex;

    options->bottomBC = options->symmetry;
    options->bottomBCValue = options->ex;
  } else if (options->testCase == options->NohTestCase ||
             options->testCase == options->BiNohTestCase) {
    // const ℕ leftBC = symmetry; const ℝ[2] leftBCValue = ey;
    // const ℕ rightBC = imposedVelocity; const ℝ[2] rightBCValue = zeroVect;
    // const ℕ topBC = imposedVelocity; const ℝ[2] topBCValue = zeroVect;
    // const ℕ bottomBC = symmetry; const ℝ[2] bottomBCValue = ex;

    options->leftBC = options->symmetry;
    options->leftBCValue = options->ey;

    options->rightBC = options->imposedVelocity;
    options->rightBCValue = options->zeroVect;

    options->topBC = options->imposedVelocity;
    options->topBCValue = options->zeroVect;

    options->bottomBC = options->symmetry;
    options->bottomBCValue = options->ex;
  }
  options->FluxBC = options->leftFluxBC + options->rightFluxBC +
                    options->bottomFluxBC + options->topFluxBC;
}
/**
 * Job initMeshGeometryForCells called @1.0 in simulate method.
 * In variables: X
 * Out variables: Xc, Xc_x, Xc_y, perim, v
 */
void EucclhydRemap::initMeshGeometryForCells() noexcept {
  Kokkos::parallel_for(
      "initMeshGeometryForCells", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        int cId(cCells);
        double reduction11 = 0.0;
        {
          auto nodesOfCellC(mesh->getNodesOfCell(cId));
          for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
               pNodesOfCellC++) {
            int pId(nodesOfCellC[pNodesOfCellC]);
            int pPlus1Id(nodesOfCellC[(pNodesOfCellC + 1 + nbNodesOfCell) %
                                      nbNodesOfCell]);
            int pNodes(pId);
            int pPlus1Nodes(pPlus1Id);
            reduction11 =
                reduction11 + (crossProduct2d(X(pNodes), X(pPlus1Nodes)));
          }
        }
        double vol = 0.5 * reduction11;
        RealArray1D<dim> reduction12 = options->zeroVect;
        {
          auto nodesOfCellC(mesh->getNodesOfCell(cId));
          for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
               pNodesOfCellC++) {
            int pId(nodesOfCellC[pNodesOfCellC]);
            int pPlus1Id(nodesOfCellC[(pNodesOfCellC + 1 + nbNodesOfCell) %
                                      nbNodesOfCell]);
            int pNodes(pId);
            int pPlus1Nodes(pPlus1Id);
            reduction12 = ArrayOperations::plus(
                reduction12,
                (ArrayOperations::multiply(
                    crossProduct2d(X(pNodes), X(pPlus1Nodes)),
                    (ArrayOperations::plus(X(pNodes), X(pPlus1Nodes))))));
          }
        }
        RealArray1D<dim> xc =
            ArrayOperations::multiply(1.0 / (6.0 * vol), reduction12);
        Xc(cCells) = xc;
        Xc_x(cCells) = xc[0];
        Xc_y(cCells) = xc[1];
        v(cCells) = vol;
        double reduction13 = 0.0;
        {
          auto nodesOfCellC(mesh->getNodesOfCell(cId));
          for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
               pNodesOfCellC++) {
            int pId(nodesOfCellC[pNodesOfCellC]);
            int pPlus1Id(nodesOfCellC[(pNodesOfCellC + 1 + nbNodesOfCell) %
                                      nbNodesOfCell]);
            int pNodes(pId);
            int pPlus1Nodes(pPlus1Id);
            reduction13 =
                reduction13 + (MathFunctions::norm(ArrayOperations::minus(
                                  X(pNodes), X(pPlus1Nodes))));
          }
        }
        perim(cCells) = reduction13;
      });
}

/**
 * Job initVpAndFpc called @1.0 in simulate method.
 * In variables: zeroVect
 * Out variables: F_n0, Vnode_n0
 */
void EucclhydRemap::initVpAndFpc() noexcept {
  Kokkos::parallel_for(
      "initVpAndFpc", nbNodes, KOKKOS_LAMBDA(const int& pNodes) {
        int pId(pNodes);
        {
          auto cellsOfNodeP(mesh->getCellsOfNode(pId));
          for (int cCellsOfNodeP = 0; cCellsOfNodeP < cellsOfNodeP.size();
               cCellsOfNodeP++) {
            Vnode_n0(pNodes) = options->zeroVect;
            F_n0(pNodes, cCellsOfNodeP) = options->zeroVect;
          }
        }
      });
}

/**
 * Job initCellInternalEnergy called @2.0 in simulate method.
 * In variables: NohTestCase, SedovTestCase, SodCase, TriplePoint, X,
 * X_EDGE_LENGTH, Xc, Y_EDGE_LENGTH, gamma, p0, rho0, testCase, threshold Out
 * variables: eps_n0
 */
void EucclhydRemap::initCellInternalEnergy() noexcept {
  Kokkos::parallel_for("initCellInternalEnergy", nbCells,
                       KOKKOS_LAMBDA(const int& cCells) {
                         for (imat = 0; imat < nbmatmax; imat++) {
                           epsp_n0(cCells)[imat] = 0.0;
                         }
                       });

  if (options->testCase == options->NohTestCase ||
      options->testCase == options->BiNohTestCase ||
      options->testCase == options->UnitTestCase ||
      options->testCase == options->BiUnitTestCase) {
    double eps0 = options->p0 / ((options->gamma - 1.0) * options->rho0);
    Kokkos::parallel_for("initCellInternalEnergy", nbCells,
                         KOKKOS_LAMBDA(const int& cCells) {
                           eps_n0(cCells) = eps0;
                           for (imat = 0; imat < nbmatmax; imat++)
                             epsp_n0(cCells)[imat] = eps0;
                         });
  } else if (options->testCase == options->SodCase) {
    Kokkos::parallel_for("initCellInternalEnergy", nbCells,
                         KOKKOS_LAMBDA(const int& cCells) {
                           double pInit;
                           double rhoInit;
                           double epsInit;
                           if (Xc(cCells)[0] <= 0.5) {
                             pInit = 1.0;
                             rhoInit = 1.0;
                           } else {
                             pInit = 0.1;
                             rhoInit = 0.125;
                           }
                           epsInit = pInit / ((options->gamma - 1.0) * rhoInit);
                           eps_n0(cCells) = epsInit;
                           epsp_n0(cCells)[0] = epsInit;
                         });
  } else if (options->testCase == options->BiSodCase) {
    Kokkos::parallel_for(
        "initCellInternalEnergy", nbCells, KOKKOS_LAMBDA(const int& cCells) {
          double pInit;
          double rhoInit;
          double epsInit;
          double r = 0.;
          // r= sqrt(Xc(cCells)[0]*Xc(cCells)[0]+Xc(cCells)[1]*Xc(cCells)[1]);
          r = Xc(cCells)[0];  // sod en X
          // r= Xc(cCells)[1]; //
          if (r <= 0.5) {
            pInit = 1.0;
            rhoInit = 1.0;
            epsInit = pInit / ((options->gammap[0] - 1.0) * rhoInit);
            epsp_n0(cCells)[0] = epsInit;
            epsp_n0(cCells)[1] = 0.;
            eps_n0(cCells) = epsInit;
            // std::cout << " cell " << cCells << "  e1= " << epsp_n0(cCells)[0]
            // << "  e2= " << epsp_n0(cCells)[0] << std::endl;
          } else {
            pInit = 0.1;
            rhoInit = 0.125;
            epsInit = pInit / ((options->gammap[1] - 1.0) * rhoInit);
            epsp_n0(cCells)[0] = 0.;
            epsp_n0(cCells)[1] = epsInit;
            eps_n0(cCells) = epsInit;
            // std::cout << " cell " << cCells << "  e1= " << epsp_n0(cCells)[0]
            // << "  e2= " << epsp_n0(cCells)[0] << std::endl;
          }
        });
  } else if (options->testCase == options->BiShockBubble) {
    Kokkos::parallel_for(
        "initCellInternalEnergy", nbCells, KOKKOS_LAMBDA(const int& cCells) {
          double pInit;
          double rhoInit;
          double epsInit;
          // Air partout
          rhoInit = 1.0;
          pInit = 1.e5;
          epsInit = pInit / ((options->gammap[0] - 1.0) * rhoInit);
          epsp_n0(cCells)[0] = epsInit;
          epsp_n0(cCells)[1] = 0.;
          eps_n0(cCells) = epsInit;
          // bulle surchargera l'aire
          // centre de la bulle
          RealArray1D<dim> Xb = {{0.320, 0.04}};
          double rb = 0.025;
          double r = sqrt((Xc(cCells)[0] - Xb[0]) * (Xc(cCells)[0] - Xb[0]) +
                          (Xc(cCells)[1] - Xb[1]) * (Xc(cCells)[1] - Xb[1]));

          if (r <= rb) {
            pInit = 1.e5;
            rhoInit = 0.182;
            epsInit = pInit / ((options->gammap[0] - 1.0) * rhoInit);
            epsp_n0(cCells)[0] = 0.;
            epsp_n0(cCells)[1] = epsInit;
            eps_n0(cCells) = epsInit;
          }
        });
  } else if (options->testCase == options->SedovTestCase ||
             options->testCase == options->BiSedovTestCase) {
    double eps1 = options->p0 / ((options->gamma - 1.0) * options->rho0);
    Kokkos::parallel_for(
        "initCellInternalEnergy", nbCells, KOKKOS_LAMBDA(const int& cCells) {
          int cId(cCells);
          bool isCenterCell = false;
          {
            auto nodesOfCellC(mesh->getNodesOfCell(cId));
            for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
                 pNodesOfCellC++) {
              int pId(nodesOfCellC[pNodesOfCellC]);
              int pNodes(pId);
              if (MathFunctions::norm(X(pNodes)) < options->threshold)
                isCenterCell = true;
            }
          }
          if (isCenterCell) {
            double total_energy_deposit = 0.979264;
            double dx = options->X_EDGE_LENGTH;
            double dy = options->Y_EDGE_LENGTH;
            eps_n0(cCells) = (eps1 + total_energy_deposit / (4.0 * dx * dy));
          } else {
            eps_n0(cCells) = eps1;
          }
        });
  } else if (options->testCase == options->TriplePoint) {
    Kokkos::parallel_for("initCellInternalEnergy", nbCells,
                         KOKKOS_LAMBDA(const int& cCells) {
                           double pInit;
                           double rhoInit;
                           double epsInit;
                           if (Xc(cCells)[0] <= 0.01) {
                             pInit = 1.0;
                             rhoInit = 1.0;
                           } else {
                             if (Xc(cCells)[1] <= 0.015) {
                               pInit = 0.1;
                               rhoInit = 1.;
                             } else {
                               pInit = 0.1;
                               rhoInit = 0.1;
                             }
                           }
                           epsInit = pInit / ((options->gamma - 1.0) * rhoInit);
                           eps_n0(cCells) = epsInit;
                           epsp_n0(cCells)[0] = epsInit;
                         });
  } else if (options->testCase == options->BiTriplePoint) {
    Kokkos::parallel_for(
        "initCellInternalEnergy", nbCells, KOKKOS_LAMBDA(const int& cCells) {
          double pInit;
          double rhoInit;
          double epsInit;
          if (Xc(cCells)[0] <= 0.01) {
            pInit = 1.0;  // 1.e5; // 1.0;
            rhoInit = 1.0;
            epsInit = pInit / ((options->gamma - 1.0) * rhoInit);
            eps_n0(cCells) = epsInit;
            epsp_n0(cCells)[0] = epsInit;
          } else {
            if (Xc(cCells)[1] <= 0.015) {
              pInit = 0.1;  // 1.e4; // 0.1;
              rhoInit = 1.;
              epsInit = pInit / ((options->gamma - 1.0) * rhoInit);
              eps_n0(cCells) = epsInit;
              epsp_n0(cCells)[1] = epsInit;
            } else {
              pInit = 0.1;  // 1.e4; // 0.1;
              rhoInit = 0.1;
              epsInit = pInit / ((options->gamma - 1.0) * rhoInit);
              eps_n0(cCells) = epsInit;
              epsp_n0(cCells)[2] = epsInit;
            }
          }
        });
  }
}

/**
 * Job initCellVelocity called @2.0 in simulate method.
 * In variables: NohTestCase, SedovTestCase, SodCase, Xc, testCase, u0
 * Out variables: V_n0
 */
void EucclhydRemap::initCellVelocity() noexcept {
  Kokkos::parallel_for(
      "initCellVelocity", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        if (options->testCase == options->NohTestCase ||
            options->testCase == options->BiNohTestCase) {
          double n1 = Xc(cCells)[0];
          double n2 = Xc(cCells)[1];
          double normVect = MathFunctions::sqrt(n1 * n1 + n2 * n2);
          V_n0(cCells)[0] = -options->u0 * n1 / normVect;
          V_n0(cCells)[1] = -options->u0 * n2 / normVect;
        } else if (options->testCase == options->SedovTestCase ||
                   options->testCase == options->SodCase ||
                   options->testCase == options->TriplePoint ||
                   options->testCase == options->BiSedovTestCase ||
                   options->testCase == options->BiSodCase ||
                   options->testCase == options->BiTriplePoint) {
          V_n0(cCells)[0] = 0.0;
          V_n0(cCells)[1] = 0.0;
        } else if (options->testCase == options->BiShockBubble) {
          if (Xc(cCells)[0] >= 0.60) {
            V_n0(cCells)[0] = -124.824;
            V_n0(cCells)[1] = 0.0;
          }
        } else if (options->testCase == options->UnitTestCase ||
                   options->testCase == options->BiUnitTestCase) {
          V_n0(cCells)[0] = 1.0;
          V_n0(cCells)[1] = 0.0;
        }
      });
}

/**
 * Job initDensity called @2.0 in simulate method.
 * In variables: SodCase, Xc, rho0, testCase
 * Out variables: rho_n0
 */
void EucclhydRemap::initDensity() noexcept {
  Kokkos::parallel_for("initDensity", nbCells,
                       KOKKOS_LAMBDA(const int& cCells) {
                         for (imat = 0; imat < nbmatmax; imat++) {
                           fracvol(cCells)[imat] = 0.0;
                           fracmass(cCells)[imat] = 0.0;
                           rhop_n0(cCells)[imat] = 0.0;
                         }
                       });
  if (options->testCase == options->UnitTestCase) {
    Kokkos::parallel_for("initDensity", nbCells,
                         KOKKOS_LAMBDA(const int& cCells) {
                           if (Xc(cCells)[0] <= 0.5) {
                             fracvol(cCells)[0] = 1.;
                             fracvol(cCells)[1] = 0.;
                             fracmass(cCells)[0] = 1.;
                             fracmass(cCells)[1] = 0.;
                             rho_n0(cCells) = 1.0;
                             rhop_n0(cCells)[0] = 1.0;
                             rhop_n0(cCells)[1] = 0.;
                           } else {
                             fracvol(cCells)[0] = 1.;
                             fracvol(cCells)[1] = 0.;
                             fracmass(cCells)[0] = 1.;
                             fracmass(cCells)[1] = 0.;
                             rho_n0(cCells) = 1.;
                             rhop_n0(cCells)[0] = 1.;
                             rhop_n0(cCells)[1] = 0.;
                           }
                         });
  } else if (options->testCase == options->BiUnitTestCase) {
    Kokkos::parallel_for("initDensity", nbCells,
                         KOKKOS_LAMBDA(const int& cCells) {
                           if (Xc(cCells)[0] <= 0.5) {
                             fracvol(cCells)[0] = 1.;
                             fracvol(cCells)[1] = 0.;
                             fracmass(cCells)[0] = 1.;
                             fracmass(cCells)[1] = 0.;
                             rho_n0(cCells) = 1.0;
                             rhop_n0(cCells)[0] = 1.0;
                             rhop_n0(cCells)[1] = 0.;
                           } else {
                             fracvol(cCells)[0] = 0.;
                             fracvol(cCells)[1] = 1.;
                             fracmass(cCells)[0] = 0.;
                             fracmass(cCells)[1] = 1.;
                             rho_n0(cCells) = 1.;
                             rhop_n0(cCells)[0] = 0.;
                             rhop_n0(cCells)[1] = 1.0;
                           }
                         });
  } else if (options->testCase == options->SodCase) {
    Kokkos::parallel_for("initDensity", nbCells,
                         KOKKOS_LAMBDA(const int& cCells) {
                           fracvol(cCells)[0] = 1.;
                           fracvol(cCells)[1] = 0.;
                           fracmass(cCells)[0] = 1.;
                           fracmass(cCells)[1] = 0.;
                           if (Xc(cCells)[0] <= 0.5) {
                             rho_n0(cCells) = 1.0;
                             rhop_n0(cCells)[0] = 1.0;
                             rhop_n0(cCells)[1] = 1.0;
                           } else {
                             rho_n0(cCells) = 0.125;
                             rhop_n0(cCells)[0] = 0.125;
                             rhop_n0(cCells)[1] = 0.125;
                           }
                         });
  } else if (options->testCase == options->BiSodCase) {
    Kokkos::parallel_for("initDensity", nbCells,
                         KOKKOS_LAMBDA(const int& cCells) {
                           double r = 0.;
                           // r=
                           // sqrt(Xc(cCells)[0]*Xc(cCells)[0]+Xc(cCells)[1]*Xc(cCells)[1]);
                           r = Xc(cCells)[0];  // sod en X
                           // r= Xc(cCells)[1]; // sod en Y

                           if (r <= 0.5) {
                             fracvol(cCells)[0] = 1.;
                             fracvol(cCells)[1] = 0.;
                             fracmass(cCells)[0] = 1.;
                             fracmass(cCells)[1] = 0.;
                             rho_n0(cCells) = 1.0;
                             rhop_n0(cCells)[0] = 1.0;
                             rhop_n0(cCells)[1] = 0.;
                           } else {
                             fracvol(cCells)[0] = 0.;
                             fracvol(cCells)[1] = 1.;
                             fracmass(cCells)[0] = 0.;
                             fracmass(cCells)[1] = 1.;
                             rho_n0(cCells) = 0.125;
                             rhop_n0(cCells)[0] = 0.;
                             rhop_n0(cCells)[1] = 0.125;
                           }
                         });
  } else if (options->testCase == options->BiShockBubble) {
    Kokkos::parallel_for(
        "initDensity", nbCells, KOKKOS_LAMBDA(const int& cCells) {
          double pInit;
          double rhoInit;
          double epsInit;
          // Air partout
          rhop_n0(cCells)[0] = 1.0;
          rho_n0(cCells) = 1.0;
          fracvol(cCells)[0] = 1.;
          fracvol(cCells)[1] = 0.;

          fracmass(cCells)[0] = 1.;
          fracmass(cCells)[1] = 0.;
          // bulle surchargera l'aire
          // centre de la bulle
          RealArray1D<dim> Xb = {{0.320, 0.04}};
          double rb = 0.025;
          double r = sqrt((Xc(cCells)[0] - Xb[0]) * (Xc(cCells)[0] - Xb[0]) +
                          (Xc(cCells)[1] - Xb[1]) * (Xc(cCells)[1] - Xb[1]));
          if (r <= rb) {
            rhop_n0(cCells)[0] = 0.0;
            rhop_n0(cCells)[1] = 0.182;
            rho_n0(cCells) = 0.182;

            fracvol(cCells)[0] = 0.;
            fracvol(cCells)[1] = 1.;

            fracmass(cCells)[0] = 0.;
            fracmass(cCells)[1] = 1.;
          }
        });
  }

  else if (options->testCase == options->TriplePoint) {
    Kokkos::parallel_for("initDensity", nbCells,
                         KOKKOS_LAMBDA(const int& cCells) {
                           fracvol(cCells)[0] = 1.;
                           fracvol(cCells)[1] = 0.;

                           fracmass(cCells)[0] = 1.;
                           fracmass(cCells)[1] = 0.;
                           if (Xc(cCells)[0] <= 0.01) {
                             // std::cout << " cell " << cCells << "  x= " <<
                             // Xc(cCells)[0] << " y= " << Xc(cCells)[1] <<
                             // std::endl;
                             rho_n0(cCells) = 1.0;
                           } else {
                             if (Xc(cCells)[1] <= 0.015) {
                               // std::cout << " cell cas 2  " << cCells << " x=
                               // " << Xc(cCells)[0] << "  y= " << Xc(cCells)[1]
                               // << std::endl;
                               rho_n0(cCells) = 1.0;
                             } else {
                               // std::cout << " cell cas 3  " << cCells << " x=
                               // " << Xc(cCells)[0] << "  y= " << Xc(cCells)[1]
                               // << std::endl;
                               rho_n0(cCells) = 0.1;
                             }
                           }
                         });
  } else if (options->testCase == options->BiTriplePoint) {
    Kokkos::parallel_for(
        "initDensity", nbCells, KOKKOS_LAMBDA(const int& cCells) {
          if (Xc(cCells)[0] <= 0.01) {
            std::cout << " cell cas 1 " << cCells << "  x= " << Xc(cCells)[0]
                      << "  y= " << Xc(cCells)[1] << std::endl;
            rho_n0(cCells) = 1.0;
            fracvol(cCells)[0] = 1.;
            fracmass(cCells)[0] = 1.;
            rhop_n0(cCells)[0] = 1.0;

          } else {
            if (Xc(cCells)[1] <= 0.015) {
              std::cout << " cell cas 2  " << cCells << "  x= " << Xc(cCells)[0]
                        << "  y= " << Xc(cCells)[1] << std::endl;
              rho_n0(cCells) = 1.0;
              fracvol(cCells)[1] = 1.;
              fracmass(cCells)[1] = 1.;
              rhop_n0(cCells)[1] = 1.0;
            } else {
              std::cout << " cell cas 3  " << cCells << "  x= " << Xc(cCells)[0]
                        << "  y= " << Xc(cCells)[1] << std::endl;
              rho_n0(cCells) = 0.1;
              fracvol(cCells)[2] = 1.;
              fracmass(cCells)[2] = 1.;
              rhop_n0(cCells)[2] = 0.1;
            }
          }
        });
  } else {
    Kokkos::parallel_for("initDensity", nbCells,
                         KOKKOS_LAMBDA(const int& cCells) {
                           fracvol(cCells)[0] = 1.;
                           fracvol(cCells)[1] = 0.;

                           fracmass(cCells)[0] = 1.;
                           fracmass(cCells)[1] = 0.;
                           rho_n0(cCells) = options->rho0;
                         });
  }
  Kokkos::parallel_for("initDensity", nbCells,
                       KOKKOS_LAMBDA(const int& cCells) {
                         // pour les sorties au temps 0:
                         fracvol1(cCells) = fracvol(cCells)[0];
                         fracvol2(cCells) = fracvol(cCells)[1];
                         fracvol3(cCells) = fracvol(cCells)[2];
                         Vxc(cCells) = V_n0(cCells)[0];
                         Vyc(cCells) = V_n0(cCells)[1];
                         // indicateur mailles mixtes
                         int matcell(0);
                         int imatpure(-1);
                         for (int imat = 0; imat < nbmatmax; imat++)
                           if (fracvol(cCells)[imat] > options->threshold) {
                             matcell++;
                             imatpure = imat;
                           }

                         if (matcell > 1) {
                           mixte(cCells) = 1;
                           pure(cCells) = -1;
                         } else {
                           mixte(cCells) = 0;
                           pure(cCells) = imatpure;
                         }
                       });
}

/**
 * Job initMeshGeometryForFaces called @2.0 in simulate method.
 * In variables: X, Xc, ex, ey, threshold
 * Out variables: Xf, faceLength, faceNormal, outerFaceNormal
 */
void EucclhydRemap::initMeshGeometryForFaces() noexcept {
  auto faces(mesh->getFaces());
  Kokkos::parallel_for(
      "initMeshGeometryForFaces", nbFaces, KOKKOS_LAMBDA(const int& fFaces) {
        int fId(faces[fFaces]);
        int n1FirstNodeOfFaceF(mesh->getFirstNodeOfFace(fId));
        int n1Id(n1FirstNodeOfFaceF);
        int n1Nodes(n1Id);
        int n2SecondNodeOfFaceF(mesh->getSecondNodeOfFace(fId));
        int n2Id(n2SecondNodeOfFaceF);
        int n2Nodes(n2Id);
        RealArray1D<dim> X_face = ArrayOperations::multiply(
            0.5, (ArrayOperations::plus(X(n1Nodes), X(n2Nodes))));
        RealArray1D<dim> face_vec =
            ArrayOperations::minus(X(n2Nodes), X(n1Nodes));
        Xf(fFaces) = X_face;
        faceLength(fFaces) = MathFunctions::norm(face_vec);
        {
          auto cellsOfFaceF(mesh->getCellsOfFace(fId));
          for (int cCellsOfFaceF = 0; cCellsOfFaceF < cellsOfFaceF.size();
               cCellsOfFaceF++) {
            int cId(cellsOfFaceF[cCellsOfFaceF]);
            int cCells(cId);
            int fFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), fId));
            outerFaceNormal(cCells, fFacesOfCellC) = ArrayOperations::divide(
                (ArrayOperations::minus(X_face, Xc(cCells))),
                MathFunctions::norm(
                    ArrayOperations::minus(X_face, Xc(cCells))));
            if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                 cCells == dbgcell1)) {
              std::cout << " cell   " << cCells << std::endl;
              std::cout << "fId " << fId
                        << " outerFaceNormal(cCells,fFacesOfCellC) "
                        << outerFaceNormal(cCells, fFacesOfCellC)
                        << " fFaceOfCellC " << fFacesOfCellC << std::endl;
            }
          }
        }
        RealArray1D<dim> face_normal;
        if (MathFunctions::fabs(MathFunctions::dot(face_vec, options->ex)) <
            options->threshold)
          face_normal = options->ex;
        else
          face_normal = options->ey;
        faceNormal(fFaces) = face_normal;
      });
}

void EucclhydRemap::initPart() noexcept {
  Kokkos::parallel_for("initPart", nbPart, KOKKOS_LAMBDA(const int& ipart) {
    Vpart_n0(ipart) = options->zeroVect;

    rpart(ipart) = 5.e-6;
    rhopart(ipart) = 1110.;
    wpart(ipart) = 1.e7;
    vpart(ipart) = 4. * Pi * rpart(ipart) * rpart(ipart) * rpart(ipart) / 3.;
    mpart(ipart) = rhopart(ipart) * vpart(ipart);
  });

  if (options->testCase == options->BiSodCase ||
      options->testCase == options->SodCase)
    Kokkos::parallel_for("initPart", nbPart, KOKKOS_LAMBDA(const int& ipart) {
      Xpart_n0(ipart)[1] = (ipart % 10) * 0.01 + ipart * 0.0001;

      if (ipart < nbPart / 6)
        Xpart_n0(ipart)[0] = 0.545;
      else if (ipart >= nbPart / 6 && ipart < 2 * nbPart / 6)
        Xpart_n0(ipart)[0] = 0.595;
      else if (ipart >= 2 * nbPart / 6 && ipart < 3 * nbPart / 6)
        Xpart_n0(ipart)[0] = 0.645;
      else if (ipart >= 3 * nbPart / 6 && ipart < 4 * nbPart / 6)
        Xpart_n0(ipart)[0] = 0.695;
      else if (ipart >= 4 * nbPart / 6 && ipart < 5 * nbPart / 6)
        Xpart_n0(ipart)[0] = 0.745;
      else if (ipart >= 5 * nbPart / 6)
        Xpart_n0(ipart)[0] = 0.795;

      int icell = MathFunctions::max(
          floor(Xpart_n0(ipart)[0] / options->X_EDGE_LENGTH), 0);
      int jcell = MathFunctions::max(
          floor(Xpart_n0(ipart)[1] / options->Y_EDGE_LENGTH), 0);
      ICellp(ipart) = jcell * options->X_EDGE_ELEMS + icell;
      listpart(ICellp(ipart)).push_back(ipart);
      IMatp(ipart) = 1;
      std::cout << " Part  " << ipart << "  xp= " << Xpart_n0(ipart)[0]
                << "  yp= " << Xpart_n0(ipart)[1] << "  " << ICellp(ipart)
                << std::endl;
    });
  else if (options->testCase == options->BiTriplePoint)
    Kokkos::parallel_for("initPart", nbPart, KOKKOS_LAMBDA(const int& ipart) {
      Xpart_n0(ipart)[1] = (ipart % 20) * 0.00075 + 0.00060;
      int section = nbPart / 20;
      for (int i = 0; i < section; i++)
        if (ipart >= i * nbPart / section && ipart < (i + 1) * nbPart / section)
          Xpart_n0(ipart)[0] = 0.01 + 0.0005 * (i + 1);

      int icell = MathFunctions::max(
          floor(Xpart_n0(ipart)[0] / options->X_EDGE_LENGTH), 0);
      int jcell = MathFunctions::max(
          floor(Xpart_n0(ipart)[1] / options->Y_EDGE_LENGTH), 0);
      ICellp(ipart) = jcell * options->X_EDGE_ELEMS + icell;
      listpart(ICellp(ipart)).push_back(ipart);
      IMatp(ipart) = 1;
      std::cout << " Part  " << ipart << "  xp= " << Xpart_n0(ipart)[0]
                << "  yp= " << Xpart_n0(ipart)[1] << "  " << ICellp(ipart)
                << std::endl;
    });
}
/**
 * Job setUpTimeLoopN called @3.0 in simulate method.
 * In variables: F_n0, V_n0, Vnode_n0, eps_n0, rho_n0
 * Out variables: F_n, V_n, Vnode_n, eps_n, rho_n
 */
void EucclhydRemap::setUpTimeLoopN() noexcept {
  deep_copy(Vnode_n, Vnode_n0);
  deep_copy(rho_n, rho_n0);
  deep_copy(rhop_n, rhop_n0);
  deep_copy(V_n, V_n0);
  deep_copy(eps_n, eps_n0);
  deep_copy(epsp_n, epsp_n0);
  deep_copy(F_n, F_n0);
  deep_copy(Xpart_n, Xpart_n0);

  if (options->testCase == options->SodCase ||
      options->testCase == options->BiSodCase) {
    // const ℝ δt_init = 1.0e-4;
    options->deltat_init = 1.0e-4;
    deltat_n = options->deltat_init;
  } else if (options->testCase == options->BiShockBubble) {
    options->deltat_init = 1.e-7;
    deltat_n = 1.0e-7;
  } else if (options->testCase == options->SedovTestCase ||
             options->testCase == options->BiSedovTestCase) {
    // const ℝ δt_init = 1.0e-4;
    options->deltat_init = 1.0e-4;
    deltat_n = 1.0e-4;
  } else if (options->testCase == options->NohTestCase ||
             options->testCase == options->BiNohTestCase) {
    // const ℝ δt_init = 1.0e-4;
    options->deltat_init = 1.0e-4;
    deltat_n = 1.0e-4;
  } else if (options->testCase == options->TriplePoint ||
             options->testCase == options->BiTriplePoint) {
    // const ℝ δt_init = 1.0e-5; avec donnees adimensionnées
    options->deltat_init = 1.0e-5;  // avec pression de 1.e5 / 1.e-8
    deltat_n = 1.0e-5;
  }
  ETOTALE_0 = 0.;
  Kokkos::parallel_for(
      "initETOTALE", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        ETOT_0(cCells) =
            (rho_n0(cCells) * v(cCells)) *
            (eps_n0(cCells) + 0.5 * (V_n0(cCells)[0] * V_n0(cCells)[0] +
                                     V_n0(cCells)[1] * V_n0(cCells)[1]));
        MTOT_0(cCells) = (rho_n0(cCells) * v(cCells));
      });
  double reductionE(0.), reductionM(0.);
  {
    Kokkos::Sum<double> reducerE(reductionE);
    Kokkos::parallel_reduce("reductionE", nbCells,
                            KOKKOS_LAMBDA(const int& cCells, double& x) {
                              reducerE.join(x, ETOT_0(cCells));
                            },
                            reducerE);
    Kokkos::Sum<double> reducerM(reductionM);
    Kokkos::parallel_reduce("reductionM", nbCells,
                            KOKKOS_LAMBDA(const int& cCells, double& x) {
                              reducerM.join(x, MTOT_0(cCells));
                            },
                            reducerM);
  }
  ETOTALE_0 = reductionE;
  MASSET_0 = reductionM;
}
