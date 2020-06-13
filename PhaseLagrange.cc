#include <stdlib.h>                    // for exit
#include <Kokkos_Core.hpp>           // for KOKKOS_LAMBDA
#include <algorithm>                   // for equal, copy
#include <array>                       // for array, operator!=
#include <iostream>                    // for operator<<, basic_ostream::ope...
#include <limits>                      // for numeric_limits
#include <vector>                      // for vector, allocator
#include "EucclhydRemap.h"             // for EucclhydRemap, EucclhydRemap::...
#include "Utiles-Impl.h"               // for EucclhydRemap::tensProduct
#include "mesh/CartesianMesh2D.h"      // for CartesianMesh2D
#include "types/ArrayOperations.h"     // for plus, multiply, minus, divide
#include "types/MathFunctions.h"       // for max, min, dot, matVectProduct
#include "types/MultiArray.h"          // for operator<<
#include "utils/Utils.h"               // for indexOf


/**
 * Job computeCornerNormal called @1.0 in simulate method.
 * In variables: X
 * Out variables: lminus, lpc_n, lplus, nminus, nplus
 */
void EucclhydRemap::computeCornerNormal() noexcept {
  Kokkos::parallel_for(
      "computeCornerNormal", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        int cId(cCells);
        {
          auto nodesOfCellC(mesh->getNodesOfCell(cId));
          for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
               pNodesOfCellC++) {
            int pId(nodesOfCellC[pNodesOfCellC]);
            int pMinus1Id(nodesOfCellC[(pNodesOfCellC - 1 + nbNodesOfCell) %
                                       nbNodesOfCell]);
            int pPlus1Id(nodesOfCellC[(pNodesOfCellC + 1 + nbNodesOfCell) %
                                      nbNodesOfCell]);
            int cCellsOfNodeP(utils::indexOf(mesh->getCellsOfNode(pId), cId));
            int pMinus1Nodes(pMinus1Id);
            int pNodes(pId);
            int pPlus1Nodes(pPlus1Id);
            RealArray1D<dim> xp = X(pNodes);
            RealArray1D<dim> xpPlus = X(pPlus1Nodes);
            RealArray1D<dim> xpMinus = X(pMinus1Nodes);
            RealArray1D<dim> npc_plus;
            npc_plus[0] = 0.5 * (xpPlus[1] - xp[1]);
            npc_plus[1] = 0.5 * (xp[0] - xpPlus[0]);
            double lpc_plus = MathFunctions::norm(npc_plus);
            npc_plus = ArrayOperations::divide(npc_plus, lpc_plus);
            nplus(pNodes, cCellsOfNodeP) = npc_plus;
            lplus(pNodes, cCellsOfNodeP) = lpc_plus;
            RealArray1D<dim> npc_minus;
            npc_minus[0] = 0.5 * (xp[1] - xpMinus[1]);
            npc_minus[1] = 0.5 * (xpMinus[0] - xp[0]);
            double lpc_minus = MathFunctions::norm(npc_minus);
            npc_minus = ArrayOperations::divide(npc_minus, lpc_minus);
            nminus(pNodes, cCellsOfNodeP) = npc_minus;
            lminus(pNodes, cCellsOfNodeP) = lpc_minus;
            lpc_n(pNodes, cCellsOfNodeP) = ArrayOperations::plus(
                ArrayOperations::multiply(lpc_plus, npc_plus),
                ArrayOperations::multiply(lpc_minus, npc_minus));
          }
        }
      });
}

/**
 * Job computeEOS called @1.0 in executeTimeLoopN method.
 * In variables: eos, eosPerfectGas, eps_n, gamma, rho_n
 * Out variables: c, p
 */
void EucclhydRemap::computeEOS() noexcept {
  if (options->eos == options->eosPerfectGas)
    Kokkos::parallel_for(
        "computeEOS", nbCells, KOKKOS_LAMBDA(const int& cCells) {
          p(cCells) = 0.;
          for (imat = 0; imat < nbmatmax; imat++) {
            pp(cCells)[imat] = (options->gammap[imat] - 1.0) *
                               rhop_n(cCells)[imat] * epsp_n(cCells)[imat];
            p(cCells) += fracvol(cCells)[imat] * pp(cCells)[imat];
            if (rhop_n(cCells)[imat] > 0.) {
              vitsonp(cCells)[imat] =
                  MathFunctions::sqrt(options->gammap[imat] * pp(cCells)[imat] /
                                      rhop_n(cCells)[imat]);
            } else
              vitsonp(cCells)[imat] = 0.;
          }
          if (rho_n(cCells) > 0.) {
            vitson(cCells) =
                MathFunctions::sqrt(options->gamma * p(cCells) / rho_n(cCells));
          } else {
            std::cout << " cell " << cCells
                      << " densité moyenne négative ou nulle " << rho_n(cCells)
                      << std::endl;
          }
        });
}
/**
 * Job computeGradients called @1.0 in executeTimeLoopN method.
 * In variables: F_n, Vnode_n, lpc_n, spaceOrder, v
 * Out variables: gradV, gradp, gradp1, gradp2, gradp3
 */
void EucclhydRemap::computeGradients() noexcept {
  Kokkos::parallel_for(
      "computeDissipationMatrix", nbNodes, KOKKOS_LAMBDA(const int& pNodes) {
        int pId(pNodes);
        {
          for (imat = 0; imat < nbmatmax; imat++)
            fracvolnode(pNodes)[imat] = 0.;
          auto cellsOfNodeP(mesh->getCellsOfNode(pId));
          for (int cCellsOfNodeP = 0; cCellsOfNodeP < cellsOfNodeP.size();
               cCellsOfNodeP++) {
            int cId(cellsOfNodeP[cCellsOfNodeP]);
            int cCells(cId);
            for (imat = 0; imat < nbmatmax; imat++)
              fracvolnode(pNodes)[imat] += fracvol(cCells)[imat] * 0.25;
          }
        }
      });
  Kokkos::parallel_for(
      "computeGradients", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        int cId(cCells);
        RealArray1D<dim> reductionF1 = options->zeroVect;
        RealArray1D<dim> reductionF2 = options->zeroVect;
        RealArray1D<dim> reductionF3 = options->zeroVect;
        {
          auto nodesOfCellC(mesh->getNodesOfCell(cId));
          for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
               pNodesOfCellC++) {
            int pId(nodesOfCellC[pNodesOfCellC]);
            int cCellsOfNodeP(utils::indexOf(mesh->getCellsOfNode(pId), cId));
            int pNodes(pId);
            reductionF1 = ArrayOperations::plus(
                reductionF1,
                ArrayOperations::multiply(fracvolnode(pNodes)[0],
                                          lpc_n(pNodes, cCellsOfNodeP)));
            reductionF2 = ArrayOperations::plus(
                reductionF2,
                ArrayOperations::multiply(fracvolnode(pNodes)[1],
                                          lpc_n(pNodes, cCellsOfNodeP)));
            reductionF3 = ArrayOperations::plus(
                reductionF3,
                ArrayOperations::multiply(fracvolnode(pNodes)[2],
                                          lpc_n(pNodes, cCellsOfNodeP)));
          }
        }
        gradf1(cCells) = ArrayOperations::divide(reductionF1, v(cCells));
        gradf2(cCells) = ArrayOperations::divide(reductionF2, v(cCells));
        gradf3(cCells) = ArrayOperations::divide(reductionF3, v(cCells));
      });
  if (options->spaceOrder == 2)
    Kokkos::parallel_for(
        "computeGradients", nbCells, KOKKOS_LAMBDA(const int& cCells) {
          int cId(cCells);
          RealArray2D<dim, dim> reduction14 = options->zeroMat;
          {
            auto nodesOfCellC(mesh->getNodesOfCell(cId));
            for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
                 pNodesOfCellC++) {
              int pId(nodesOfCellC[pNodesOfCellC]);
              int cCellsOfNodeP(utils::indexOf(mesh->getCellsOfNode(pId), cId));
              int pNodes(pId);
              reduction14 = ArrayOperations::plus(
                  reduction14,
                  (tensProduct(Vnode_n(pNodes), lpc_n(pNodes, cCellsOfNodeP))));
            }
          }
          gradV(cCells) = ArrayOperations::divide(reduction14, v(cCells));
          RealArray1D<dim> reduction15 = options->zeroVect;
          RealArray1D<dim> reduction15a = options->zeroVect;
          RealArray1D<dim> reduction15b = options->zeroVect;
          RealArray1D<dim> reduction15c = options->zeroVect;
          {
            auto nodesOfCellC(mesh->getNodesOfCell(cId));
            for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
                 pNodesOfCellC++) {
              int pId(nodesOfCellC[pNodesOfCellC]);
              int cCellsOfNodeP(utils::indexOf(mesh->getCellsOfNode(pId), cId));
              int pNodes(pId);
              reduction15 = ArrayOperations::plus(reduction15,
                                                  (F_n(pNodes, cCellsOfNodeP)));

              reduction15a = ArrayOperations::plus(
                  reduction15a, (F1_n(pNodes, cCellsOfNodeP)));

              reduction15b = ArrayOperations::plus(
                  reduction15b, (F2_n(pNodes, cCellsOfNodeP)));

              reduction15c = ArrayOperations::plus(
                  reduction15c, (F3_n(pNodes, cCellsOfNodeP)));
            }
          }
          gradp(cCells) = ArrayOperations::divide(reduction15, v(cCells));
          gradp1(cCells) = ArrayOperations::divide(reduction15a, v(cCells));
          gradp2(cCells) = ArrayOperations::divide(reduction15b, v(cCells));
          gradp3(cCells) = ArrayOperations::divide(reduction15c, v(cCells));
        });
}

/**
 * Job computeMass called @1.0 in executeTimeLoopN method.
 * In variables: rho_n, v
 * Out variables: m
 */
void EucclhydRemap::computeMass() noexcept {
  Kokkos::parallel_for(
      "computeMass", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        m(cCells) = rho_n(cCells) * v(cCells);
        for (imat = 0; imat < nbmatmax; imat++)
          mp(cCells)[imat] = fracmass(cCells)[imat] * m(cCells);
      });
}
/**
 * Job computeDissipationMatrix called @2.0 in executeTimeLoopN method.
 * In variables: c, lminus, lplus, nminus, nplus, rho_n
 * Out variables: M
 */
void EucclhydRemap::computeDissipationMatrix() noexcept {
  Kokkos::parallel_for(
      "computeDissipationMatrix", nbNodes, KOKKOS_LAMBDA(const int& pNodes) {
        int pId(pNodes);
        {
          auto cellsOfNodeP(mesh->getCellsOfNode(pId));
          for (int cCellsOfNodeP = 0; cCellsOfNodeP < cellsOfNodeP.size();
               cCellsOfNodeP++) {
            int cId(cellsOfNodeP[cCellsOfNodeP]);
            int cCells(cId);
            RealArray2D<dim, dim> cornerMatrix = ArrayOperations::plus(
                ArrayOperations::multiply(
                    lplus(pNodes, cCellsOfNodeP),
                    tensProduct(nplus(pNodes, cCellsOfNodeP),
                                nplus(pNodes, cCellsOfNodeP))),
                ArrayOperations::multiply(
                    lminus(pNodes, cCellsOfNodeP),
                    tensProduct(nminus(pNodes, cCellsOfNodeP),
                                nminus(pNodes, cCellsOfNodeP))));
            M(pNodes, cCellsOfNodeP) = ArrayOperations::multiply(
                rho_n(cCells) * vitson(cCells), cornerMatrix);

            M1(pNodes, cCellsOfNodeP) = ArrayOperations::multiply(
                rhop_n(cCells)[0] * vitson(cCells), cornerMatrix);
            M2(pNodes, cCellsOfNodeP) = ArrayOperations::multiply(
                rhop_n(cCells)[1] * vitson(cCells), cornerMatrix);
            M3(pNodes, cCellsOfNodeP) = ArrayOperations::multiply(
                rhop_n(cCells)[2] * vitson(cCells), cornerMatrix);
          }
        }
      });
}
/**
 * Job computedeltatc called @2.0 in executeTimeLoopN method.
 * In variables: V_n, c, perim, v
 * Out variables: deltatc
 */
void EucclhydRemap::computedeltatc() noexcept {
  Kokkos::parallel_for(
      "computedeltatc", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        deltatc(cCells) =
            v(cCells) / (perim(cCells) *
                         (MathFunctions::norm(V_n(cCells)) + vitson(cCells)));
      });
}
/**
 * Job extrapolateValue called @2.0 in executeTimeLoopN method.
 * In variables: V_n, X, Xc, gradV, gradp, p, spaceOrder
 * Out variables: V_extrap, p_extrap, pp_extrap
 */
void EucclhydRemap::extrapolateValue() noexcept {
  if (options->spaceOrder == 1) {
    Kokkos::parallel_for(
        "extrapolateValue", nbCells, KOKKOS_LAMBDA(const int& cCells) {
          int cId(cCells);
          {
            auto nodesOfCellC(mesh->getNodesOfCell(cId));
            for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
                 pNodesOfCellC++) {
              V_extrap(cCells, pNodesOfCellC) = V_n(cCells);
              p_extrap(cCells, pNodesOfCellC) = p(cCells);
              for (imat = 0; imat < nbmatmax; imat++)
                pp_extrap(cCells, pNodesOfCellC)[imat] = pp(cCells)[imat];
            }
          }
        });
  } else {
    Kokkos::parallel_for(
        "extrapolateValue", nbNodes, KOKKOS_LAMBDA(const int& pNodes) {
          int pId(pNodes);
          double reduction16 = numeric_limits<double>::max();
          double reduction16a = numeric_limits<double>::max();
          double reduction16b = numeric_limits<double>::max();
          double reduction16c = numeric_limits<double>::max();
          {
            auto cellsOfNodeP(mesh->getCellsOfNode(pId));
            for (int dCellsOfNodeP = 0; dCellsOfNodeP < cellsOfNodeP.size();
                 dCellsOfNodeP++) {
              int dId(cellsOfNodeP[dCellsOfNodeP]);
              int dCells(dId);
              reduction16 = MathFunctions::min(reduction16, p(dCells));
              reduction16a = MathFunctions::min(reduction16a, pp(dCells)[0]);
              reduction16b = MathFunctions::min(reduction16b, pp(dCells)[1]);
              reduction16c = MathFunctions::min(reduction16c, pp(dCells)[2]);
            }
          }
          double minP = reduction16;
          double minP1 = reduction16a;
          double minP2 = reduction16b;
          double minP3 = reduction16c;
          double reduction17 = numeric_limits<double>::min();
          double reduction17a = numeric_limits<double>::min();
          double reduction17b = numeric_limits<double>::min();
          double reduction17c = numeric_limits<double>::min();
          {
            auto cellsOfNodeP(mesh->getCellsOfNode(pId));
            for (int dCellsOfNodeP = 0; dCellsOfNodeP < cellsOfNodeP.size();
                 dCellsOfNodeP++) {
              int dId(cellsOfNodeP[dCellsOfNodeP]);
              int dCells(dId);
              reduction17 = MathFunctions::max(reduction17, p(dCells));
              reduction17a = MathFunctions::max(reduction17a, pp(dCells)[0]);
              reduction17b = MathFunctions::max(reduction17b, pp(dCells)[1]);
              reduction17c = MathFunctions::max(reduction17c, pp(dCells)[2]);
            }
          }
          double maxP = reduction17;
          double maxP1 = reduction17a;
          double maxP2 = reduction17b;
          double maxP3 = reduction17c;
          double reduction18 = numeric_limits<double>::max();
          {
            auto cellsOfNodeP(mesh->getCellsOfNode(pId));
            for (int dCellsOfNodeP = 0; dCellsOfNodeP < cellsOfNodeP.size();
                 dCellsOfNodeP++) {
              int dId(cellsOfNodeP[dCellsOfNodeP]);
              int dCells(dId);
              reduction18 = MathFunctions::min(reduction18, V_n(dCells)[0]);
            }
          }
          double minVx = reduction18;
          double reduction19 = numeric_limits<double>::min();
          {
            auto cellsOfNodeP(mesh->getCellsOfNode(pId));
            for (int dCellsOfNodeP = 0; dCellsOfNodeP < cellsOfNodeP.size();
                 dCellsOfNodeP++) {
              int dId(cellsOfNodeP[dCellsOfNodeP]);
              int dCells(dId);
              reduction19 = MathFunctions::max(reduction19, V_n(dCells)[0]);
            }
          }
          double maxVx = reduction19;
          double reduction20 = numeric_limits<double>::max();
          {
            auto cellsOfNodeP(mesh->getCellsOfNode(pId));
            for (int dCellsOfNodeP = 0; dCellsOfNodeP < cellsOfNodeP.size();
                 dCellsOfNodeP++) {
              int dId(cellsOfNodeP[dCellsOfNodeP]);
              int dCells(dId);
              reduction20 = MathFunctions::min(reduction20, V_n(dCells)[1]);
            }
          }
          double minVy = reduction20;
          double reduction21 = numeric_limits<double>::min();
          {
            auto cellsOfNodeP(mesh->getCellsOfNode(pId));
            for (int dCellsOfNodeP = 0; dCellsOfNodeP < cellsOfNodeP.size();
                 dCellsOfNodeP++) {
              int dId(cellsOfNodeP[dCellsOfNodeP]);
              int dCells(dId);
              reduction21 = MathFunctions::max(reduction21, V_n(dCells)[1]);
            }
          }
          double maxVy = reduction21;
          {
            auto cellsOfNodeP(mesh->getCellsOfNode(pId));
            for (int cCellsOfNodeP = 0; cCellsOfNodeP < cellsOfNodeP.size();
                 cCellsOfNodeP++) {
              int cId(cellsOfNodeP[cCellsOfNodeP]);
              int cCells(cId);
              int pNodesOfCellC(utils::indexOf(mesh->getNodesOfCell(cId), pId));

              // double ptmp = p(cCells) + MathFunctions::dot(gradp(cCells),
              // ArrayOperations::minus(X(pNodes), Xc(cCells)));
              // p_extrap(cCells,pNodesOfCellC) =
              // MathFunctions::max(MathFunctions::min(maxP, ptmp), minP);

              // pour chaque matériau,
              double ptmp1 = pp(cCells)[0] +
                             MathFunctions::dot(
                                 gradp1(cCells),
                                 ArrayOperations::minus(X(pNodes), Xc(cCells)));
              pp_extrap(cCells, pNodesOfCellC)[0] =
                  MathFunctions::max(MathFunctions::min(maxP1, ptmp1), minP1);
              double ptmp2 = pp(cCells)[1] +
                             MathFunctions::dot(
                                 gradp2(cCells),
                                 ArrayOperations::minus(X(pNodes), Xc(cCells)));
              pp_extrap(cCells, pNodesOfCellC)[1] =
                  MathFunctions::max(MathFunctions::min(maxP2, ptmp2), minP2);
              double ptmp3 = pp(cCells)[2] +
                             MathFunctions::dot(
                                 gradp3(cCells),
                                 ArrayOperations::minus(X(pNodes), Xc(cCells)));
              pp_extrap(cCells, pNodesOfCellC)[2] =
                  MathFunctions::max(MathFunctions::min(maxP3, ptmp3), minP3);

              p_extrap(cCells, pNodesOfCellC) = 0.;
              // et on recalcule la moyenne
              for (imat = 0; imat < nbmatmax; imat++)
                p_extrap(cCells, pNodesOfCellC) +=
                    fracvol(cCells)[imat] *
                    pp_extrap(cCells, pNodesOfCellC)[imat];

              RealArray1D<dim> Vtmp = ArrayOperations::plus(
                  V_n(cCells), MathFunctions::matVectProduct(
                                   gradV(cCells), ArrayOperations::minus(
                                                      X(pNodes), Xc(cCells))));
              V_extrap(cCells, pNodesOfCellC)[0] =
                  MathFunctions::max(MathFunctions::min(maxVx, Vtmp[0]), minVx);
              V_extrap(cCells, pNodesOfCellC)[1] =
                  MathFunctions::max(MathFunctions::min(maxVy, Vtmp[1]), minVy);
            }
          }
        });
  }
}
/**
 * Job computeG called @3.0 in executeTimeLoopN method.
 * In variables: M, V_extrap, lpc_n, p_extrap
 * Out variables: G
 */
void EucclhydRemap::computeG() noexcept {
  Kokkos::parallel_for(
      "computeG", nbNodes, KOKKOS_LAMBDA(const int& pNodes) {
        int pId(pNodes);
        RealArray1D<dim> reduction1 = options->zeroVect;
        {
          auto cellsOfNodeP(mesh->getCellsOfNode(pId));
          for (int cCellsOfNodeP = 0; cCellsOfNodeP < cellsOfNodeP.size();
               cCellsOfNodeP++) {
            int cId(cellsOfNodeP[cCellsOfNodeP]);
            int cCells(cId);
            int pNodesOfCellC(utils::indexOf(mesh->getNodesOfCell(cId), pId));
            reduction1 = ArrayOperations::plus(
                reduction1,
                (ArrayOperations::plus(
                    MathFunctions::matVectProduct(
                        M(pNodes, cCellsOfNodeP),
                        V_extrap(cCells, pNodesOfCellC)),
                    ArrayOperations::multiply(p_extrap(cCells, pNodesOfCellC),
                                              lpc_n(pNodes, cCellsOfNodeP)))));
          }
        }
        G(pNodes) = reduction1;
      });
}

/**
 * Job computeNodeDissipationMatrixAndG called @3.0 in executeTimeLoopN method.
 * In variables: M
 * Out variables: Mnode
 */
void EucclhydRemap::computeNodeDissipationMatrixAndG() noexcept {
  Kokkos::parallel_for(
      "computeNodeDissipationMatrixAndG", nbNodes,
      KOKKOS_LAMBDA(const int& pNodes) {
        int pId(pNodes);
        RealArray2D<dim, dim> reduction0 = options->zeroMat;
        {
          auto cellsOfNodeP(mesh->getCellsOfNode(pId));
          for (int cCellsOfNodeP = 0; cCellsOfNodeP < cellsOfNodeP.size();
               cCellsOfNodeP++) {
            reduction0 =
                ArrayOperations::plus(reduction0, (M(pNodes, cCellsOfNodeP)));
          }
        }
        Mnode(pNodes) = reduction0;
      });
}
/**
 * Job computeNodeVelocity called @4.0 in executeTimeLoopN method.
 * In variables: G, Mnode
 * Out variables: Vnode_nplus1
 */
void EucclhydRemap::computeNodeVelocity() noexcept {
  auto innerNodes(mesh->getInnerNodes());
  Kokkos::parallel_for(
      "computeNodeVelocity", nbInnerNodes,
      KOKKOS_LAMBDA(const int& pInnerNodes) {
        int pId(innerNodes[pInnerNodes]);
        int pNodes(pId);
        Vnode_nplus1(pNodes) =
            MathFunctions::matVectProduct(inverse(Mnode(pNodes)), G(pNodes));
      });
}
/**
 * Job computeFaceVelocity called @5.0 in executeTimeLoopN method.
 * In variables: Vnode_nplus1, faceNormal
 * Out variables: faceNormalVelocity
 */
void EucclhydRemap::computeFaceVelocity() noexcept {
  auto faces(mesh->getFaces());
  Kokkos::parallel_for(
      "computeFaceVelocity", nbFaces, KOKKOS_LAMBDA(const int& fFaces) {
        int fId(faces[fFaces]);
        RealArray1D<dim> reduction5 = options->zeroVect;
        {
          auto nodesOfFaceF(mesh->getNodesOfFace(fId));
          for (int pNodesOfFaceF = 0; pNodesOfFaceF < nodesOfFaceF.size();
               pNodesOfFaceF++) {
            int pId(nodesOfFaceF[pNodesOfFaceF]);
            int pNodes(pId);
            reduction5 =
                ArrayOperations::plus(reduction5, (Vnode_nplus1(pNodes)));
          }
        }
        faceNormalVelocity(fFaces) = MathFunctions::dot(
            ArrayOperations::multiply(0.5, reduction5), faceNormal(fFaces));
      });
}

/**
 * Job computeLagrangePosition called @5.0 in executeTimeLoopN method.
 * In variables: Vnode_nplus1, X, deltat_n
 * Out variables: XLagrange
 */
void EucclhydRemap::computeLagrangePosition() noexcept {
  Kokkos::parallel_for(
      "computeLagrangePosition", nbNodes, KOKKOS_LAMBDA(const int& pNodes) {
        XLagrange(pNodes) = ArrayOperations::plus(
            X(pNodes),
            ArrayOperations::multiply(Vnode_nplus1(pNodes), deltat_n));
      });
  auto faces(mesh->getFaces());
  Kokkos::parallel_for(
      "computeLagrangePosition", nbFaces, KOKKOS_LAMBDA(const int& fFaces) {
        int fId(faces[fFaces]);
        int n1FirstNodeOfFaceF(mesh->getFirstNodeOfFace(fId));
        int n1Id(n1FirstNodeOfFaceF);
        int n1Nodes(n1Id);
        int n2SecondNodeOfFaceF(mesh->getSecondNodeOfFace(fId));
        int n2Id(n2SecondNodeOfFaceF);
        int n2Nodes(n2Id);
        RealArray1D<dim> X_face = ArrayOperations::multiply(
            0.5,
            (ArrayOperations::plus(XLagrange(n1Nodes), XLagrange(n2Nodes))));
        RealArray1D<dim> face_vec =
            ArrayOperations::minus(XLagrange(n2Nodes), XLagrange(n1Nodes));
        XfLagrange(fFaces) = X_face;
        faceLengthLagrange(fFaces) = MathFunctions::norm(face_vec);
      });
}

/**
 * Job computeSubCellForce called @5.0 in executeTimeLoopN method.
 * In variables: M, V_extrap, Vnode_nplus1, lpc_n, p_extrap
 * Out variables: F_nplus1
 */
void EucclhydRemap::computeSubCellForce() noexcept {
  Kokkos::parallel_for(
      "computeSubCellForce", nbNodes, KOKKOS_LAMBDA(const int& pNodes) {
        int pId(pNodes);
        {
          auto cellsOfNodeP(mesh->getCellsOfNode(pId));
          for (int cCellsOfNodeP = 0; cCellsOfNodeP < cellsOfNodeP.size();
               cCellsOfNodeP++) {
            int cId(cellsOfNodeP[cCellsOfNodeP]);
            int cCells(cId);
            int pNodesOfCellC(utils::indexOf(mesh->getNodesOfCell(cId), pId));
            F_nplus1(pNodes, cCellsOfNodeP) = ArrayOperations::plus(
                ArrayOperations::multiply(-p_extrap(cCells, pNodesOfCellC),
                                          lpc_n(pNodes, cCellsOfNodeP)),
                MathFunctions::matVectProduct(
                    M(pNodes, cCellsOfNodeP),
                    ArrayOperations::minus(Vnode_nplus1(pNodes),
                                           V_extrap(cCells, pNodesOfCellC))));

            F1_nplus1(pNodes, cCellsOfNodeP) = ArrayOperations::plus(
                ArrayOperations::multiply(-pp_extrap(cCells, pNodesOfCellC)[0],
                                          lpc_n(pNodes, cCellsOfNodeP)),
                MathFunctions::matVectProduct(
                    M1(pNodes, cCellsOfNodeP),
                    ArrayOperations::minus(Vnode_nplus1(pNodes),
                                           V_extrap(cCells, pNodesOfCellC))));

            F2_nplus1(pNodes, cCellsOfNodeP) = ArrayOperations::plus(
                ArrayOperations::multiply(-pp_extrap(cCells, pNodesOfCellC)[1],
                                          lpc_n(pNodes, cCellsOfNodeP)),
                MathFunctions::matVectProduct(
                    M2(pNodes, cCellsOfNodeP),
                    ArrayOperations::minus(Vnode_nplus1(pNodes),
                                           V_extrap(cCells, pNodesOfCellC))));

            F3_nplus1(pNodes, cCellsOfNodeP) = ArrayOperations::plus(
                ArrayOperations::multiply(-pp_extrap(cCells, pNodesOfCellC)[2],
                                          lpc_n(pNodes, cCellsOfNodeP)),
                MathFunctions::matVectProduct(
                    M3(pNodes, cCellsOfNodeP),
                    ArrayOperations::minus(Vnode_nplus1(pNodes),
                                           V_extrap(cCells, pNodesOfCellC))));
          }
        }
      });
}
/**
 * Job computeLagrangeVolumeAndCenterOfGravity called @6.0 in executeTimeLoopN
 * method. In variables: XLagrange Out variables: XcLagrange, vLagrange
 */
void EucclhydRemap::computeLagrangeVolumeAndCenterOfGravity() noexcept {
  Kokkos::parallel_for(
      "computeLagrangeVolumeAndCenterOfGravity", nbCells,
      KOKKOS_LAMBDA(const int& cCells) {
        int cId(cCells);
        double reduction6 = 0.0;
        {
          auto nodesOfCellC(mesh->getNodesOfCell(cId));
          for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
               pNodesOfCellC++) {
            int pId(nodesOfCellC[pNodesOfCellC]);
            int pPlus1Id(nodesOfCellC[(pNodesOfCellC + 1 + nbNodesOfCell) %
                                      nbNodesOfCell]);
            int pNodes(pId);
            int pPlus1Nodes(pPlus1Id);
            reduction6 = reduction6 + (crossProduct2d(XLagrange(pNodes),
                                                      XLagrange(pPlus1Nodes)));
          }
        }
        double vol = 0.5 * reduction6;
        vLagrange(cCells) = vol;
        for (imat = 0; imat < nbmatmax; imat++)
          vpLagrange(cCells)[imat] = fracvol(cCells)[imat] * vol;
        RealArray1D<dim> reduction7 = options->zeroVect;
        {
          auto nodesOfCellC(mesh->getNodesOfCell(cId));
          for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
               pNodesOfCellC++) {
            int pId(nodesOfCellC[pNodesOfCellC]);
            int pPlus1Id(nodesOfCellC[(pNodesOfCellC + 1 + nbNodesOfCell) %
                                      nbNodesOfCell]);
            int pNodes(pId);
            int pPlus1Nodes(pPlus1Id);
            reduction7 = ArrayOperations::plus(
                reduction7,
                (ArrayOperations::multiply(
                    crossProduct2d(XLagrange(pNodes), XLagrange(pPlus1Nodes)),
                    (ArrayOperations::plus(XLagrange(pNodes),
                                           XLagrange(pPlus1Nodes))))));
          }
        }
        XcLagrange(cCells) =
            ArrayOperations::multiply(1.0 / (6.0 * vol), reduction7);
      });
}
/**
 * Job computeFacedeltaxLagrange called @7.0 in executeTimeLoopN method.
 * In variables: XcLagrange, faceNormal
 * Out variables: deltaxLagrange
 */
void EucclhydRemap::computeFacedeltaxLagrange() noexcept {
  auto faces(mesh->getFaces());
  Kokkos::parallel_for(
      "computeFacedeltaxLagrange", nbFaces, KOKKOS_LAMBDA(const int& fFaces) {
        int fId(faces[fFaces]);
        int cfFrontCellF(mesh->getFrontCell(fId));
        int cfId(cfFrontCellF);
        int cfCells(cfId);
        int cbBackCellF(mesh->getBackCell(fId));
        int cbId(cbBackCellF);
        int cbCells(cbId);
        deltaxLagrange(fFaces) = MathFunctions::dot(
            ArrayOperations::minus(XcLagrange(cfCells), XcLagrange(cbCells)),
            faceNormal(fFaces));
      });
}

/**
 * Job updateCellCenteredLagrangeVariables called @7.0 in executeTimeLoopN
 * method. In variables: F_nplus1, V_n, Vnode_nplus1, deltat_n, eps_n, lpc_n, m,
 * rho_n, vLagrange Out variables: ULagrange
 */
void EucclhydRemap::updateCellCenteredLagrangeVariables() noexcept {
  Kokkos::parallel_for(
      "updateCellCenteredLagrangeVariables", nbCells,
      KOKKOS_LAMBDA(const int& cCells) {
        int cId(cCells);
        double reduction2 = 0.0;
        {
          auto nodesOfCellC(mesh->getNodesOfCell(cId));
          for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
               pNodesOfCellC++) {
            int pId(nodesOfCellC[pNodesOfCellC]);
            int cCellsOfNodeP(utils::indexOf(mesh->getCellsOfNode(pId), cId));
            int pNodes(pId);
            reduction2 =
                reduction2 + (MathFunctions::dot(lpc_n(pNodes, cCellsOfNodeP),
                                                 Vnode_nplus1(pNodes)));
          }
        }
        double rhoLagrange =
            1 / (1 / rho_n(cCells) + deltat_n / m(cCells) * reduction2);
        RealArray1D<dim> reduction3 = options->zeroVect;
        {
          auto nodesOfCellC(mesh->getNodesOfCell(cId));
          for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
               pNodesOfCellC++) {
            int pId(nodesOfCellC[pNodesOfCellC]);
            int cCellsOfNodeP(utils::indexOf(mesh->getCellsOfNode(pId), cId));
            int pNodes(pId);
            reduction3 = ArrayOperations::plus(
                reduction3, (F_nplus1(pNodes, cCellsOfNodeP)));
          }
        }
        ForceGradp(cCells) = ArrayOperations::divide(reduction3, v(cCells));
        RealArray1D<dim> VLagrange = ArrayOperations::plus(
            V_n(cCells),
            ArrayOperations::multiply(deltat_n / m(cCells), reduction3));
        double reduction4 = 0.0;
        RealArray1D<nbmatmax> preduction4 = options->zeroVectmat;
        {
          auto nodesOfCellC(mesh->getNodesOfCell(cId));
          for (int pNodesOfCellC = 0; pNodesOfCellC < nodesOfCellC.size();
               pNodesOfCellC++) {
            int pId(nodesOfCellC[pNodesOfCellC]);
            int cCellsOfNodeP(utils::indexOf(mesh->getCellsOfNode(pId), cId));
            int pNodes(pId);
            reduction4 =
                reduction4 + (MathFunctions::dot(
                                 F_nplus1(pNodes, cCellsOfNodeP),
                                 ArrayOperations::minus(
                                     Vnode_nplus1(pNodes),
                                     ArrayOperations::multiply(
                                         0.5, (ArrayOperations::plus(
                                                  V_n(cCells), VLagrange))))));
            preduction4[0] =
                preduction4[0] +
                (MathFunctions::dot(
                    F1_nplus1(pNodes, cCellsOfNodeP),
                    ArrayOperations::minus(
                        Vnode_nplus1(pNodes),
                        ArrayOperations::multiply(
                            0.5,
                            (ArrayOperations::plus(V_n(cCells), VLagrange))))));
            preduction4[1] =
                preduction4[1] +
                (MathFunctions::dot(
                    F2_nplus1(pNodes, cCellsOfNodeP),
                    ArrayOperations::minus(
                        Vnode_nplus1(pNodes),
                        ArrayOperations::multiply(
                            0.5,
                            (ArrayOperations::plus(V_n(cCells), VLagrange))))));
            preduction4[2] =
                preduction4[2] +
                (MathFunctions::dot(
                    F3_nplus1(pNodes, cCellsOfNodeP),
                    ArrayOperations::minus(
                        Vnode_nplus1(pNodes),
                        ArrayOperations::multiply(
                            0.5,
                            (ArrayOperations::plus(V_n(cCells), VLagrange))))));
            if (F1_nplus1(pNodes, cCellsOfNodeP) !=
                F1_nplus1(pNodes, cCellsOfNodeP))
              std::cout << " cell   " << cCells
                        << " F1_nplus1(pNodes,cCellsOfNodeP) "
                        << F1_nplus1(pNodes, cCellsOfNodeP) << std::endl;
            if (F2_nplus1(pNodes, cCellsOfNodeP) !=
                F2_nplus1(pNodes, cCellsOfNodeP))
              std::cout << " cell   " << cCells
                        << " F2_nplus1(pNodes,cCellsOfNodeP) "
                        << F2_nplus1(pNodes, cCellsOfNodeP) << std::endl;
            if (F3_nplus1(pNodes, cCellsOfNodeP) !=
                F3_nplus1(pNodes, cCellsOfNodeP))
              std::cout << " cell   " << cCells
                        << " F3_nplus1(pNodes,cCellsOfNodeP) "
                        << F3_nplus1(pNodes, cCellsOfNodeP) << std::endl;
            if (V_n(cCells) != V_n(cCells))
              std::cout << " cell   " << cCells << " V_n(cCells) "
                        << V_n(cCells) << std::endl;
            if (Vnode_nplus1(pNodes) != Vnode_nplus1(pNodes))
              std::cout << " cell   " << cCells << " Vnode_nplus1(pNodes) "
                        << Vnode_nplus1(pNodes) << std::endl;
            if (VLagrange != VLagrange)
              std::cout << " cell   " << cCells << " VLagrange " << VLagrange
                        << std::endl;
          }
        }

        double epsLagrange = eps_n(cCells) + deltat_n / m(cCells) * reduction4;
        RealArray1D<nbmatmax> pepsLagrange;
        RealArray1D<nbmatmax> pepsLagrangec;
        for (imat = 0; imat < nbmatmax; imat++) {
          pepsLagrange[imat] = 0.;
          pepsLagrangec[imat] = 0.;
          if (fracvol(cCells)[imat] > options->threshold &&
              mp(cCells)[imat] != 0.)
            pepsLagrange[imat] =
                epsp_n(cCells)[imat] + fracvol(cCells)[imat] * deltat_n /
                                           mp(cCells)[imat] * preduction4[imat];
        }
        for (imat = 0; imat < nbmatmax; imat++) {
          ULagrange(cCells)[imat] = vpLagrange(cCells)[imat];

          ULagrange(cCells)[nbmatmax + imat] =
              fracmass(cCells)[imat] * vLagrange(cCells) * rhoLagrange;

          ULagrange(cCells)[2 * nbmatmax + imat + 2] =
              fracmass(cCells)[imat] * vLagrange(cCells) * rhoLagrange *
              pepsLagrange[imat];
        }

        ULagrange(cCells)[2 * nbmatmax] =
            vLagrange(cCells) * rhoLagrange * VLagrange[0];
        ULagrange(cCells)[2 * nbmatmax + 1] =
            vLagrange(cCells) * rhoLagrange * VLagrange[1];
        // projection de l'energie cinétique
        if (options->projectionConservative == 1)
          ULagrange(cCells)[3 * nbmatmax + 2] =
              0.5 * vLagrange(cCells) * rhoLagrange *
              (VLagrange[0] * VLagrange[0] + VLagrange[1] * VLagrange[1]);

        if (options->projectionAvecPlateauPente == 1) {
          // option ou on ne regarde pas la variation de rho, V et e
          // phi = (f1, f2, rho1, rho2, Vx, Vy, e1, e2
          // ce qui permet d'ecrire le flux telque
          // Flux = (dv1 = f1dv, dv2=f2*dv, dm1=rho1*df1, dm2=rho2*df2, d(mV) =
          // V*(dm1+dm2), d(m1e1) = e1*dm1,  d(m2e2) = e2*dm2 dans
          // computeIntersectionPP

          double somme_volume = 0.;
          for (imat = 0; imat < nbmatmax; imat++) {
            somme_volume += ULagrange(cCells)[imat];
          }
          // Phi volume
          double somme_masse = 0.;
          for (imat = 0; imat < nbmatmax; imat++) {
            Phi(cCells)[imat] = ULagrange(cCells)[imat] / somme_volume;

            // Phi masse
            if (ULagrange(cCells)[imat] != 0.)
              Phi(cCells)[nbmatmax + imat] =
                  ULagrange(cCells)[nbmatmax + imat] / ULagrange(cCells)[imat];
            else
              Phi(cCells)[nbmatmax + imat] = 0.;
            somme_masse += ULagrange(cCells)[nbmatmax + imat];
          }
          // Phi Vitesse
          Phi(cCells)[2 * nbmatmax] =
              ULagrange(cCells)[2 * nbmatmax] / somme_masse;
          Phi(cCells)[2 * nbmatmax + 1] =
              ULagrange(cCells)[2 * nbmatmax + 1] / somme_masse;
          // Phi energie
          for (imat = 0; imat < nbmatmax; imat++) {
            if (ULagrange(cCells)[nbmatmax + imat] != 0.)
              Phi(cCells)[2 * nbmatmax + imat + 2] =
                  ULagrange(cCells)[2 * nbmatmax + imat + 2] /
                  ULagrange(cCells)[nbmatmax + imat];
            else
              Phi(cCells)[2 * nbmatmax + imat + 2] = 0.;
          }
          // Phi energie cinétique
          if (options->projectionConservative == 1)
            Phi(cCells)[3 * nbmatmax + 2] =
                ULagrange(cCells)[3 * nbmatmax + 2] / somme_masse;

        } else {
          Phi(cCells) =
              ArrayOperations::divide(ULagrange(cCells), vLagrange(cCells));
        }

        if ((cCells == dbgcell3 || cCells == dbgcell2 || cCells == dbgcell1) &&
            test_debug == 1) {
          std::cout << " Apres Phase Lagrange cell   " << cCells << "Phi"
                    << Phi(cCells) << std::endl;
          std::cout << " cell   " << cCells << "ULagrange " << ULagrange(cCells)
                    << std::endl;
        }
        if (ULagrange(cCells) != ULagrange(cCells)) {
          std::cout << " cell   " << cCells << " Ulagrange "
                    << ULagrange(cCells) << std::endl;
          std::cout << " cell   " << cCells << " f1 " << fracvol(cCells)[0]
                    << " f2 " << fracvol(cCells)[1] << " f3 "
                    << fracvol(cCells)[2] << std::endl;
          std::cout << " cell   " << cCells << " c1 " << fracmass(cCells)[0]
                    << " c2 " << fracmass(cCells)[1] << " c3 "
                    << fracmass(cCells)[2] << std::endl;
          std::cout << " cell   " << cCells << " m1 " << mp(cCells)[0] << " m2 "
                    << mp(cCells)[1] << " m3 " << mp(cCells)[2] << std::endl;
          std::cout << " cell   " << cCells << " pepsLagrange[0] "
                    << pepsLagrange[0] << " preduction4[0] " << preduction4[0]
                    << " epsp_n[0] " << epsp_n(cCells)[0] << std::endl;
          std::cout << " cell   " << cCells << " pepsLagrange[1] "
                    << pepsLagrange[1] << " preduction4[1] " << preduction4[1]
                    << " epsp_n[1] " << epsp_n(cCells)[1] << std::endl;
          std::cout << " cell   " << cCells << " pepsLagrange[2] "
                    << pepsLagrange[2] << " preduction4[2] " << preduction4[2]
                    << " epsp_n[2] " << epsp_n(cCells)[2] << std::endl;
          std::cout << " densites  1 " << rhop_n(cCells)[0] << " 2 "
                    << rhop_n(cCells)[1] << " 3 " << rhop_n(cCells)[2]
                    << std::endl;
          exit(1);
        }

        // ETOT_L(cCells) = (rhoLagrange * vLagrange(cCells)) * (epsLagrange +
        // 0.5 * (VLagrange[0] * VLagrange[0] + VLagrange[1] * VLagrange[1]));
        MTOT_L(cCells) = 0.;
        ETOT_L(cCells) =
            (rhoLagrange * vLagrange(cCells)) *
            (fracmass(cCells)[0] * pepsLagrange[0] +
             fracmass(cCells)[1] * pepsLagrange[1] +
             fracmass(cCells)[2] * pepsLagrange[2] +
             0.5 * (VLagrange[0] * VLagrange[0] + VLagrange[1] * VLagrange[1]));
        for (imat = 0; imat < nbmatmax; imat++) {
          MTOT_L(cCells) +=
              fracmass(cCells)[imat] * (rhoLagrange * vLagrange(cCells));
        }
      });
  double reductionE(0.), reductionM(0.);
  {
    Kokkos::Sum<double> reducerE(reductionE);
    Kokkos::parallel_reduce(
        "reductionE", nbCells,
        KOKKOS_LAMBDA(const int& cCells, double& x) {
          reducerE.join(x, ETOT_L(cCells));
        },
        reducerE);
    Kokkos::Sum<double> reducerM(reductionM);
    Kokkos::parallel_reduce(
        "reductionM", nbCells,
        KOKKOS_LAMBDA(const int& cCells, double& x) {
          reducerM.join(x, MTOT_L(cCells));
        },
        reducerM);
  }
  ETOTALE_L = reductionE;
  MASSET_L = reductionM;
}
