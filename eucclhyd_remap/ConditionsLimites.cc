#include <Kokkos_Core.hpp>
#include <algorithm>                // for copy
#include <array>                    // for array
#include <iostream>                 // for operator<<, basic_ostream::operat...
#include <vector>                   // for allocator, vector
#include "EucclhydRemap.h"          // for EucclhydRemap, EucclhydRemap::Opt...
#include "UtilesRemap-Impl.h"       // for EucclhydRemap::computeRemapFlux
#include "mesh/CartesianMesh2D.h"   // for CartesianMesh2D
#include "types/ArrayOperations.h"  // for multiply, plus
#include "types/MathFunctions.h"    // for dot, matVectProduct, norm
#include "types/MultiArray.h"       // for operator<<
#include "utils/Utils.h"            // for indexOf

/**
 * Job computeBoundaryNodeVelocities called @4.0 in executeTimeLoopN method.
 * In variables: G, Mnode, bottomBC, bottomBCValue, leftBC, leftBCValue,
 * rightBC, rightBCValue, topBC, topBCValue Out variables: Vnode_nplus1
 */
void EucclhydRemap::computeBoundaryNodeVelocities() noexcept {
  auto leftNodes(mesh->getLeftNodes());
  Kokkos::parallel_for("computeBoundaryNodeVelocities", nbLeftNodes,
                       KOKKOS_LAMBDA(const int& pLeftNodes) {
                         int pId(leftNodes[pLeftNodes]);
                         int pNodes(pId);
                         Vnode_nplus1(pNodes) = nodeVelocityBoundaryCondition(
                             options->leftBC, options->leftBCValue,
                             Mnode(pNodes), G(pNodes));
                       });
  auto rightNodes(mesh->getRightNodes());
  Kokkos::parallel_for("computeBoundaryNodeVelocities", nbRightNodes,
                       KOKKOS_LAMBDA(const int& pRightNodes) {
                         int pId(rightNodes[pRightNodes]);
                         int pNodes(pId);
                         Vnode_nplus1(pNodes) = nodeVelocityBoundaryCondition(
                             options->rightBC, options->rightBCValue,
                             Mnode(pNodes), G(pNodes));
                       });
  auto topNodes(mesh->getTopNodes());
  Kokkos::parallel_for("computeBoundaryNodeVelocities", nbTopNodes,
                       KOKKOS_LAMBDA(const int& pTopNodes) {
                         int pId(topNodes[pTopNodes]);
                         int pNodes(pId);
                         Vnode_nplus1(pNodes) = nodeVelocityBoundaryCondition(
                             options->topBC, options->topBCValue, Mnode(pNodes),
                             G(pNodes));
                       });
  auto bottomNodes(mesh->getBottomNodes());
  Kokkos::parallel_for("computeBoundaryNodeVelocities", nbBottomNodes,
                       KOKKOS_LAMBDA(const int& pBottomNodes) {
                         int pId(bottomNodes[pBottomNodes]);
                         int pNodes(pId);
                         Vnode_nplus1(pNodes) = nodeVelocityBoundaryCondition(
                             options->bottomBC, options->bottomBCValue,
                             Mnode(pNodes), G(pNodes));
                       });
  auto topLeftNode(mesh->getTopLeftNode());
  Kokkos::parallel_for(
      "computeBoundaryNodeVelocities", nbTopLeftNode,
      KOKKOS_LAMBDA(const int& pTopLeftNode) {
        int pId(topLeftNode[pTopLeftNode]);
        int pNodes(pId);
        Vnode_nplus1(pNodes) = nodeVelocityBoundaryConditionCorner(
            options->topBC, options->topBCValue, options->leftBC,
            options->leftBCValue, Mnode(pNodes), G(pNodes));
      });
  auto topRightNode(mesh->getTopRightNode());
  Kokkos::parallel_for(
      "computeBoundaryNodeVelocities", nbTopRightNode,
      KOKKOS_LAMBDA(const int& pTopRightNode) {
        int pId(topRightNode[pTopRightNode]);
        int pNodes(pId);
        Vnode_nplus1(pNodes) = nodeVelocityBoundaryConditionCorner(
            options->topBC, options->topBCValue, options->rightBC,
            options->rightBCValue, Mnode(pNodes), G(pNodes));
      });
  auto bottomLeftNode(mesh->getBottomLeftNode());
  Kokkos::parallel_for(
      "computeBoundaryNodeVelocities", nbBottomLeftNode,
      KOKKOS_LAMBDA(const int& pBottomLeftNode) {
        int pId(bottomLeftNode[pBottomLeftNode]);
        int pNodes(pId);
        Vnode_nplus1(pNodes) = nodeVelocityBoundaryConditionCorner(
            options->bottomBC, options->bottomBCValue, options->leftBC,
            options->leftBCValue, Mnode(pNodes), G(pNodes));
      });
  auto bottomRightNode(mesh->getBottomRightNode());
  Kokkos::parallel_for(
      "computeBoundaryNodeVelocities", nbBottomRightNode,
      KOKKOS_LAMBDA(const int& pBottomRightNode) {
        int pId(bottomRightNode[pBottomRightNode]);
        int pNodes(pId);
        Vnode_nplus1(pNodes) = nodeVelocityBoundaryConditionCorner(
            options->bottomBC, options->bottomBCValue, options->rightBC,
            options->rightBCValue, Mnode(pNodes), G(pNodes));
      });
}
KOKKOS_INLINE_FUNCTION
RealArray1D<EucclhydRemap::dim> EucclhydRemap::nodeVelocityBoundaryCondition(
    int BC, RealArray1D<EucclhydRemap::dim> BCValue,
    RealArray2D<EucclhydRemap::dim, EucclhydRemap::dim> Mp,
    RealArray1D<EucclhydRemap::dim> Gp) {
  if (BC == 200)
    return ArrayOperations::multiply(
        MathFunctions::dot(Gp, BCValue) /
            (MathFunctions::dot(MathFunctions::matVectProduct(Mp, BCValue),
                                BCValue)),
        BCValue);
  else if (BC == 201)
    return BCValue;
  else if (BC == 202)
    return MathFunctions::matVectProduct(inverse(Mp), Gp);

  return options
      ->zeroVect;  // inutile juste pour eviter le warning de compilation
}

KOKKOS_INLINE_FUNCTION
RealArray1D<EucclhydRemap::dim>
EucclhydRemap::nodeVelocityBoundaryConditionCorner(
    int BC1, RealArray1D<EucclhydRemap::dim> BCValue1, int BC2,
    RealArray1D<EucclhydRemap::dim> BCValue2,
    RealArray2D<EucclhydRemap::dim, EucclhydRemap::dim> Mp,
    RealArray1D<EucclhydRemap::dim> Gp) {
  if (BC1 == 200 && BC2 == 200) {
    if (MathFunctions::fabs(
            MathFunctions::fabs(MathFunctions::dot(BCValue1, BCValue2)) -
            MathFunctions::norm(BCValue1) * MathFunctions::norm(BCValue2)) <
        1.0E-8)
      return ArrayOperations::multiply(
          MathFunctions::dot(Gp, BCValue1) /
              (MathFunctions::dot(MathFunctions::matVectProduct(Mp, BCValue1),
                                  BCValue1)),
          BCValue1);
    else {
      return options->zeroVect;
    }
  } else if (BC1 == 201 && BC2 == 201) {
    return ArrayOperations::multiply(
        0.5, (ArrayOperations::plus(BCValue1, BCValue2)));
  } else if (BC1 == 202 && BC2 == 202) {
    return MathFunctions::matVectProduct(inverse(Mp), Gp);
  } else {
    { return options->zeroVect; }
  }
}

RealArray1D<EucclhydRemap::nbequamax> EucclhydRemap::computeBoundaryFluxes(
    int proj, int cCells, RealArray1D<EucclhydRemap::dim> exy) {
  RealArray1D<nbequamax> phiFace_fFaces = options->Uzero;
  int nbCellX = options->X_EDGE_ELEMS;
  int nbCellY = options->Y_EDGE_ELEMS;
  if (options->bottomFluxBC == 1 && cCells < nbCellX && exy[1] == 1) {
    // cellules Bottom
    int cId(cCells);
    int fbBottomFaceOfCellC(mesh->getBottomFaceOfCell(cId));
    int fbId(fbBottomFaceOfCellC);
    int fbFaces(utils::indexOf(mesh->getFaces(), fbId));
    int fbFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), fbId));

    int ftTopFaceOfCellC(mesh->getTopFaceOfCell(cId));
    int ftId(ftTopFaceOfCellC);
    int ftFaces(utils::indexOf(mesh->getFaces(), ftId));
    int ftFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), ftId));

    std::cout << " Bottom cell   " << cCells << std::endl;
    if (proj == 1) phiFace_fFaces = phiFace1(ftFaces);
    if (proj == 2) phiFace_fFaces = phiFace2(ftFaces);
    return computeRemapFlux(
        options->projectionAvecPlateauPente, faceNormalVelocity(fbFaces),
        faceNormal(fbFaces), faceLength(fbFaces), phiFace_fFaces,
        outerFaceNormal(cCells, fbFacesOfCellC), exy, deltat_n);
  }
  if (options->topFluxBC == 1 && cCells <= nbCellX * (nbCellY - 1) &&
      cCells < nbCellX * nbCellY && exy[1] == 1) {
    // cellules top
    int cId(cCells);
    int fbBottomFaceOfCellC(mesh->getBottomFaceOfCell(cId));
    int fbId(fbBottomFaceOfCellC);
    int fbFaces(utils::indexOf(mesh->getFaces(), fbId));
    int fbFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), fbId));

    int ftTopFaceOfCellC(mesh->getTopFaceOfCell(cId));
    int ftId(ftTopFaceOfCellC);
    int ftFaces(utils::indexOf(mesh->getFaces(), ftId));
    int ftFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), ftId));
    std::cout << " Top cell   " << cCells << std::endl;

    if (proj == 1) phiFace_fFaces = phiFace1(fbFaces);
    if (proj == 2) phiFace_fFaces = phiFace2(fbFaces);
    return computeRemapFlux(
        options->projectionAvecPlateauPente, faceNormalVelocity(ftFaces),
        faceNormal(ftFaces), faceLength(ftFaces), phiFace_fFaces,
        outerFaceNormal(cCells, ftFacesOfCellC), exy, deltat_n);
  }
  if (options->leftFluxBC == 1 && exy[0] == 1) {
    // cellules de gauche - a optimiser
    for (int icCells = 0; icCells < nbCellX * nbCellY;
         icCells = icCells + nbCellX) {
      if (icCells == cCells) {
        int cId(cCells);
        int frRightFaceOfCellC(mesh->getRightFaceOfCell(cId));
        int frId(frRightFaceOfCellC);
        int frFaces(utils::indexOf(mesh->getFaces(), frId));
        int frFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), frId));

        int flLeftFaceOfCellC(mesh->getLeftFaceOfCell(cId));
        int flId(flLeftFaceOfCellC);
        int flFaces(utils::indexOf(mesh->getFaces(), flId));
        int flFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), flId));

        std::cout << " Left cell   " << cCells << std::endl;

        if (proj == 1) phiFace_fFaces = phiFace1(frFaces);
        if (proj == 2) phiFace_fFaces = phiFace2(frFaces);
        return computeRemapFlux(
            options->projectionAvecPlateauPente, faceNormalVelocity(flFaces),
            faceNormal(flFaces), faceLength(flFaces), phiFace_fFaces,
            outerFaceNormal(cCells, flFacesOfCellC), exy, deltat_n);
      }
    }
  }
  if (options->rightFluxBC == 1 && exy[0] == 1) {
    // cellules de droite
    for (int icCells = nbCellX - 1; icCells < nbCellX * nbCellY;
         icCells = icCells + nbCellX) {
      if (icCells == cCells) {
        int cId(cCells);
        int frRightFaceOfCellC(mesh->getRightFaceOfCell(cId));
        int frId(frRightFaceOfCellC);
        int frFaces(utils::indexOf(mesh->getFaces(), frId));
        int frFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), frId));

        int flLeftFaceOfCellC(mesh->getLeftFaceOfCell(cId));
        int flId(flLeftFaceOfCellC);
        int flFaces(utils::indexOf(mesh->getFaces(), flId));
        int flFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), flId));

        if (proj == 1) phiFace_fFaces = phiFace1(flFaces);
        if (proj == 2) phiFace_fFaces = phiFace2(flFaces);

        if (proj == 1 &&
            (cCells == dbgcell1 || cCells == dbgcell2 || cCells == dbgcell3)) {
          std::cout << " AP 1 Right cell   " << exy << " " << cCells << " "
                    << flFaces << " " << phiFace1(flFaces) << " " << frFaces
                    << " " << phiFace_fFaces << std::endl;
          std::cout << " faceNormalVelocity " << faceNormalVelocity(flFaces)
                    << " " << faceNormalVelocity(frFaces) << std::endl;
          std::cout << " faceLength " << faceLength(flFaces) << " "
                    << faceLength(frFaces) << std::endl;
          std::cout << " outerFaceNormal(cCells,frFacesOfCellC) "
                    << outerFaceNormal(cCells, frFacesOfCellC)
                    << " frFacesOfCellC " << frFacesOfCellC << std::endl;
        }

        if (proj == 2 &&
            (cCells == dbgcell1 || cCells == dbgcell2 || cCells == dbgcell3)) {
          std::cout << " AP 2 Right cell   " << exy << " " << cCells << " "
                    << flFaces << " " << phiFace2(flFaces) << " " << frFaces
                    << " " << phiFace_fFaces << std::endl;
          std::cout << " faceNormalVelocity " << faceNormalVelocity(flFaces)
                    << " " << faceNormalVelocity(frFaces) << std::endl;
          std::cout << " faceLength " << faceLength(flFaces) << " "
                    << faceLength(frFaces) << std::endl;
          std::cout << " outerFaceNormal(cCells,frFacesOfCellC) "
                    << outerFaceNormal(cCells, frFacesOfCellC)
                    << " frFacesOfCellC " << frFacesOfCellC << std::endl;
        }
        //
        return computeRemapFlux(
            options->projectionAvecPlateauPente, faceNormalVelocity(frFaces),
            faceNormal(frFaces), faceLength(frFaces), phiFace_fFaces,
            outerFaceNormal(cCells, frFacesOfCellC), exy, deltat_n);
      }
    }
  }
  return phiFace_fFaces;  // options->Uzero;
}
