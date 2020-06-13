#include <math.h>  // for sqrt
#include <Kokkos_Core.hpp>
#include <algorithm>                // for copy
#include <array>                    // for array
#include <iostream>                 // for operator<<, basic_ostream::operat...
#include <vector>                   // for allocator, vector
#include "EucclhydRemap.h"          // for EucclhydRemap, EucclhydRemap::Opt...
#include "UtilesRemap-Impl.h"       // for EucclhydRemap::computeIntersectionPP
#include "mesh/CartesianMesh2D.h"   // for CartesianMesh2D
#include "types/ArrayOperations.h"  // for divide, minus, plus
#include "types/MathFunctions.h"    // for dot
#include "types/MultiArray.h"       // for operator<<
#include "utils/Utils.h"            // for indexOf

/**
 * Job computeGradPhiFace2 called @12.0 in executeTimeLoopN method.
 * In variables: Uremap1, deltaxLagrange, projectionOrder, vLagrange, x_then_y_n
 * Out variables: gradPhiFace2
 */
void EucclhydRemap::computeGradPhiFace2() noexcept {
  if (options->projectionOrder > 1) {
    if (x_then_y_n) {
      auto innerHorizontalFaces(mesh->getInnerHorizontalFaces());
      Kokkos::parallel_for(
          "computeGradPhiFace2", nbInnerHorizontalFaces,
          KOKKOS_LAMBDA(const int& fInnerHorizontalFaces) {
            int fId(innerHorizontalFaces[fInnerHorizontalFaces]);
            int fFaces(utils::indexOf(mesh->getFaces(), fId));
            int cfFrontCellF(mesh->getFrontCell(fId));
            int cfId(cfFrontCellF);
            int cfCells(cfId);
            int cbBackCellF(mesh->getBackCell(fId));
            int cbId(cbBackCellF);
            int cbCells(cbId);
            // gradPhiFace2(fFaces) =
            // ArrayOperations::divide((ArrayOperations::minus(ArrayOperations::divide(Uremap1(cfCells),
            // vLagrange(cfCells)),
            // ArrayOperations::divide(Uremap1(cbCells),
            // vLagrange(cbCells)))),
            // deltaxLagrange(fFaces));
            gradPhiFace2(fFaces) = ArrayOperations::divide(
                ArrayOperations::minus(Phi(cfCells), Phi(cbCells)),
                deltaxLagrange(fFaces));
            // std::cout << "gradphiface2 H " << gradPhiFace2(fFaces) <<
            // std::endl; Phi(cfCells) =
            // ArrayOperations::divide(Uremap1(cfCells), vLagrange(cfCells));
            // Phi(cbCells) = ArrayOperations::divide(Uremap1(cbCells),
            // vLagrange(cbCells));
            int n1FirstNodeOfFaceF(mesh->getFirstNodeOfFace(fId));
            int n1Id(n1FirstNodeOfFaceF);
            int n1Nodes(n1Id);
            int n2SecondNodeOfFaceF(mesh->getSecondNodeOfFace(fId));
            int n2Id(n2SecondNodeOfFaceF);
            int n2Nodes(n2Id);
            LfLagrange(fFaces) =
                sqrt((XLagrange(n1Nodes)[0] - XLagrange(n2Nodes)[0]) *
                         (XLagrange(n1Nodes)[0] - XLagrange(n2Nodes)[0]) +
                     (XLagrange(n1Nodes)[1] - XLagrange(n2Nodes)[1]) *
                         (XLagrange(n1Nodes)[1] - XLagrange(n2Nodes)[1]));
            HvLagrange(cfCells) = vLagrange(cfCells) / LfLagrange(fFaces);
            HvLagrange(cbCells) = vLagrange(cbCells) / LfLagrange(fFaces);
            // seconde methode
            HvLagrange(cfCells) = 0.;
            HvLagrange(cbCells) = 0.;
            int cbfrRightFaceOfCellC(mesh->getRightFaceOfCell(cbId));
            int cbfrId(cbfrRightFaceOfCellC);
            int cbfrFaces(utils::indexOf(mesh->getFaces(), cbfrId));

            HvLagrange(cbCells) += 0.5 * faceLengthLagrange(cbfrFaces);
            int cbflLeftFaceOfCellC(mesh->getLeftFaceOfCell(cbId));
            int cbflId(cbflLeftFaceOfCellC);
            int cbflFaces(utils::indexOf(mesh->getFaces(), cbflId));
            HvLagrange(cbCells) += 0.5 * faceLengthLagrange(cbflFaces);

            int frRightFaceOfCellC(mesh->getRightFaceOfCell(cfId));
            int frId(frRightFaceOfCellC);
            int frFaces(utils::indexOf(mesh->getFaces(), frId));
            HvLagrange(cfCells) += 0.5 * faceLengthLagrange(frFaces);

            int flLeftFaceOfCellC(mesh->getLeftFaceOfCell(cbId));
            int flId(flLeftFaceOfCellC);
            int flFaces(utils::indexOf(mesh->getFaces(), flId));
            HvLagrange(cfCells) += 0.5 * faceLengthLagrange(flFaces);

            // // troisieme methode (longueur constante - face de droite par
            // exemple) HvLagrange(cfCells) = 0.; HvLagrange(cbCells) = 0.; int
            // cbfrRightFaceOfCellC(mesh->getRightFaceOfCell(cbId)); int
            // cbfrId(cbfrRightFaceOfCellC); int
            // cbfrFaces(utils::indexOf(mesh->getFaces(),cbfrId));
            // HvLagrange(cbCells) = faceLength(cbfrFaces);
            // int frRightFaceOfCellC(mesh->getRightFaceOfCell(cfId));
            // int frId(frRightFaceOfCellC);
            // int frFaces(utils::indexOf(mesh->getFaces(),frId));
            // //std::cout << " FaceLengthLagrange " <<
            // faceLengthLagrange(frFaces) << std::endl; HvLagrange(cfCells) =
            // faceLength(frFaces);

            if (fFaces == face_debug1 || fFaces == face_debug2) {
              std::cout << " Pour Proj Verti fFaces " << fFaces << " Longeur L"
                        << LfLagrange(fFaces) << std::endl;
              std::cout << " FaceLengthLagrange " << faceLengthLagrange(fFaces)
                        << std::endl;
              std::cout << " HvLagrange(cfCells) " << frFaces << " "
                        << HvLagrange(cfCells) << std::endl;
              std::cout << " HvLagrange(cbCells) " << cbfrFaces << " "
                        << HvLagrange(cbCells) << std::endl;
            }
          });
    } else {
      auto innerVerticalFaces(mesh->getInnerVerticalFaces());
      Kokkos::parallel_for(
          "computeGradPhiFace2", nbInnerVerticalFaces,
          KOKKOS_LAMBDA(const int& fInnerVerticalFaces) {
            int fId(innerVerticalFaces[fInnerVerticalFaces]);
            int fFaces(utils::indexOf(mesh->getFaces(), fId));
            int cfFrontCellF(mesh->getFrontCell(fId));
            int cfId(cfFrontCellF);
            int cfCells(cfId);
            int cbBackCellF(mesh->getBackCell(fId));
            int cbId(cbBackCellF);
            int cbCells(cbId);
            // gradPhiFace2(fFaces) =
            // ArrayOperations::divide((ArrayOperations::minus(ArrayOperations::divide(Uremap1(cfCells),
            // vLagrange(cfCells)),
            // ArrayOperations::divide(Uremap1(cbCells),
            // vLagrange(cbCells)))),
            // deltaxLagrange(fFaces));
            gradPhiFace2(fFaces) = ArrayOperations::divide(
                ArrayOperations::minus(Phi(cfCells), Phi(cbCells)),
                deltaxLagrange(fFaces));
            // std::cout << "gradphiface2 V" << gradPhiFace2(fFaces) <<
            // std::endl; Phi(cfCells) =
            // ArrayOperations::divide(Uremap1(cfCells), vLagrange(cfCells));
            // Phi(cbCells) = ArrayOperations::divide(Uremap1(cbCells),
            // vLagrange(cbCells));
            int n1FirstNodeOfFaceF(mesh->getFirstNodeOfFace(fId));
            int n1Id(n1FirstNodeOfFaceF);
            int n1Nodes(n1Id);
            int n2SecondNodeOfFaceF(mesh->getSecondNodeOfFace(fId));
            int n2Id(n2SecondNodeOfFaceF);
            int n2Nodes(n2Id);
            LfLagrange(fFaces) =
                sqrt((XLagrange(n1Nodes)[0] - XLagrange(n2Nodes)[0]) *
                         (XLagrange(n1Nodes)[0] - XLagrange(n2Nodes)[0]) +
                     (XLagrange(n1Nodes)[1] - XLagrange(n2Nodes)[1]) *
                         (XLagrange(n1Nodes)[1] - XLagrange(n2Nodes)[1]));
            HvLagrange(cfCells) = vLagrange(cfCells) / LfLagrange(fFaces);
            HvLagrange(cbCells) = vLagrange(cbCells) / LfLagrange(fFaces);

            // seconde methode
            HvLagrange(cfCells) = 0.;
            HvLagrange(cbCells) = 0.;
            int cbfbBottomFaceOfCellC(mesh->getBottomFaceOfCell(cbId));
            int cbfbId(cbfbBottomFaceOfCellC);
            int cbfbFaces(utils::indexOf(mesh->getFaces(), cbfbId));

            HvLagrange(cbCells) += 0.5 * faceLengthLagrange(cbfbFaces);
            int cbftTopFaceOfCellC(mesh->getTopFaceOfCell(cbId));
            int cbftId(cbftTopFaceOfCellC);
            int cbftFaces(utils::indexOf(mesh->getFaces(), cbftId));
            HvLagrange(cbCells) += 0.5 * faceLengthLagrange(cbftFaces);

            int fbBottomFaceOfCellC(mesh->getBottomFaceOfCell(cfId));
            int fbId(fbBottomFaceOfCellC);
            int fbFaces(utils::indexOf(mesh->getFaces(), fbId));
            HvLagrange(cfCells) += 0.5 * faceLengthLagrange(fbFaces);

            int ftTopFaceOfCellC(mesh->getTopFaceOfCell(cfId));
            int ftId(ftTopFaceOfCellC);
            int ftFaces(utils::indexOf(mesh->getFaces(), ftId));
            HvLagrange(cfCells) += 0.5 * faceLengthLagrange(ftFaces);

            // // troisieme methode (longueur constante - face du bas par
            // exemple) HvLagrange(cfCells) = 0.; HvLagrange(cbCells) = 0.; int
            // cbfbBottomFaceOfCellC(mesh->getBottomFaceOfCell(cbId)); int
            // cbfbId(cbfbBottomFaceOfCellC); int
            // cbfbFaces(utils::indexOf(mesh->getFaces(),cbfbId));
            // HvLagrange(cbCells) = faceLength(cbfbFaces);
            // int fbBottomFaceOfCellC(mesh->getBottomFaceOfCell(cfId));
            // int fbId(fbBottomFaceOfCellC);
            // int fbFaces(utils::indexOf(mesh->getFaces(),fbId));
            // HvLagrange(cfCells) = faceLength(fbFaces);

            if (fFaces == face_debug1 || fFaces == face_debug2) {
              std::cout << " Pour Proj Horizontal fFaces " << fFaces
                        << " Longeur L" << LfLagrange(fFaces) << std::endl;
              std::cout << " FaceLengthLagrange " << faceLengthLagrange(fFaces)
                        << std::endl;
              std::cout << " HvLagrange(cfCells) " << fbFaces << " "
                        << HvLagrange(cfCells) << std::endl;
              std::cout << " HvLagrange(cbCells) " << cbfbFaces << " "
                        << HvLagrange(cbCells) << std::endl;
            }
          });
    }
  }
}

/**
 * Job computeGradPhi2 called @13.0 in executeTimeLoopN method.
 * In variables: gradPhiFace2, projectionLimiterId, projectionOrder, x_then_y_n
 * Out variables: gradPhi2
 */
void EucclhydRemap::computeGradPhi2() noexcept {
  if (options->projectionOrder > 1) {
    if (x_then_y_n) {
      // std::cout << " Phase 2 Verticale computeGradPhi2 " << std::endl;
      Kokkos::parallel_for(
          "computeGradPhi2", nbCells, KOKKOS_LAMBDA(const int& cCells) {
            int cId(cCells);
            int fbBottomFaceOfCellC(mesh->getBottomFaceOfCell(cId));
            int fbId(fbBottomFaceOfCellC);
            int fbFaces(utils::indexOf(mesh->getFaces(), fbId));
            int ftTopFaceOfCellC(mesh->getTopFaceOfCell(cId));
            int ftId(ftTopFaceOfCellC);
            int ftFaces(utils::indexOf(mesh->getFaces(), ftId));
            // maille dessous
            int cfFrontCellF(mesh->getFrontCell(fbId));
            int cfId(cfFrontCellF);
            int cfCells(cfId);
            int fbFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), fbId));
            // maille dessus
            int cbBackCellF(mesh->getBackCell(ftId));
            int cbId(cbBackCellF);
            int cbCells(cbId);
            int ftFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), ftId));

            if (cbCells == -1) cbCells = cCells;
            if (cfCells == -1) cfCells = cCells;
            bool voisinage_pure =
                (options->projectionLimiteurMixte == 1) &&
                (mixte(cCells) == 0 && mixte(cfCells) == 0 &&
                 mixte(cbCells) == 0 && pure(cCells) == pure(cfCells) &&
                 pure(cCells) == pure(cbCells));

            int limiter = options->projectionLimiterId;
            if ((options->projectionAvecPlateauPente == 1) && voisinage_pure)
              limiter = options->projectionLimiterIdPure;
            // if (cCells == 4005 || cCells == 3864) std::cout << cCells << "
            // Phase 2 V voisinage pure =  " << voisinage_pure << std::endl;
            gradPhi2(cCells) = computeAndLimitGradPhi(
                limiter, gradPhiFace2(fbFaces), gradPhiFace2(ftFaces),
                Phi(cCells), Phi(cbCells), Phi(cfCells), HvLagrange(cCells),
                HvLagrange(cbCells), HvLagrange(cfCells));
            if (cCells == dbgcell3 || cCells == dbgcell2 ||
                cCells == dbgcell1) {
              std::cout << " Verti cell   " << cCells << "gradphi1 "
                        << gradPhi2(cCells) << std::endl;
              std::cout << " Grad Phi  " << gradPhiFace1(fbFaces) << " "
                        << gradPhiFace1(ftFaces) << std::endl;
            }

            if (options->projectionAvecPlateauPente == 1) {
              // std::cout << "gradphiface2 " << gradPhiFace2(fbFaces) <<
              // gradPhiFace2(ftFaces) << std::endl; std::cout << "gradphi2 "
              // <<gradPhi2(cCells) << std::endl;

              // std::cout << " Phase 2 Verticale " << std::endl;
              RealArray1D<dim> exy = {{0.0, 1.0}};
              // std::cout << " cfCells " << cfCells << " cCells " << cCells <<
              // " cbCells " << cbCells << std::endl; std::cout << " fbFaces "
              // << fbFaces << " fbId" << fbId << std::endl;

              double Flux_sortant_av =
                  MathFunctions::dot(outerFaceNormal(cCells, fbFacesOfCellC),
                                     exy) *
                  faceNormalVelocity(fbFaces);
              if (voisinage_pure)
                deltaPhiFaceAr(cCells) = computeIntersectionPPPure(
                    gradPhi2(cCells), Phi(cCells), Phi(cbCells), Phi(cfCells),
                    HvLagrange(cCells), HvLagrange(cbCells),
                    HvLagrange(cfCells), Flux_sortant_av, deltat_n, 0, cCells,
                    options->threshold);
              else
                deltaPhiFaceAr(cCells) = computeIntersectionPP(
                    gradPhi2(cCells), Phi(cCells), Phi(cbCells), Phi(cfCells),
                    HvLagrange(cCells), HvLagrange(cbCells),
                    HvLagrange(cfCells), Flux_sortant_av, deltat_n, 0, cCells,
                    options->threshold);

              // std::cout << " fbFaces " << fbFaces << " deltaphiface " <<
              // deltaPhiFaceAr(cCells) <<std::endl; std::cout << " ftFaces " <<
              // ftFaces << " ftId" << ftId << std::endl;
              double Flux_sortant_ar =
                  MathFunctions::dot(outerFaceNormal(cCells, ftFacesOfCellC),
                                     exy) *
                  faceNormalVelocity(ftFaces);
              if (voisinage_pure)
                deltaPhiFaceAv(cCells) = computeIntersectionPPPure(
                    gradPhi2(cCells), Phi(cCells), Phi(cbCells), Phi(cfCells),
                    HvLagrange(cCells), HvLagrange(cbCells),
                    HvLagrange(cfCells), Flux_sortant_ar, deltat_n, 1, cCells,
                    options->threshold);
              else
                deltaPhiFaceAv(cCells) = computeIntersectionPP(
                    gradPhi2(cCells), Phi(cCells), Phi(cbCells), Phi(cfCells),
                    HvLagrange(cCells), HvLagrange(cbCells),
                    HvLagrange(cfCells), Flux_sortant_ar, deltat_n, 1, cCells,
                    options->threshold);

              // std::cout << " ftFaces " << ftFaces << " deltaphiface " <<
              // deltaPhiFaceAv(cCells) <<std::endl;
            }
          });
    } else {
      // std::cout << " Phase 2 Horizontale computeGradPhi2 " << std::endl;
      Kokkos::parallel_for(
          "computeGradPhi2", nbCells, KOKKOS_LAMBDA(const int& cCells) {
            int cId(cCells);
            int frRightFaceOfCellC(mesh->getRightFaceOfCell(cId));
            int frId(frRightFaceOfCellC);
            int frFaces(utils::indexOf(mesh->getFaces(), frId));
            int flLeftFaceOfCellC(mesh->getLeftFaceOfCell(cId));
            int flId(flLeftFaceOfCellC);
            int flFaces(utils::indexOf(mesh->getFaces(), flId));

            // maille devant
            int cfFrontCellF(mesh->getFrontCell(frId));
            int cfId(cfFrontCellF);
            int cfCells(cfId);
            int frFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), frId));
            // maille deriere
            int cbBackCellF(mesh->getBackCell(flId));
            int cbId(cbBackCellF);
            int cbCells(cbId);
            int flFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), flId));

            if (cbCells == -1) cbCells = cCells;
            if (cfCells == -1) cfCells = cCells;
            bool voisinage_pure =
                (options->projectionLimiteurMixte == 1) &&
                (mixte(cCells) == 0 && mixte(cfCells) == 0 &&
                 mixte(cbCells) == 0 && pure(cCells) == pure(cfCells) &&
                 pure(cCells) == pure(cbCells));

            int limiter = options->projectionLimiterId;
            if ((options->projectionAvecPlateauPente == 1) && voisinage_pure)
              limiter = options->projectionLimiterIdPure;
            // if (cCells == 4005 || cCells == 3864) std::cout << cCells << "
            // Phase 2 H voisinage pure =  " << voisinage_pure << std::endl;
            gradPhi2(cCells) = computeAndLimitGradPhi(
                limiter, gradPhiFace2(frFaces), gradPhiFace2(flFaces),
                Phi(cCells), Phi(cfCells), Phi(cbCells), HvLagrange(cCells),
                HvLagrange(cfCells), HvLagrange(cbCells));
            if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                 cCells == dbgcell1) &&
                test_debug == 1) {
              std::cout << " Hori cell   " << cCells << "gradphi2 "
                        << gradPhi2(cCells) << std::endl;
              std::cout << " Grad Phi  " << gradPhiFace1(flFaces) << " "
                        << gradPhiFace1(frFaces) << std::endl;
            }
            if (options->projectionAvecPlateauPente == 1) {
              // std::cout << " Phase 2 Horizontale " << std::endl;
              RealArray1D<dim> exy = {{1.0, 0.0}};
              // std::cout << " bCells " << cbCells << " Cells " << cCells << "
              // fCells " << cfCells << std::endl; std::cout << " flFaces " <<
              // flFaces << " flId" << flId <<std::endl;
              if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                   cCells == dbgcell1) &&
                  test_debug == 1) {
                std::cout << " cell   " << cCells << std::endl;
              }
              double Flux_sortant_ar =
                  MathFunctions::dot(outerFaceNormal(cCells, flFacesOfCellC),
                                     exy) *
                  faceNormalVelocity(flFaces);
              if (voisinage_pure)
                deltaPhiFaceAr(cCells) = computeIntersectionPPPure(
                    gradPhi2(cCells), Phi(cCells), Phi(cfCells), Phi(cbCells),
                    HvLagrange(cCells), HvLagrange(cfCells),
                    HvLagrange(cbCells), Flux_sortant_ar, deltat_n, 0, cCells,
                    options->threshold);
              else
                deltaPhiFaceAr(cCells) = computeIntersectionPP(
                    gradPhi2(cCells), Phi(cCells), Phi(cfCells), Phi(cbCells),
                    HvLagrange(cCells), HvLagrange(cfCells),
                    HvLagrange(cbCells), Flux_sortant_ar, deltat_n, 0, cCells,
                    options->threshold);
              if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                   cCells == dbgcell1) &&
                  test_debug == 1) {
                std::cout << " flFaces " << flFaces << " deltaphiface "
                          << deltaPhiFaceAr(cCells) << std::endl;
                std::cout << " frFaces " << frFaces << " frId" << frId
                          << std::endl;
              }
              double Flux_sortant_av =
                  MathFunctions::dot(outerFaceNormal(cCells, frFacesOfCellC),
                                     exy) *
                  faceNormalVelocity(frFaces);
              if (voisinage_pure)
                deltaPhiFaceAv(cCells) = computeIntersectionPPPure(
                    gradPhi2(cCells), Phi(cCells), Phi(cfCells), Phi(cbCells),
                    HvLagrange(cCells), HvLagrange(cfCells),
                    HvLagrange(cbCells), Flux_sortant_av, deltat_n, 1, cCells,
                    options->threshold);
              else
                deltaPhiFaceAv(cCells) = computeIntersectionPP(
                    gradPhi2(cCells), Phi(cCells), Phi(cfCells), Phi(cbCells),
                    HvLagrange(cCells), HvLagrange(cfCells),
                    HvLagrange(cbCells), Flux_sortant_av, deltat_n, 1, cCells,
                    options->threshold);

              if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                   cCells == dbgcell1) &&
                  test_debug == 1) {
                std::cout << " frFaces " << frFaces << " deltaphiface "
                          << deltaPhiFaceAv(cCells) << std::endl;
              }
            }
          });
    }
  }
}

/**
 * Job computeUpwindFaceQuantitiesForProjection2 called @14.0 in
 * executeTimeLoopN method. In variables: Uremap1, XcLagrange, Xf,
 * deltaxLagrange, faceNormal, faceNormalVelocity, gradPhi2, vLagrange,
 * x_then_y_n Out variables: phiFace2
 */
void EucclhydRemap::computeUpwindFaceQuantitiesForProjection2() noexcept {
  if (x_then_y_n) {
    // std::cout << " Phase Projection 2 Verticale " << std::endl;
    auto innerHorizontalFaces(mesh->getInnerHorizontalFaces());
    Kokkos::parallel_for(
        "computeUpwindFaceQuantitiesForProjection2", nbInnerHorizontalFaces,
        KOKKOS_LAMBDA(const int& fInnerHorizontalFaces) {
          int fId(innerHorizontalFaces[fInnerHorizontalFaces]);
          int fFaces(utils::indexOf(mesh->getFaces(), fId));
          int cfFrontCellF(mesh->getFrontCell(fId));
          int cfId(cfFrontCellF);
          int cfCells(cfId);
          int cbBackCellF(mesh->getBackCell(fId));
          int cbId(cbBackCellF);
          int cbCells(cbId);
          if (options->projectionAvecPlateauPente == 0) {
            phiFace2(fFaces) = computeUpwindFaceQuantities(
                faceNormal(fFaces), faceNormalVelocity(fFaces),
                deltaxLagrange(fFaces), Xf(fFaces),
                ArrayOperations::divide(Uremap1(cbCells), vLagrange(cbCells)),
                gradPhi2(cbCells), XcLagrange(cbCells),
                ArrayOperations::divide(Uremap1(cfCells), vLagrange(cfCells)),
                gradPhi2(cfCells), XcLagrange(cfCells));
          } else {
            phiFace2(fFaces) = ArrayOperations::minus(deltaPhiFaceAv(cfCells),
                                                      deltaPhiFaceAr(cbCells));
            // cfCells maille arriere, cbCells maille devant
            // front et back inversé ?
            // std::cout << " Phase Projection 2 Verticale " << std::endl;
            // std::cout << " cbCells " << cbCells << " cfCells " << cfCells <<
            // std::endl; std::cout << " fFaces " << fFaces << " fId" << fId <<
            // std::endl; std::cout << " phiFace2(fFaces) " << phiFace2(fFaces)
            // <<  std::endl; std::cout << " deltaPhiFaceAv(cfCells) " <<
            // deltaPhiFaceAv(cfCells) <<  std::endl; std::cout << "
            // deltaPhiFaceAr(cbCells) " << deltaPhiFaceAr(cbCells) <<
            // std::endl;
          }
          if (fFaces == face_debug1 || fFaces == face_debug2) {
            std::cout << " Phase Projection 2 Verticale " << std::endl;
            std::cout << " cbCells " << cbCells << " cfCells " << cfCells
                      << std::endl;
            std::cout << " fFaces " << fFaces << " fId" << fId << std::endl;
            std::cout << " phiFace2(fFaces) " << phiFace2(fFaces) << std::endl;
            std::cout << " deltaPhiFaceAv(cfCells) " << deltaPhiFaceAv(cfCells)
                      << std::endl;
            std::cout << " deltaPhiFaceAr(cbCells) " << deltaPhiFaceAr(cbCells)
                      << std::endl;
          }
        });
  } else {
    // std::cout << " Phase Projection 2 Horizontale " << std::endl;
    auto innerVerticalFaces(mesh->getInnerVerticalFaces());
    Kokkos::parallel_for(
        "computeUpwindFaceQuantitiesForProjection2", nbInnerVerticalFaces,
        KOKKOS_LAMBDA(const int& fInnerVerticalFaces) {
          int fId(innerVerticalFaces[fInnerVerticalFaces]);
          int fFaces(utils::indexOf(mesh->getFaces(), fId));
          int cfFrontCellF(mesh->getFrontCell(fId));
          int cfId(cfFrontCellF);
          int cfCells(cfId);
          int cbBackCellF(mesh->getBackCell(fId));
          int cbId(cbBackCellF);
          int cbCells(cbId);
          if (options->projectionAvecPlateauPente == 0) {
            phiFace2(fFaces) = computeUpwindFaceQuantities(
                faceNormal(fFaces), faceNormalVelocity(fFaces),
                deltaxLagrange(fFaces), Xf(fFaces),
                ArrayOperations::divide(Uremap1(cbCells), vLagrange(cbCells)),
                gradPhi2(cbCells), XcLagrange(cbCells),
                ArrayOperations::divide(Uremap1(cfCells), vLagrange(cfCells)),
                gradPhi2(cfCells), XcLagrange(cfCells));
          } else {
            phiFace2(fFaces) = ArrayOperations::minus(deltaPhiFaceAv(cbCells),
                                                      deltaPhiFaceAr(cfCells));
            // std::cout << " Phase Projection 2 Horizontale " << std::endl;
            // std::cout << " cbCells " << cbCells << " cfCells " << cfCells <<
            // std::endl; std::cout << " fFaces " << fFaces << " fId" << fId <<
            // std::endl; std::cout << " phiFace2(fFaces) " << phiFace2(fFaces)
            // <<  std::endl; std::cout << " deltaPhiFaceAv(cbCells) " <<
            // deltaPhiFaceAv(cbCells) <<  std::endl; std::cout << "
            // deltaPhiFaceAr(cfCells) " << deltaPhiFaceAr(cfCells) <<
            // std::endl;
          }
          if (fFaces == face_debug1 || fFaces == face_debug2) {
            std::cout << " Phase Projection 2 Horizontale " << std::endl;
            std::cout << " cbCells " << cbCells << " cfCells " << cfCells
                      << std::endl;
            std::cout << " fFaces " << fFaces << " fId" << fId << std::endl;
            std::cout << " phiFace2(fFaces) " << phiFace2(fFaces) << std::endl;
            std::cout << " deltaPhiFaceAv(cbCells) " << deltaPhiFaceAv(cbCells)
                      << std::endl;
            std::cout << " deltaPhiFaceAr(cfCells) " << deltaPhiFaceAr(cfCells)
                      << std::endl;
          }
        });
  }
}

/**
 * Job computeUremap2 called @15.0 in executeTimeLoopN method.
 * In variables: Uremap1, deltat_n, faceLength, faceNormal, faceNormalVelocity,
 * outerFaceNormal, phiFace2, x_then_y_n Out variables: Uremap2
 */
void EucclhydRemap::computeUremap2() noexcept {
  RealArray1D<dim> exy = xThenYToDirection(!x_then_y_n);
  Kokkos::parallel_for(
      "computeUremap2", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        int cId(cCells);
        RealArray1D<nbequamax> reduction9 = options->Uzero;
        {
          auto neighbourCellsC(mesh->getNeighbourCells(cId));
          for (int dNeighbourCellsC = 0;
               dNeighbourCellsC < neighbourCellsC.size(); dNeighbourCellsC++) {
            int dId(neighbourCellsC[dNeighbourCellsC]);
            int fCommonFaceCD(mesh->getCommonFace(cId, dId));
            int fId(fCommonFaceCD);
            int fFaces(utils::indexOf(mesh->getFaces(), fId));
            int fFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), fId));
            reduction9 = ArrayOperations::plus(
                reduction9,
                (computeRemapFlux(
                    options->projectionAvecPlateauPente,
                    faceNormalVelocity(fFaces), faceNormal(fFaces),
                    faceLength(fFaces), phiFace2(fFaces),
                    outerFaceNormal(cCells, fFacesOfCellC), exy, deltat_n)));
            //
            if (cCells == dbgcell3 || cCells == dbgcell2 ||
                cCells == dbgcell1) {
              std::cout << " --------- flux face interieur "
                           "----------------------------------"
                        << std::endl;
              std::cout << " cell   " << cCells << "reduction9 " << reduction9
                        << std::endl;
              std::cout << " face " << fFaces << " de longueur "
                        << faceLength(fFaces) << " "
                        << computeRemapFlux(
                               options->projectionAvecPlateauPente,
                               faceNormalVelocity(fFaces), faceNormal(fFaces),
                               faceLength(fFaces), phiFace2(fFaces),
                               outerFaceNormal(cCells, fFacesOfCellC), exy,
                               deltat_n)
                        << std::endl;
              std::cout << " --------- ----------------------------------"
                        << std::endl;
            }
            //
          }
          if (options->FluxBC > 0) {
            // flux exterieur
            if ((cCells == dbgcell3 || cCells == dbgcell3 ||
                 cCells == dbgcell1)) {
              std::cout << " --------- flux face exterieur "
                           "----------------------------------"
                        << std::endl;
              std::cout << exy << " AV cell   " << cCells << "reduction9 "
                        << reduction9 << std::endl;
            }

            reduction9 = ArrayOperations::plus(
                reduction9, (computeBoundaryFluxes(2, cCells, exy)));

            if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                 cCells == dbgcell1)) {
              std::cout << exy << " cell   " << cCells << "reduction9 "
                        << reduction9 << std::endl;
              std::cout << computeBoundaryFluxes(2, cCells, exy) << std::endl;
              std::cout << " --------- ----------------------------------"
                        << std::endl;
            }
          }
        }

        Uremap2(cCells) = ArrayOperations::minus(Uremap1(cCells), reduction9);
      });
}
