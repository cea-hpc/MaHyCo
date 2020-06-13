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
 * Job computeGradPhiFace1 called @8.0 in executeTimeLoopN method.
 * In variables: ULagrange, deltaxLagrange, projectionOrder, vLagrange,
 * x_then_y_n Out variables: gradPhiFace1
 */
void EucclhydRemap::computeGradPhiFace1() noexcept {
  if (options->projectionOrder > 1) {
    if (x_then_y_n) {
      auto innerVerticalFaces(mesh->getInnerVerticalFaces());
      Kokkos::parallel_for(
          "computeGradPhiFace1", nbInnerVerticalFaces,
          KOKKOS_LAMBDA(const int& fInnerVerticalFaces) {
            int fId(innerVerticalFaces[fInnerVerticalFaces]);
            int fFaces(utils::indexOf(mesh->getFaces(), fId));
            int cfFrontCellF(mesh->getFrontCell(fId));
            int cfId(cfFrontCellF);
            int cfCells(cfId);
            int cbBackCellF(mesh->getBackCell(fId));
            int cbId(cbBackCellF);
            int cbCells(cbId);
            // gradPhiFace1(fFaces) =
            // ArrayOperations::divide((ArrayOperations::minus(ArrayOperations::divide(ULagrange(cfCells),
            // vLagrange(cfCells)),
            // ArrayOperations::divide(ULagrange(cbCells),
            // vLagrange(cbCells)))), deltaxLagrange(fFaces));
            gradPhiFace1(fFaces) = ArrayOperations::divide(
                ArrayOperations::minus(Phi(cfCells), Phi(cbCells)),
                deltaxLagrange(fFaces));
            // Phi(cfCells) = ArrayOperations::divide(ULagrange(cfCells),
            // vLagrange(cfCells)); Phi(cbCells) =
            // ArrayOperations::divide(ULagrange(cbCells), vLagrange(cbCells));
            // std::cout << " cbCells " << cbCells << " Phi(cbCells) " <<
            // Phi(cbCells) << std::endl; std::cout << " cfCells " << cfCells <<
            // " Phi(cfCells) " << Phi(cfCells) << std::endl;
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
            //
            // troisieme methode (longueur constante - face du bas par exemple)
            // HvLagrange(cfCells) = 0.;
            // HvLagrange(cbCells) = 0.;
            // int cbfbBottomFaceOfCellC(mesh->getBottomFaceOfCell(cbId));
            // int cbfbId(cbfbBottomFaceOfCellC);
            // int cbfbFaces(utils::indexOf(mesh->getFaces(),cbfbId));
            // HvLagrange(cbCells) = faceLength(cbfbFaces);
            // int fbBottomFaceOfCellC(mesh->getBottomFaceOfCell(cfId));
            // int fbId(fbBottomFaceOfCellC);
            // int fbFaces(utils::indexOf(mesh->getFaces(),fbId));
            // HvLagrange(cfCells) = faceLength(fbFaces);

            if (fFaces == face_debug1 || fFaces == face_debug2) {
              std::cout << " Pour Proj Hori fFaces " << fFaces << " Longeur L"
                        << LfLagrange(fFaces) << std::endl;
              std::cout << " FaceLengthLagrange " << faceLengthLagrange(fFaces)
                        << std::endl;
              std::cout << " HvLagrange(cfCells) " << fbFaces << " "
                        << HvLagrange(cfCells) << std::endl;
              std::cout << " HvLagrange(cbCells) " << cbfbFaces << " "
                        << HvLagrange(cbCells) << std::endl;
            }
          });
    } else {
      auto innerHorizontalFaces(mesh->getInnerHorizontalFaces());
      Kokkos::parallel_for(
          "computeGradPhiFace1", nbInnerHorizontalFaces,
          KOKKOS_LAMBDA(const int& fInnerHorizontalFaces) {
            int fId(innerHorizontalFaces[fInnerHorizontalFaces]);
            int fFaces(utils::indexOf(mesh->getFaces(), fId));
            int cfFrontCellF(mesh->getFrontCell(fId));
            int cfId(cfFrontCellF);
            int cfCells(cfId);
            int cbBackCellF(mesh->getBackCell(fId));
            int cbId(cbBackCellF);
            int cbCells(cbId);
            // gradPhiFace1(fFaces) =
            // ArrayOperations::divide((ArrayOperations::minus(ArrayOperations::divide(ULagrange(cfCells),
            // vLagrange(cfCells)),
            // ArrayOperations::divide(ULagrange(cbCells),
            // vLagrange(cbCells)))), deltaxLagrange(fFaces));
            gradPhiFace1(fFaces) = ArrayOperations::divide(
                ArrayOperations::minus(Phi(cfCells), Phi(cbCells)),
                deltaxLagrange(fFaces));
            // Phi(cfCells) = ArrayOperations::divide(ULagrange(cfCells),
            // vLagrange(cfCells)); Phi(cbCells) =
            // ArrayOperations::divide(ULagrange(cbCells), vLagrange(cbCells));
            // std::cout << " cbCells " << cbCells << " Phi(cbCells) " <<
            // Phi(cbCells) << std::endl; std::cout << " cfCells " << cfCells <<
            // " Phi(cfCells) " << Phi(cfCells) << std::endl;
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
            // HvLagrange(cfCells) =  faceLength(frFaces);

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
    }
  }
}
/**
 * Job computeGradPhi1 called @9.0 in executeTimeLoopN method.
 * In variables: gradPhiFace1, projectionLimiterId, projectionOrder, x_then_y_n
 * Out variables: gradPhi1
 */
void EucclhydRemap::computeGradPhi1() noexcept {
  if (options->projectionOrder > 1) {
    if (x_then_y_n) {
      // std::cout << " Phase 1 Horizontale computeGradPhi1 " << std::endl;
      Kokkos::parallel_for(
          "computeGradPhi1", nbCells, KOKKOS_LAMBDA(const int& cCells) {
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

            // std::cout << " Phase 1 Horizontale " << std::endl;
            RealArray1D<dim> exy = {{1.0, 0.0}};

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
            // Phase 1 H voisinage pure =  " << voisinage_pure << std::endl;
            gradPhi1(cCells) = computeAndLimitGradPhi(
                limiter, gradPhiFace1(frFaces), gradPhiFace1(flFaces),
                Phi(cCells), Phi(cfCells), Phi(cbCells), HvLagrange(cCells),
                HvLagrange(cfCells), HvLagrange(cbCells));
            if (cCells == dbgcell3 || cCells == dbgcell2 ||
                cCells == dbgcell1) {
              std::cout << " Hori cell   " << cCells << "gradphi1 "
                        << gradPhi1(cCells) << std::endl;
              std::cout << " Grad Phi  " << gradPhiFace1(frFaces) << " "
                        << gradPhiFace1(flFaces) << std::endl;
            }

            if (options->projectionAvecPlateauPente == 1) {
              // std::cout << " cbCells " << cbCells << " Cells " << cCells << "
              // cfCells " << cfCells << std::endl; std::cout << " flFaces " <<
              // flFaces << " flId" << flId << std::endl; std::cout << " flFaces
              // " << flFaces << " faceNormalVerticalVelo " <<
              // faceNormalVelocity(flFaces) << " et" <<
              // outerFaceNormal(cCells,flFacesOfCellC) << std::endl;
              double Flux_sortant_ar =
                  MathFunctions::dot(outerFaceNormal(cCells, flFacesOfCellC),
                                     exy) *
                  faceNormalVelocity(flFaces);
              if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                   cCells == dbgcell1) &&
                  test_debug == 1) {
                std::cout << " cell   " << cCells << std::endl;
              }

              if (voisinage_pure)
                deltaPhiFaceAr(cCells) = computeIntersectionPPPure(
                    gradPhi1(cCells), Phi(cCells), Phi(cfCells), Phi(cbCells),
                    HvLagrange(cCells), HvLagrange(cfCells),
                    HvLagrange(cbCells), Flux_sortant_ar, deltat_n, 0, cCells,
                    options->threshold);
              else
                deltaPhiFaceAr(cCells) = computeIntersectionPP(
                    gradPhi1(cCells), Phi(cCells), Phi(cfCells), Phi(cbCells),
                    HvLagrange(cCells), HvLagrange(cfCells),
                    HvLagrange(cbCells), Flux_sortant_ar, deltat_n, 0, cCells,
                    options->threshold);
              if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                   cCells == dbgcell1) &&
                  test_debug == 1) {
                std::cout << " flFaces " << flFaces << " deltaphifaceAR "
                          << deltaPhiFaceAr(cCells) << std::endl;
                std::cout << " flFaces " << flFaces << " Flux_sortant_ar "
                          << Flux_sortant_ar << std::endl;
                std::cout << " frFaces " << frFaces << " frId" << frId
                          << std::endl;
              }

              // std::cout << " frFaces " << fbFaces << " faceNormalVerticalVelo
              // " <<  faceNormalVelocity(frFaces) << " et" <<
              // outerFaceNormal(cCells,frFacesOfCellC) << std::endl;
              double Flux_sortant_av =
                  MathFunctions::dot(outerFaceNormal(cCells, frFacesOfCellC),
                                     exy) *
                  faceNormalVelocity(frFaces);
              if (voisinage_pure)
                deltaPhiFaceAv(cCells) = computeIntersectionPPPure(
                    gradPhi1(cCells), Phi(cCells), Phi(cfCells), Phi(cbCells),
                    HvLagrange(cCells), HvLagrange(cfCells),
                    HvLagrange(cbCells), Flux_sortant_av, deltat_n, 1, cCells,
                    options->threshold);
              else
                deltaPhiFaceAv(cCells) = computeIntersectionPP(
                    gradPhi1(cCells), Phi(cCells), Phi(cfCells), Phi(cbCells),
                    HvLagrange(cCells), HvLagrange(cfCells),
                    HvLagrange(cbCells), Flux_sortant_av, deltat_n, 1, cCells,
                    options->threshold);
              if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                   cCells == dbgcell1) &&
                  test_debug == 1) {
                std::cout << " frFaces " << frFaces << " deltaphifaceAV "
                          << deltaPhiFaceAv(cCells) << std::endl;
                std::cout << " frFaces " << frFaces << " Flux_sortant_av "
                          << Flux_sortant_av << std::endl;
              }
            }
          });
    } else {
      // std::cout << " Phase 1 Verticale computeGradPhi1 " << std::endl;
      Kokkos::parallel_for(
          "computeGradPhi1", nbCells, KOKKOS_LAMBDA(const int& cCells) {
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
            // Phase 1 V voisinage pure =  " << voisinage_pure << std::endl;
            gradPhi1(cCells) = computeAndLimitGradPhi(
                limiter, gradPhiFace1(fbFaces), gradPhiFace1(ftFaces),
                Phi(cCells), Phi(cbCells), Phi(cfCells), HvLagrange(cCells),
                HvLagrange(cbCells), HvLagrange(cfCells));
            if (cCells == dbgcell3 || cCells == dbgcell2 ||
                cCells == dbgcell1) {
              std::cout << " Verti cell   " << cCells << "gradphi1 "
                        << gradPhi1(cCells) << std::endl;
              std::cout << " Grad Phi  " << gradPhiFace1(fbFaces) << " "
                        << gradPhiFace1(ftFaces) << std::endl;
            }
            if (options->projectionAvecPlateauPente == 1) {
              // std::cout << " Phase 1 Verticale " << std::endl;
              RealArray1D<dim> exy = {{0.0, 1.0}};
              // std::cout << " cfCells " << cfCells << " Cells " << cCells << "
              // cbCells " << cbCells << std::endl; std::cout << " fbFaces " <<
              // fbFaces << " faceNormalVerticalVelo " <<
              // faceNormalVelocity(fbFaces) << " et" <<
              // outerFaceNormal(cCells,fbFacesOfCellC) << std::endl;
              double Flux_sortant_av =
                  MathFunctions::dot(outerFaceNormal(cCells, fbFacesOfCellC),
                                     exy) *
                  faceNormalVelocity(fbFaces);
              if (voisinage_pure)
                deltaPhiFaceAr(cCells) = computeIntersectionPPPure(
                    gradPhi1(cCells), Phi(cCells), Phi(cbCells), Phi(cfCells),
                    HvLagrange(cCells), HvLagrange(cbCells),
                    HvLagrange(cfCells), Flux_sortant_av, deltat_n, 0, cCells,
                    options->threshold);
              else
                deltaPhiFaceAr(cCells) = computeIntersectionPP(
                    gradPhi1(cCells), Phi(cCells), Phi(cbCells), Phi(cfCells),
                    HvLagrange(cCells), HvLagrange(cbCells),
                    HvLagrange(cfCells), Flux_sortant_av, deltat_n, 0, cCells,
                    options->threshold);

              // std::cout << " ftFaces " << ftFaces << " ftId" << ftId <<
              // std::endl; std::cout << " ftFaces " << ftFaces << "
              // faceNormalVerticalVelo " <<  faceNormalVelocity(ftFaces) << "
              // et" << outerFaceNormal(cCells,ftFacesOfCellC) << std::endl;
              double Flux_sortant_ar =
                  MathFunctions::dot(outerFaceNormal(cCells, ftFacesOfCellC),
                                     exy) *
                  faceNormalVelocity(ftFaces);

              if (voisinage_pure)
                deltaPhiFaceAv(cCells) = computeIntersectionPPPure(
                    gradPhi1(cCells), Phi(cCells), Phi(cbCells), Phi(cfCells),
                    HvLagrange(cCells), HvLagrange(cbCells),
                    HvLagrange(cfCells), Flux_sortant_ar, deltat_n, 1, cCells,
                    options->threshold);
              else
                deltaPhiFaceAv(cCells) = computeIntersectionPP(
                    gradPhi1(cCells), Phi(cCells), Phi(cbCells), Phi(cfCells),
                    HvLagrange(cCells), HvLagrange(cbCells),
                    HvLagrange(cfCells), Flux_sortant_ar, deltat_n, 1, cCells,
                    options->threshold);
            }
          });
    }
  }
}

/**
 * Job computeUpwindFaceQuantitiesForProjection1 called @10.0 in
 * executeTimeLoopN method. In variables: ULagrange, XcLagrange, Xf,
 * deltaxLagrange, faceNormal, faceNormalVelocity, gradPhi1, vLagrange,
 * x_then_y_n Out variables: phiFace1
 */
void EucclhydRemap::computeUpwindFaceQuantitiesForProjection1() noexcept {
  if (x_then_y_n) {
    // std::cout << " Phase Projection 1 Horizontale " << std::endl;
    auto innerVerticalFaces(mesh->getInnerVerticalFaces());
    Kokkos::parallel_for(
        "computeUpwindFaceQuantitiesForProjection1", nbInnerVerticalFaces,
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
            phiFace1(fFaces) = computeUpwindFaceQuantities(
                faceNormal(fFaces), faceNormalVelocity(fFaces),
                deltaxLagrange(fFaces), Xf(fFaces),
                ArrayOperations::divide(ULagrange(cbCells), vLagrange(cbCells)),
                gradPhi1(cbCells), XcLagrange(cbCells),
                ArrayOperations::divide(ULagrange(cfCells), vLagrange(cfCells)),
                gradPhi1(cfCells), XcLagrange(cfCells));
          } else {
            phiFace1(fFaces) = ArrayOperations::minus(deltaPhiFaceAv(cbCells),
                                                      deltaPhiFaceAr(cfCells));

            // std::cout << " ----------------------------------" << std::endl;
            // std::cout << " Phase Projection 1 Horizontal " << std::endl;
            // std::cout << " cbCells " << cbCells << " cfCells " << cfCells <<
            // std::endl; std::cout << " fFaces " << fFaces << " fId" << fId <<
            // std::endl; std::cout << " phiFace1(fFaces) " << phiFace1(fFaces)
            // <<  std::endl; std::cout << " deltaPhiFaceAv(cbCells) " <<
            // deltaPhiFaceAv(cbCells) <<  std::endl; std::cout << "
            // deltaPhiFaceAr(cfCells) " << deltaPhiFaceAr(cfCells) <<
            // std::endl;
          }
          if (fFaces == face_debug1 || fFaces == face_debug2) {
            std::cout << " Phase Projection 1 Horizontale " << std::endl;
            std::cout << " cbCells " << cbCells << " cfCells " << cfCells
                      << std::endl;
            std::cout << " fFaces " << fFaces << " fId" << fId << std::endl;
            std::cout << " phiFace1(fFaces) " << phiFace1(fFaces) << std::endl;
            std::cout << " deltaPhiFaceAv(cbCells) " << deltaPhiFaceAv(cbCells)
                      << std::endl;
            std::cout << " deltaPhiFaceAr(cfCells) " << deltaPhiFaceAr(cfCells)
                      << std::endl;
            std::cout << " cbcell   " << cbCells << " Ulagrange "
                      << ULagrange(cbCells) << std::endl;
            std::cout << " cfcell   " << cfCells << " Ulagrange "
                      << ULagrange(cfCells) << std::endl;
          }
        });
  } else {
    // std::cout << " Phase Projection 1 Verticale " << std::endl;
    auto innerHorizontalFaces(mesh->getInnerHorizontalFaces());
    Kokkos::parallel_for(
        "computeUpwindFaceQuantitiesForProjection1", nbInnerHorizontalFaces,
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
            phiFace1(fFaces) = computeUpwindFaceQuantities(
                faceNormal(fFaces), faceNormalVelocity(fFaces),
                deltaxLagrange(fFaces), Xf(fFaces),
                ArrayOperations::divide(ULagrange(cbCells), vLagrange(cbCells)),
                gradPhi1(cbCells), XcLagrange(cbCells),
                ArrayOperations::divide(ULagrange(cfCells), vLagrange(cfCells)),
                gradPhi1(cfCells), XcLagrange(cfCells));
          } else {
            phiFace1(fFaces) = ArrayOperations::minus(deltaPhiFaceAv(cfCells),
                                                      deltaPhiFaceAr(cbCells));
            // cfCells maille arriere, cbCells maille devant
            // std::cout << " ----------------------------------" << std::endl;
            // std::cout << " Phase Projection 1 Verticale " << std::endl;
            // std::cout << " cfCells " << cfCells << " cbCells " << cbCells <<
            // std::endl; std::cout << " fFaces " << fFaces << " fId" << fId <<
            // std::endl; std::cout << " phiFace1(fFaces) " << phiFace1(fFaces)
            // <<  std::endl; std::cout << " deltaPhiFaceAv(cfCells) " <<
            // deltaPhiFaceAv(cfCells) <<  std::endl; std::cout << "
            // deltaPhiFaceAr(cbCells) " << deltaPhiFaceAr(cbCells) <<
            // std::endl;
          }
          if (fFaces == face_debug1 || fFaces == face_debug2) {
            std::cout << " Phase Projection 1 Verticale " << std::endl;
            std::cout << " cfCells " << cfCells << " cbCells " << cbCells
                      << std::endl;
            std::cout << " fFaces " << fFaces << " fId" << fId << std::endl;
            std::cout << " phiFace1(fFaces) " << phiFace1(fFaces) << std::endl;
            std::cout << " deltaPhiFaceAv(cfCells) " << deltaPhiFaceAv(cfCells)
                      << std::endl;
            std::cout << " deltaPhiFaceAr(cbCells) " << deltaPhiFaceAr(cbCells)
                      << std::endl;
          }
        });
  }
}
/**
 * Job computeUremap1 called @11.0 in executeTimeLoopN method.
 * In variables: ULagrange, deltat_n, faceLength, faceNormal,
 * faceNormalVelocity, outerFaceNormal, phiFace1, x_then_y_n Out variables:
 * Uremap1
 */
void EucclhydRemap::computeUremap1() noexcept {
  RealArray1D<dim> exy = xThenYToDirection(x_then_y_n);
  Kokkos::parallel_for(
      "computeUremap1", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        int cId(cCells);

        // std::cout << " cCells " << cCells << std::endl;
        // std::cout << "ULagrange " << ULagrange(cCells) << std::endl;
        RealArray1D<nbequamax> reduction8 = options->Uzero;
        {
          auto neighbourCellsC(mesh->getNeighbourCells(cId));
          for (int dNeighbourCellsC = 0;
               dNeighbourCellsC < neighbourCellsC.size(); dNeighbourCellsC++) {
            int dId(neighbourCellsC[dNeighbourCellsC]);
            int fCommonFaceCD(mesh->getCommonFace(cId, dId));
            int fId(fCommonFaceCD);
            int fFaces(utils::indexOf(mesh->getFaces(), fId));
            int fFacesOfCellC(utils::indexOf(mesh->getFacesOfCell(cId), fId));
            reduction8 = ArrayOperations::plus(
                reduction8,
                (computeRemapFlux(
                    options->projectionAvecPlateauPente,
                    faceNormalVelocity(fFaces), faceNormal(fFaces),
                    faceLength(fFaces), phiFace1(fFaces),
                    outerFaceNormal(cCells, fFacesOfCellC), exy, deltat_n)));
            //
            if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                 cCells == dbgcell1)) {
              std::cout << " --------- flux face interieur "
                           "----------------------------------"
                        << std::endl;
              std::cout << " cell   " << cCells << "reduction8 " << reduction8
                        << std::endl;
              std::cout << " face " << fFaces << " phiFace1(fFaces) "
                        << phiFace1(fFaces) << " de longueur "
                        << faceLength(fFaces) << " "
                        << computeRemapFlux(
                               options->projectionAvecPlateauPente,
                               faceNormalVelocity(fFaces), faceNormal(fFaces),
                               faceLength(fFaces), phiFace1(fFaces),
                               outerFaceNormal(cCells, fFacesOfCellC), exy,
                               deltat_n)
                        << " outer " << outerFaceNormal(cCells, fFacesOfCellC)
                        << "facenormal" << faceNormalVelocity(fFaces) << " dot"
                        << MathFunctions::dot(
                               outerFaceNormal(cCells, fFacesOfCellC), exy)
                        << std::endl;
              std::cout << " --------- ----------------------------------"
                        << std::endl;
            }
            //
          }
          if (options->FluxBC > 0) {
            // flux exterieur eventuel
            if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                 cCells == dbgcell1)) {
              std::cout << " --------- flux face exterieur "
                           "----------------------------------"
                        << std::endl;
              std::cout << exy << " AV cell   " << cCells << "reduction8 "
                        << reduction8 << std::endl;
            }
            //
            reduction8 = ArrayOperations::plus(
                reduction8, (computeBoundaryFluxes(1, cCells, exy)));
            //
            if ((cCells == dbgcell3 || cCells == dbgcell2 ||
                 cCells == dbgcell1)) {
              std::cout << exy << " AP cell   " << cCells << "reduction8 "
                        << reduction8 << std::endl;
              std::cout << computeBoundaryFluxes(1, cCells, exy) << std::endl;
              std::cout << " --------- ----------------------------------"
                        << std::endl;
            }
          }
        }
        Uremap1(cCells) = ArrayOperations::minus(ULagrange(cCells), reduction8);

        if ((cCells == dbgcell3 || cCells == dbgcell2 || cCells == dbgcell1)) {
          std::cout << " cell   " << cCells << "ULagrange " << ULagrange(cCells)
                    << std::endl;
          std::cout << " cell   " << cCells << "reduction8 " << reduction8
                    << std::endl;
          std::cout << " cell   " << cCells << "Uremap1 " << Uremap1(cCells)
                    << std::endl;
          std::cout << " cell   " << cCells << "Phi" << Phi(cCells)
                    << std::endl;
        }

        if (options->projectionAvecPlateauPente == 1) {
          // option ou on ne regarde pas la variation de rho, V et e
          // phi = (f1, f2, rho1*f1, rho2*f2, Vx, Vy, e1, e2
          // ce qui permet d'ecrire le flux telque
          // Flux = (dv1 = f1dv, dv2=f2*dv, dm1=rho1*df1, dm2=rho2*df2, d(mVx) =
          // Vx*(dm1+dm2), d(mVy) = Vy*(dm1+dm2), d(m1e1) = e1*dm1,  d(m2e2) =
          // e2*dm2 dans computeIntersectionPP

          double somme_volume = 0.;
          for (imat = 0; imat < nbmatmax; imat++) {
            somme_volume += Uremap1(cCells)[imat];
          }
          // Phi volume
          double somme_masse = 0.;
          for (imat = 0; imat < nbmatmax; imat++) {
            Phi(cCells)[imat] = Uremap1(cCells)[imat] / somme_volume;
            // Phi masse
            if (Uremap1(cCells)[imat] != 0.)
              Phi(cCells)[nbmatmax + imat] =
                  Uremap1(cCells)[nbmatmax + imat] / (Uremap1(cCells)[imat]);
            else
              Phi(cCells)[nbmatmax + imat] = 0.;
            somme_masse += Uremap1(cCells)[nbmatmax + imat];
          }
          // Phi Vitesse
          Phi(cCells)[2 * nbmatmax] =
              Uremap1(cCells)[2 * nbmatmax] / somme_masse;
          Phi(cCells)[2 * nbmatmax + 1] =
              Uremap1(cCells)[2 * nbmatmax + 1] / somme_masse;
          // Phi energie
          for (imat = 0; imat < nbmatmax; imat++) {
            if (Uremap1(cCells)[nbmatmax + imat] != 0.)
              Phi(cCells)[2 * nbmatmax + imat + 2] =
                  Uremap1(cCells)[2 * nbmatmax + imat + 2] /
                  Uremap1(cCells)[nbmatmax + imat];
            else
              Phi(cCells)[2 * nbmatmax + imat + 2] = 0.;
          }
          // Phi energie cinétique
          if (options->projectionConservative == 1)
            Phi(cCells)[3 * nbmatmax + 2] =
                Uremap1(cCells)[3 * nbmatmax + 2] / somme_masse;

        } else {
          Phi(cCells) =
              ArrayOperations::divide(Uremap1(cCells), vLagrange(cCells));
        }

        if ((cCells == dbgcell3 || cCells == dbgcell2 || cCells == dbgcell1)) {
          std::cout << " cell   " << cCells << "Phi apres Uremap1"
                    << Phi(cCells) << std::endl;
        }
        // Mises à jour de l'indicateur mailles mixtes
        int matcell(0);
        int imatpure(-1);
        for (int imat = 0; imat < nbmatmax; imat++)
          if (Phi(cCells)[imat] > 0.) {
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
