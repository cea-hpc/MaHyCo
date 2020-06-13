#include <math.h>    // for pow, floor, sqrt
#include <stdlib.h>  // for exit
#include <Kokkos_Core.hpp>
#include <algorithm>                // for max, min
#include <array>                    // for array
#include <iostream>                 // for operator<<, basic_ostream, basic_...
#include <vector>                   // for allocator, vector
#include "EucclhydRemap.h"          // for EucclhydRemap, EucclhydRemap::Opt...
#include "types/ArrayOperations.h"  // for multiply, minus, plus
#include "types/MathFunctions.h"    // for max

void EucclhydRemap::updateParticlePosition() noexcept {
  Kokkos::parallel_for(
      "updateParticleCoefficient", nbCells,
      KOKKOS_LAMBDA(const int& cCells) { listpart(cCells).clear(); });
  Kokkos::parallel_for("initPart", nbPart, KOKKOS_LAMBDA(const int& ipart) {
    Xpart_nplus1(ipart) = ArrayOperations::plus(
        Xpart_n(ipart), ArrayOperations::multiply(Vpart_n(ipart), deltat_n));
    if (Xpart_nplus1(ipart)[1] < 0.) {
      Xpart_nplus1(ipart)[1] = -Xpart_nplus1(ipart)[1];
      Vpart_n(ipart)[1] = -Vpart_n(ipart)[1];
    }
    int icell = MathFunctions::max(
        floor(Xpart_nplus1(ipart)[0] / options->X_EDGE_LENGTH), 0);
    int jcell = MathFunctions::max(
        floor(Xpart_nplus1(ipart)[1] / options->Y_EDGE_LENGTH), 0);
    ICellp(ipart) = jcell * options->X_EDGE_ELEMS + icell;

    // conditions limites
    if (fracvol(ICellp(ipart))[IMatp(ipart)] < 0.25) {
      RealArray1D<dim> gradf = options->zeroVect;
      if (IMatp(ipart) == 0) gradf = gradf1(ICellp(ipart));
      if (IMatp(ipart) == 1) gradf = gradf2(ICellp(ipart));
      if (IMatp(ipart) == 2) gradf = gradf3(ICellp(ipart));

      // std::cout << " CL : Part  " << ipart << " vit " << Vpart_n(ipart)
      //	    << "  xp= " << Xpart_nplus1(ipart)[0] << "  yp= " <<
      // Xpart_nplus1(ipart)[1] << "  " << ICellp(ipart) << std::endl;

      // particule sort du materiau IMatp(ipart)
      // changement de la vitesse
      double module_vit = sqrt(Vpart_n(ipart)[0] * Vpart_n(ipart)[0] +
                               Vpart_n(ipart)[1] * Vpart_n(ipart)[1]);
      double module_gradf = sqrt(gradf[0] * gradf[0] + gradf[1] * gradf[1]);
      Vpart_n(ipart) =
          ArrayOperations::multiply(gradf, 1.1 * module_vit / module_gradf);
      Xpart_nplus1(ipart) = ArrayOperations::plus(
          Xpart_n(ipart), ArrayOperations::multiply(Vpart_n(ipart), deltat_n));
      if (Xpart_nplus1(ipart)[1] < 0.) {
        Xpart_nplus1(ipart)[1] = -Xpart_nplus1(ipart)[1];
        Vpart_n(ipart)[1] = -Vpart_n(ipart)[1];
      }
      int icell = MathFunctions::max(
          floor(Xpart_nplus1(ipart)[0] / options->X_EDGE_LENGTH), 0);
      int jcell = MathFunctions::max(
          floor(Xpart_nplus1(ipart)[1] / options->Y_EDGE_LENGTH), 0);
      ICellp(ipart) = jcell * options->X_EDGE_ELEMS + icell;

      // std::cout << " AP : Part  " << ipart << " vit " << Vpart_n(ipart)
      //	    << "  xp= " << Xpart_nplus1(ipart)[0] << "  yp= " <<
      // Xpart_nplus1(ipart)[1] << "  " << ICellp(ipart) << std::endl;
      // std::cout << " AP : Part  " << ipart << " gradf " << gradf << "
      // module_vit " << module_vit << std::endl;

      if (fracvol(ICellp(ipart))[IMatp(ipart)] < 0.05)
        std::cout << " Part  " << ipart << " non recuperee " << std::endl;
    }

    listpart(ICellp(ipart)).push_back(ipart);
    // std::cout << " Part  " << ipart << " vit " << Vpart_nplus1(ipart) <<
    // "  xp= " << Xpart_nplus1(ipart)[0] << "  yp= " <<
    // Xpart_nplus1(ipart)[1] << "  " << ICellp(ipart)
    //	  << " " << listpart(ICellp(ipart)).size() << std::endl;
  });
}

void EucclhydRemap::updateParticleCoefficients() noexcept {
  Kokkos::parallel_for(
      "updateParticleCoefficient", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        fracpart(cCells) = 1.;
        for (int ipart = 0; ipart < listpart(cCells).size(); ipart++) {
          fracpart(cCells) -= wpart(ipart) * vpart(ipart) / v(cCells);
          fracpart(cCells) = max(fracpart(cCells), 0.);
          if (fracpart(cCells) == 0.) {
            std::cout << cCells << "  " << listpart(cCells).size() << " v "
                      << v(cCells) << " " << fracpart(cCells) << std::endl;
            std::cout << " Plus de gaz dans la maille -> fin du calcul "
                      << std::endl;
            exit(1);
          }
        }
        // if (listpart(cCells).size() > 0) std::cout << cCells << "  " <<
        // listpart(cCells).size() << " v " << v(cCells) << " " <<
        // fracpart(cCells) << std::endl;
      });

  Kokkos::parallel_for(
      "updateParticleCoefficient", nbPart, KOKKOS_LAMBDA(const int& ipart) {
        int icells = ICellp(ipart);
        RealArray1D<dim> VLagrange = V_n(icells);
        double normvpvg = (Vpart_n(ipart)[0] - VLagrange[0]) *
                              (Vpart_n(ipart)[0] - VLagrange[0]) +
                          (Vpart_n(ipart)[1] - VLagrange[1]) *
                              (Vpart_n(ipart)[1] - VLagrange[1]);
        Repart(ipart) = max(2. * rho_n(icells) * rpart(ipart) *
                                MathFunctions::sqrt(normvpvg) / viscosity,
                            options->Reynolds_min);
        if (options->DragModel == options->Kliatchko) {
          if (Repart(ipart) > 1000)
            Cdpart(ipart) = 0.02 * pow(fracpart(icells), (-2.55)) +
                            0.4 * pow(fracpart(icells), (-1.78));
          else
            Cdpart(ipart) =
                24. / (Repart(ipart) * pow(fracpart(icells), (2.65))) +
                4. / (pow(Repart(ipart), (1. / 3.)) *
                      pow(fracpart(icells), (1.78)));
        } else {
          double fMac_part;
          if (Mcpart(ipart) < 0.5)
            fMac_part = 1.;
          else if (Mcpart(ipart) >= 0.5 && Mcpart(ipart) < 1.0)
            fMac_part = 2.44 * Mcpart(ipart) - 0.2222;
          else if (Mcpart(ipart) >= 1.0)
            fMac_part = 0.2222;
          Cdpart(ipart) =
              fMac_part *
              pow((24. + 4. * min(Repart(ipart), options->Reynolds_max)),
                  2.13) /
              (min(Repart(ipart), options->Reynolds_min));
        }
        // std::cout << " Part  " << ipart << " cd " << Cdpart(ipart)  << "  Re
        // " << Repart(ipart) << " " << fracpart(icells) << std::endl;
      });
}

void EucclhydRemap::updateParticleVelocity() noexcept {
  Kokkos::parallel_for("initPart", nbPart, KOKKOS_LAMBDA(const int& ipart) {
    int icells = ICellp(ipart);
    Vpart_nplus1(ipart) = ArrayOperations::minus(
        Vpart_n(ipart), ArrayOperations::multiply(deltat_n / rhopart(ipart),
                                                  ForceGradp(icells)));

    // std::cout << " Part  " << ipart << " vit " << Vpart_nplus1(ipart)  <<
    // std::endl;

    RealArray1D<dim> VLagrange = V_n(icells);
    double normvpvg =
        (Vpart_n(ipart)[0] - VLagrange[0]) *
            (Vpart_n(ipart)[0] - VLagrange[0]) +
        (Vpart_n(ipart)[1] - VLagrange[1]) * (Vpart_n(ipart)[1] - VLagrange[1]);
    double drag = (3 * rho_n(icells) * Cdpart(ipart)) /
                  (8 * rhopart(ipart) * rpart(ipart)) *
                  MathFunctions::sqrt(normvpvg);
    // drag = 2000.;

    // Vpart_nplus1(ipart) = ArrayOperations::minus(Vpart_nplus1(ipart),
    // ArrayOperations::multiply(options->Drag,
    // (ArrayOperations::minus(Vpart_n(ipart) - VLagrange))));
    Vpart_nplus1(ipart)[0] =
        Vpart_nplus1(ipart)[0] -
        drag * (Vpart_n(ipart)[0] - VLagrange[0]) * deltat_n;
    Vpart_nplus1(ipart)[1] =
        Vpart_nplus1(ipart)[1] -
        drag * (Vpart_n(ipart)[1] - VLagrange[1]) * deltat_n;
    // if (abs(Vpart_nplus1(ipart)[0]) > options->threshold) std::cout << "
    // Part  " << ipart << " vit " << Vpart_nplus1(ipart)  << std::endl;
  });
}

void EucclhydRemap::updateParticleRetroaction() noexcept {
  Kokkos::parallel_for("initPart", nbPart, KOKKOS_LAMBDA(const int& ipart) {
    int icells = ICellp(ipart);
    RealArray1D<dim> AccelerationP =
        ArrayOperations::minus(Vpart_n(ipart), Vpart_nplus1(ipart));
    V_nplus1(icells)[0] += mpart(ipart) * AccelerationP[0] / m(icells);
    V_nplus1(icells)[1] += mpart(ipart) * AccelerationP[1] / m(icells);
  });
}

void EucclhydRemap::switchalpharho_rho() noexcept {
  Kokkos::parallel_for("updateParticleCoefficient", nbCells,
                       KOKKOS_LAMBDA(const int& cCells) {
                         rho_n(cCells) /= fracpart(cCells);
                         for (imat = 0; imat < nbmatmax; imat++)
                           rhop_n(cCells)[imat] /= fracpart(cCells);
                       });
}

void EucclhydRemap::switchrho_alpharho() noexcept {
  Kokkos::parallel_for("updateParticleCoefficient", nbCells,
                       KOKKOS_LAMBDA(const int& cCells) {
                         rho_n(cCells) *= fracpart(cCells);
                         for (imat = 0; imat < nbmatmax; imat++)
                           rhop_n(cCells)[imat] *= fracpart(cCells);
                       });
}
