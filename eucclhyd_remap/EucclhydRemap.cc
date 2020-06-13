#include "EucclhydRemap.h"
#include <stdlib.h>               // for exit
#include <Kokkos_Core.hpp>        // for finalize
#include <iomanip>                // for operator<<, setw, setiosflags
#include <iostream>               // for operator<<, basic_ostream, cha...
#include <limits>                 // for numeric_limits
#include <map>                    // for map
#include <utility>                // for pair, swap
#include "types/MathFunctions.h"  // for min
#include "utils/Utils.h"          // for __RESET__, __BOLD__, __GREEN__

using namespace nablalib;

/**
 * Job dumpVariables called @2.0 in executeTimeLoopN method.
 * In variables: Xc_x, Xc_y, eps_n, m, p, rho_n, t_n, v
 * Out variables:
 */
KOKKOS_INLINE_FUNCTION
void EucclhydRemap::dumpVariables() noexcept {
  std::cout << " Deltat = " << deltat_n << std::endl;
  std::cout << " ---------------------------"
            << " Energie totale(t=0) = " << ETOTALE_0
            << " Energie totale(Lagrange) = " << ETOTALE_L
            << " Energie totale(time) = " << ETOTALE_T << std::endl;
  std::cout << " ---------------------------"
            << " Masse totale(t=0) = " << MASSET_0
            << " Masse totale(Lagrange) = " << MASSET_L
            << " Masse totale(time) = " << MASSET_T << std::endl;
  nbCalls++;
  if (!writer.isDisabled() &&
      (t_n >= lastDump + options->output_time || t_n == 0.)) {
    cpu_timer.stop();
    io_timer.start();
    std::map<string, double*> cellVariables;
    std::map<string, double*> nodeVariables;
    std::map<string, double*> partVariables;
    cellVariables.insert(pair<string, double*>("Pressure", p.data()));
    cellVariables.insert(pair<string, double*>("Density", rho_n.data()));
    cellVariables.insert(pair<string, double*>("F1", fracvol1.data()));
    cellVariables.insert(pair<string, double*>("F2", fracvol2.data()));
    cellVariables.insert(pair<string, double*>("F3", fracvol3.data()));
    cellVariables.insert(pair<string, double*>("VelocityX", Vxc.data()));
    cellVariables.insert(pair<string, double*>("VelocityY", Vyc.data()));
    cellVariables.insert(pair<string, double*>("Energy", eps_n.data()));
    partVariables.insert(pair<string, double*>("VolumePart", vpart.data()));
    partVariables.insert(pair<string, double*>("VxPart", Vpart_n[0].data()));
    partVariables.insert(pair<string, double*>("VyPart", Vpart_n[1].data()));
    auto quads = mesh->getGeometry()->getQuads();
    writer.writeFile(nbCalls, t_n, nbNodes, X.data(), nbCells, quads.data(),
                     cellVariables, nodeVariables);
    writerpart.writeFile(nbCalls, t_n, nbPart, Xpart_n.data(), 0, quads.data(),
                         cellVariables, partVariables);
    lastDump = t_n;
    std::cout << " time = " << t_n << " sortie demandée " << std::endl;
    io_timer.stop();
    cpu_timer.start();
  }
}

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
KOKKOS_INLINE_FUNCTION
void EucclhydRemap::executeTimeLoopN() noexcept {
  n = 0;
  bool continueLoop = true;
  do {
    global_timer.start();
    cpu_timer.start();
    n++;
    if (n != 1)
      std::cout << "[" << __CYAN__ << __BOLD__ << setw(3) << n
                << __RESET__ "] time = " << __BOLD__
                << setiosflags(std::ios::scientific) << setprecision(8)
                << setw(16) << t_n << __RESET__;

    if (nbPart != 0) switchalpharho_rho();
    computeEOS();  // @1.0
    if (nbPart != 0) switchrho_alpharho();
    computeGradients();                           // @1.0
    computeMass();                                // @1.0
    computeDissipationMatrix();                   // @2.0
    computedeltatc();                             // @2.0
    dumpVariables();                              // @2.0
    extrapolateValue();                           // @2.0
    computeG();                                   // @3.0
    computeNodeDissipationMatrixAndG();           // @3.0
    computedeltat();                              // @3.0
    computeBoundaryNodeVelocities();              // @4.0
    computeNodeVelocity();                        // @4.0
    updateTime();                                 // @4.0
    computeFaceVelocity();                        // @5.0
    computeLagrangePosition();                    // @5.0
    computeSubCellForce();                        // @5.0
    computeLagrangeVolumeAndCenterOfGravity();    // @6.0
    computeFacedeltaxLagrange();                  // @7.0
    updateCellCenteredLagrangeVariables();        // @7.0
    computeGradPhiFace1();                        // @8.0
    computeGradPhi1();                            // @9.0
    computeUpwindFaceQuantitiesForProjection1();  // @10.0
    computeUremap1();                             // @11.0
    computeGradPhiFace2();                        // @12.0
    computeGradPhi2();                            // @13.0
    computeUpwindFaceQuantitiesForProjection2();  // @14.0
    computeUremap2();                             // @15.0
    remapCellcenteredVariable();                  // @16.0
    if (nbPart != 0) {
      updateParticlePosition();
      updateParticleCoefficients();
      updateParticleVelocity();
      updateParticleRetroaction();
    }

    // Evaluate loop condition with variables at time n
    continueLoop = (n + 1 < options->max_time_iterations &&
                    t_nplus1 < options->final_time);

    if (continueLoop) {
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
      std::cout << " {CPU: " << __BLUE__ << cpu_timer.print(true)
                << __RESET__ ", IO: " << __BLUE__ << io_timer.print(true)
                << __RESET__ "} ";
    else
      std::cout << " {CPU: " << __BLUE__ << cpu_timer.print(true)
                << __RESET__ ", IO: " << __RED__ << "none" << __RESET__ << "} ";

    // Progress
    std::cout << utils::progress_bar(n, options->max_time_iterations, t_n,
                                     options->final_time, 30);
    std::cout << __BOLD__ << __CYAN__
              << utils::Timer::print(
                     utils::eta(n, options->max_time_iterations, t_n,
                                options->final_time, deltat_n, global_timer),
                     true)
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
void EucclhydRemap::computedeltat() noexcept {
  double reduction10(numeric_limits<double>::max());
  {
    Kokkos::Min<double> reducer(reduction10);
    Kokkos::parallel_reduce("reduction10", nbCells,
                            KOKKOS_LAMBDA(const int& cCells, double& x) {
                              reducer.join(x, deltatc(cCells));
                            },
                            reducer);
  }
  deltat_nplus1 =
      MathFunctions::min(options->cfl * reduction10, deltat_n * 1.05);
  if (deltat_nplus1 < options->deltat_min) {
    std::cerr << "Fin de la simulation par pas de temps minimum "
              << deltat_nplus1 << " < " << options->deltat_min << std::endl;
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
void EucclhydRemap::updateTime() noexcept { t_nplus1 = t_n + deltat_nplus1; }

void EucclhydRemap::simulate() {
  std::cout << "\n"
            << __BLUE_BKG__ << __YELLOW__ << __BOLD__
            << "\tStarting EucclhydRemap ..." << __RESET__ << "\n\n";

  std::cout << "[" << __GREEN__ << "MESH" << __RESET__
            << "]      X=" << __BOLD__ << options->X_EDGE_ELEMS << __RESET__
            << ", Y=" << __BOLD__ << options->Y_EDGE_ELEMS << __RESET__
            << ", X length=" << __BOLD__ << options->X_EDGE_LENGTH << __RESET__
            << ", Y length=" << __BOLD__ << options->Y_EDGE_LENGTH << __RESET__
            << std::endl;

  if (Kokkos::hwloc::available()) {
    std::cout << "[" << __GREEN__ << "TOPOLOGY" << __RESET__
              << "]  NUMA=" << __BOLD__
              << Kokkos::hwloc::get_available_numa_count() << __RESET__
              << ", Cores/NUMA=" << __BOLD__
              << Kokkos::hwloc::get_available_cores_per_numa() << __RESET__
              << ", Threads/Core=" << __BOLD__
              << Kokkos::hwloc::get_available_threads_per_core() << __RESET__
              << std::endl;
  } else {
    std::cout << "[" << __GREEN__ << "TOPOLOGY" << __RESET__
              << "]  HWLOC unavailable cannot get topological informations"
              << std::endl;
  }

  // std::cout << "[" << __GREEN__ << "KOKKOS" << __RESET__ << "]    " <<
  // __BOLD__ << (is_same<MyLayout,Kokkos::LayoutLeft>::value?"Left":"Right")"
  // << __RESET__ << " layout" << std::endl;

  if (!writer.isDisabled())
    std::cout << "[" << __GREEN__ << "OUTPUT" << __RESET__
              << "]    VTK files stored in " << __BOLD__
              << writer.outputDirectory() << __RESET__ << " directory"
              << std::endl;
  else
    std::cout << "[" << __GREEN__ << "OUTPUT" << __RESET__ << "]    "
              << __BOLD__ << "Disabled" << __RESET__ << std::endl;

  computeCornerNormal();       // @1.0
  initMeshGeometryForCells();  // @1.0
  initVpAndFpc();              // @1.0
  initBoundaryConditions();
  initCellInternalEnergy();    // @2.0
  initCellVelocity();          // @2.0
  initDensity();               // @2.0
  initMeshGeometryForFaces();  // @2.0
  if (nbPart != 0) {
    initPart();
    updateParticleCoefficients();
    switchrho_alpharho();  // on travaille avec alpharho sauf pour l'EOS
  }
  setUpTimeLoopN();    // @3.0
  executeTimeLoopN();  // @4.0
  std::cout << __YELLOW__ << "\n\tDone ! Took " << __MAGENTA__ << __BOLD__
            << global_timer.print() << __RESET__ << std::endl;
}