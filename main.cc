#include <Kokkos_Core.hpp>                  // for finalize, initialize
#include <cstdlib>                          // for atof, atoi, exit
#include <iostream>                         // for operator<<, endl, basic_o...
#include <string>                           // for string
#include "EucclhydRemap.h"                  // for EucclhydRemap::Options
#include "mesh/CartesianMesh2D.h"           // for CartesianMesh2D
#include "mesh/CartesianMesh2DGenerator.h"  // for CartesianMesh2DGenerator

using namespace nablalib;

int main(int argc, char* argv[]) {
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
    std::cerr << "[ERROR] Wrong number of arguments. Expecting 4 or 5 args: X "
                 "Y Xlength Ylength (output)."
              << std::endl;
    std::cerr << "(X=100, Y=10, Xlength=0.01, Ylength=0.01 output=current "
                 "directory with no args)"
              << std::endl;
    exit(1);
  }

  auto nm = CartesianMesh2DGenerator::generate(
      o->X_EDGE_ELEMS, o->Y_EDGE_ELEMS, o->X_EDGE_LENGTH, o->Y_EDGE_LENGTH);
  auto c = new EucclhydRemap(o, nm, output);
  c->simulate();
  delete c;
  delete nm;
  delete o;
  Kokkos::finalize();
  return 0;
}
