#ifdef MSG_PASS_HAS_MPI
#include <arcane/utils/ArcaneGlobal.h>

#include <iostream>
#include <mpi.h>

using namespace Arcane;

void msg_pass_init([[maybe_unused]] int* argc, [[maybe_unused]] char*** argv) {
//#define USE_THREAD_MULTIPLE
#ifdef USE_THREAD_MULTIPLE
  // On impose le niveau MPI_THREAD_MULTIPLE
  Integer provided;
  MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
  ARCANE_ASSERT(provided==MPI_THREAD_MULTIPLE, 
      ("Impossible d'utiliser le niveau MPI_THREAD_MULTIPLE"));
  Integer rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank==0) {
    if (provided==MPI_THREAD_SERIALIZED)
      std::cout << "MPI_THREAD_SERIALIZED" << std::endl;
    else if (provided==MPI_THREAD_MULTIPLE)
      std::cout << "MPI_THREAD_MULTIPLE" << std::endl;
    else
      std::cout << "MPI_THREAD_{SINGLE|FUNNELED}" << std::endl;
  }
  // Fin init MPI en dur
#endif
}
#else
#ifdef P4GPU_HAS_WARNING_INFO
#warning "MPI non détecté : msg_pass_init(...) vide"
#endif
void msg_pass_init(int* , char*** ) {
}
#endif
