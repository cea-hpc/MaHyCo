#if defined(MSG_PASS_HAS_MPI) && defined(MSG_PASS_HAS_MPI_EXT)
#include <mpi.h>
#include <mpi-ext.h>

/*! \brief Retourne vrai si les comms peuvent se faire en utilisant directement les adresses sur le DEVICE
 */
bool is_comm_device_aware() {
#if defined(MPIX_CUDA_AWARE_SUPPORT) 
#warning "is_comm_device_aware() utilise MPIX_CUDA_AWARE_SUPPORT()"
  return (1 == MPIX_Query_cuda_support());
#else 
#warning "is_comm_device_aware() renvoie toujours false"
  return false;
#endif 
}

#else
#warning "MPI ou mpi-ext.h non détecté : is_comm_device_aware() return false"
bool is_comm_device_aware() {
  return false;
}
#endif

