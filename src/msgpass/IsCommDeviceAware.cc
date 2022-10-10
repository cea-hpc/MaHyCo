#if defined(MSG_PASS_HAS_MPI) && defined(MSG_PASS_HAS_MPI_EXT)
#include <mpi.h>
#include <mpi-ext.h>

/*! \brief Retourne vrai si les comms peuvent se faire en utilisant directement les adresses sur le DEVICE
 */
bool is_comm_device_aware() {
#if defined(MPIX_CUDA_AWARE_SUPPORT) 
#ifdef P4GPU_HAS_WARNING_INFO
#warning "is_comm_device_aware() utilise MPIX_CUDA_AWARE_SUPPORT()"
#endif
  return (1 == MPIX_Query_cuda_support());
#else 
#ifdef P4GPU_HAS_WARNING_INFO
#warning "is_comm_device_aware() renvoie toujours false"
#endif
  return false;
#endif 
}

#elif defined(MSG_PASS_HAS_COMM_DEVICE_AWARE)
#ifdef P4GPU_HAS_WARNING_INFO
#warning "MSG_PASS_HAS_COMM_DEVICE_AWARE : is_comm_device_aware() return true"
#endif
bool is_comm_device_aware() {
  return true;
}
#else
#ifdef P4GPU_HAS_WARNING_INFO
#warning "MPI ou mpi-ext.h non détecté : is_comm_device_aware() return false"
#endif
bool is_comm_device_aware() {
  return false;
}
#endif

