#include "arcane/accelerator/Accelerator.h"

#include "accenv/ProfAcc.h"

/*---------------------------------------------------------------------------*/
/* Pour le profiling sur accélérateur                                        */
/*---------------------------------------------------------------------------*/

#if defined(ARCANE_COMPILING_CUDA) && defined(PROF_ACC)

#ifdef P4GPU_HAS_WARNING_INFO
#warning "API PROF_ACC : instrumentation avec nvtx"
#endif
#include <nvtx3/nvToolsExt.h>
#include <cuda_profiler_api.h>

void prof_acc_start_capture() {
  cudaProfilerStart();
}

void prof_acc_stop_capture() {
  cudaProfilerStop();
}

void prof_acc_begin(const char* name) {
  nvtxRangePushA(name);
}

void prof_acc_end([[maybe_unused]] const char* name) {
  nvtxRangePop();
}

#elif defined(ARCANE_COMPILING_HIP) && defined(PROF_ACC)

#ifdef P4GPU_HAS_WARNING_INFO
#warning "API PROF_ACC : instrumentation avec roctx"
#endif
#include <roctracer/roctx.h>
//#include <roctracer/roctracer_ext.h>

void prof_acc_start_capture() {
  //roctracer_start();
}

void prof_acc_stop_capture() {
  //roctracer_stop();
}

void prof_acc_begin(const char* name) {
  roctxRangePushA(name);
}

void prof_acc_end([[maybe_unused]] const char* name) {
  roctxRangePop();
}

#else

#ifdef P4GPU_HAS_WARNING_INFO
#warning "API PROF_ACC : Pas d'instrumentation"
#endif
void prof_acc_start_capture() {
  return;
}

void prof_acc_stop_capture() {
  return;
}

void prof_acc_begin([[maybe_unused]] const char* name) {
  return;
}

void prof_acc_end([[maybe_unused]] const char* name) {
  return;
}

#endif

