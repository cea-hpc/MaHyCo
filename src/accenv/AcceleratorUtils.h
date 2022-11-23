#ifndef ACCELERATOR_UTILS_H
#define ACCELERATOR_UTILS_H

#include "arcane/IApplication.h"
#include "arcane/accelerator/Reduce.h"
#include "arcane/accelerator/Runner.h"
#include "arcane/accelerator/VariableViews.h"
#include "arcane/accelerator/Accelerator.h"
#include "arcane/accelerator/RunCommandLoop.h"
#include "arcane/accelerator/RunCommandEnumerate.h"
#include "arcane/accelerator/core/RunQueueBuildInfo.h"
#include <arcane/accelerator/core/Memory.h>

#include "accenv/ProfAcc.h"

#include "arcane_version.h"

/*---------------------------------------------------------------------------*/
/* Pour les accélérateurs                                                    */
/*---------------------------------------------------------------------------*/

using namespace Arcane;
namespace ax = Arcane::Accelerator;

/*---------------------------------------------------------------------------*/
/* Pour le profiling sur accélérateur                                        */
/*---------------------------------------------------------------------------*/

#if defined(ARCANE_COMPILING_CUDA) && defined(PROF_ACC)

#ifdef P4GPU_HAS_WARNING_INFO
#warning "PROF_ACC : instrumentation avec nvtx"
#endif
#include <nvtx3/nvToolsExt.h>
#include <cuda_profiler_api.h>

#ifndef PROF_ACC_START_CAPTURE
#define PROF_ACC_START_CAPTURE cudaProfilerStart()
#endif

#ifndef PROF_ACC_STOP_CAPTURE
#define PROF_ACC_STOP_CAPTURE cudaProfilerStop()
#endif

#ifndef PROF_ACC_BEGIN
#define PROF_ACC_BEGIN(__name__) nvtxRangePushA(__name__)
#endif

#ifndef PROF_ACC_END
#define PROF_ACC_END nvtxRangePop()
#endif

#elif defined(ARCANE_COMPILING_HIP) && defined(PROF_ACC)

#ifdef P4GPU_HAS_WARNING_INFO
#warning "PROF_ACC : instrumentation avec roctx"
#endif
#include <roctracer/roctx.h>
//#include <roctracer/roctracer_ext.h>

#ifndef PROF_ACC_START_CAPTURE
#define PROF_ACC_START_CAPTURE //roctracer_start()
#endif

#ifndef PROF_ACC_STOP_CAPTURE
#define PROF_ACC_STOP_CAPTURE //roctracer_stop()
#endif

#ifndef PROF_ACC_BEGIN
#define PROF_ACC_BEGIN(__name__) roctxRangePushA(__name__)
#endif

#ifndef PROF_ACC_END
#define PROF_ACC_END roctxRangePop()
#endif

#else

//#warning "Pas d'instrumentation"
#ifndef PROF_ACC_START_CAPTURE
#define PROF_ACC_START_CAPTURE 
#endif

#ifndef PROF_ACC_STOP_CAPTURE
#define PROF_ACC_STOP_CAPTURE 
#endif

#ifndef PROF_ACC_BEGIN
#define PROF_ACC_BEGIN(__name__)
#endif

#ifndef PROF_ACC_END
#define PROF_ACC_END 
#endif

#endif


/*---------------------------------------------------------------------------*/
/* Catégorie de priorités pour création d'une queue                          */
/*---------------------------------------------------------------------------*/
enum eQueuePriority {
  QP_low = 0,  //! jamais prioritaire
  QP_default,  //! avec la même priorité qu'une queue par défaut
  QP_high      //! toujours prioritaire
};

/*---------------------------------------------------------------------------*/
/* Disponibilité d'un accélérateur associé à un runner                       */
/*---------------------------------------------------------------------------*/
class AcceleratorUtils {
 public:
  static bool isAvailable(const ax::Runner& runner) {
    return ax::impl::isAcceleratorPolicy(runner.executionPolicy());
  }

  static Integer deviceCount() {
#if defined(ARCANE_COMPILING_CUDA)
    Integer device_count=0;
    cudaGetDeviceCount(&device_count);
    return device_count;
#elif defined(ARCANE_COMPILING_HIP)
    Integer device_count=0;
    auto err = hipGetDeviceCount(&device_count);
    return device_count;
#else
    return 0;
#endif
  }

  static void setDevice([[maybe_unused]] Integer device) {
#if defined(ARCANE_COMPILING_CUDA)
    cudaSetDevice(device);
#elif defined(ARCANE_COMPILING_HIP)
    auto err = hipSetDevice(device);
#else
    return ;
#endif
  }

  /*---------------------------------------------------------------------------*/
  /* Référence sur une queue asynchrone créée avec un niveau de priorité       */
  /*---------------------------------------------------------------------------*/
  static Ref<ax::RunQueue> refQueueAsync(ax::Runner& runner, eQueuePriority qp) {
    ax::RunQueueBuildInfo bi;
    // 0 = priorité par défaut
    // Plus la valeur de priorité est faible, plus la queue sera prioritaire
    // TODO : récupérer avec Arcane les valeurs [min,max] admissibles (et non plus utiliser +-10)
    if (qp==QP_high) {
      bi.setPriority(-10); 
    } else if (qp==QP_low) {
      bi.setPriority(+10); 
    } // else, par défaut on n'affecte pas de priorité
    auto ref_queue = ax::makeQueueRef(runner, bi);
    ref_queue->setAsync(true);

    return ref_queue;
  }
};

/*---------------------------------------------------------------------------*/
/* Pour gérer un nombre dynamique de RunQueue asynchrones                    */
/*---------------------------------------------------------------------------*/
class MultiAsyncRunQueue {
 public:

  MultiAsyncRunQueue(ax::Runner& runner, Integer asked_nb_queue, 
      bool unlimited=false, eQueuePriority qp=QP_default) {
    if (unlimited) {
      m_nb_queue=asked_nb_queue;
    } else {
      // au plus 32 queues (32 = nb de kernels max exécutables simultanément)
      m_nb_queue = std::min(asked_nb_queue, 32);
    }
    m_queues.reserve(m_nb_queue);
    for(Integer iq=0 ; iq<m_nb_queue ; ++iq) {
      m_queues.add(AcceleratorUtils::refQueueAsync(runner, qp));
    }
  }

  virtual ~MultiAsyncRunQueue() {
    m_nb_queue=0;
    m_queues.clear();
  }

  // Pour récupérer la iq%nbQueue()-ième queue d'exécution
  inline ax::RunQueue& queue(Integer iq) {
    return *(m_queues[iq%m_nb_queue].get());
  }

  // Force l'attente de toutes les RunQueue
  void waitAllQueues() {
    for(Integer iq=0 ; iq<m_nb_queue ; ++iq) {
      m_queues[iq]->barrier();
    }
  }

  // Nombre de RunQueue asynchrones
  inline Integer nbQueue() const {
    ARCANE_ASSERT(m_nb_queue==m_queues.size(), ("m_nb_queue est different de m_queues.size()"));
    return m_nb_queue;
  }

 protected:
  UniqueArray< Ref<ax::RunQueue> > m_queues; //!< toutes les RunQueue
  Integer m_nb_queue=0; //!< m_queues.size()
};

/*---------------------------------------------------------------------------*/
/* "Conseiller" mémoire, permet de caractériser des accès mémoire            */
/*---------------------------------------------------------------------------*/
class AccMemAdviser {
 public:
  AccMemAdviser(bool enable=true) : 
    m_enable(enable) {
#ifdef ARCANE_COMPILING_CUDA
    if (m_enable) {
      cudaGetDevice(&m_device);
    }
#endif
  }
  ~AccMemAdviser() {}

  bool enable() const { return m_enable; }

  template<typename ViewType>
  void setReadMostly([[maybe_unused]] ViewType view) {
#ifdef ARCANE_COMPILING_CUDA
    if (m_enable && view.size()) {
      cudaMemAdvise (view.data(), view.size(), cudaMemAdviseSetReadMostly,m_device);
    }
#endif
  }
 private:
  bool m_enable=true;
  int m_device=-1;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#if ARCANE_VERSION >= 30800
#define MD_Dim1 MDDim1
#define MD_Dim2 MDDim2
#else
#define MD_Dim1 1
#define MD_Dim2 2
#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class Real3_View8
{
 public:
  Real3_View8(NumArray<Real3,MD_Dim2>& v) : m_ptr(v.to1DSpan().data()), nb_cell(v.dim2Size()) {}
 public:
  ARCCORE_HOST_DEVICE Real3& operator()(int node_index,int cell_index) const
  {
    return m_ptr[node_index*nb_cell + cell_index];
  }
 private:
  Real3* m_ptr;
  Int32 nb_cell;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
