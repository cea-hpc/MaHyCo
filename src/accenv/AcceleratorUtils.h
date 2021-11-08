#ifndef ACCELERATOR_UTILS_H
#define ACCELERATOR_UTILS_H

#include "arcane/IApplication.h"
#include "arcane/accelerator/Reduce.h"
#include "arcane/accelerator/Runner.h"
#include "arcane/accelerator/Views.h"
#include "arcane/accelerator/Accelerator.h"
#include "arcane/accelerator/RunCommandLoop.h"
#include "arcane/accelerator/RunCommandEnumerate.h"

/*---------------------------------------------------------------------------*/
/* Pour les accélérateurs                                                    */
/*---------------------------------------------------------------------------*/

using namespace Arcane;
namespace ax = Arcane::Accelerator;

/*---------------------------------------------------------------------------*/
/* Pour le profiling sur accélérateur                                        */
/*---------------------------------------------------------------------------*/

#if defined(ARCANE_COMPILING_CUDA) && defined(PROF_ACC)
#define USE_PROF_ACC
#endif

#if defined(USE_PROF_ACC)

#warning "PROF_ACC : instrumentation avec nvtx"
#include <nvtx3/nvToolsExt.h>

#ifndef PROF_ACC_BEGIN
#define PROF_ACC_BEGIN(__name__) nvtxRangePushA(__name__)
#endif

#ifndef PROF_ACC_END
#define PROF_ACC_END nvtxRangePop()
#endif

#else

//#warning "Pas d'instrumentation"
#ifndef PROF_ACC_BEGIN
#define PROF_ACC_BEGIN(__name__)
#endif

#ifndef PROF_ACC_END
#define PROF_ACC_END 
#endif

#endif

/*---------------------------------------------------------------------------*/
/* Pour gérer un nombre dynamique de RunQueue asynchrones                    */
/*---------------------------------------------------------------------------*/
class MultiAsyncRunQueue {
 public:

  MultiAsyncRunQueue(ax::Runner& runner, Integer asked_nb_queue) {
    // au plus 32 queues (32 = nb de kernels max exécutables simultanément)
    m_nb_queue = std::min(asked_nb_queue, 32);
    m_queues.resize(m_nb_queue);
    for(Integer iq=0 ; iq<m_nb_queue ; ++iq) {
      m_queues[iq] = new ax::RunQueue(runner);
      m_queues[iq]->setAsync(true);
    }
  }

  virtual ~MultiAsyncRunQueue() {
    for(auto q : m_queues) {
      delete q;
    }
    m_nb_queue=0;
    m_queues.clear();
  }

  // Pour récupérer la iq%nbQueue()-ième queue d'exécution
  inline ax::RunQueue& queue(Integer iq) {
    return *(m_queues[iq%m_nb_queue]);
  }

  // Force l'attente de toutes les RunQueue
  void waitAllQueues() {
    for(auto q : m_queues) {
      q->barrier();
    }
  }

  // Nombre de RunQueue asynchrones
  inline Integer nbQueue() const {
    ARCANE_ASSERT(m_nb_queue==m_queues.size(), ("m_nb_queue est different de m_queues.size()"));
    return m_nb_queue;
  }

 protected:
  UniqueArray<ax::RunQueue*> m_queues; //!< toutes les RunQueue
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
  void setReadMostly(ViewType view) {
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


#endif
