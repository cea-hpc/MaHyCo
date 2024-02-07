#ifndef ACCELERATOR_UTILS_H
#define ACCELERATOR_UTILS_H

#include "arcane/IApplication.h"
#include "arcane/core/ISubDomain.h"
#include "arcane/accelerator/Reduce.h"
#include "arcane/accelerator/Runner.h"
#include "arcane/accelerator/VariableViews.h"
#include "arcane/accelerator/Accelerator.h"
#include "arcane/accelerator/RunCommandLoop.h"
#include "arcane/accelerator/RunCommandEnumerate.h"
#include "arcane/accelerator/core/RunQueueBuildInfo.h"
#include "arcane/accelerator/core/Memory.h"
#include "arcane/accelerator/MaterialVariableViews.h"
#include "arcane/accelerator/RunCommandMaterialEnumerate.h"
#include "arcane/accelerator/core/IAcceleratorMng.h"
#include "arcane/accelerator/core/RunQueueEvent.h"

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

  static bool isAvailable(ISubDomain* sd) {
    return ax::impl::isAcceleratorPolicy(sd->acceleratorMng()->defaultRunner()->executionPolicy());
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

  // Vue sur les queues asynchrones
  ArrayView< Ref<ax::RunQueue> > view() {
    return m_queues.view();
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
/* Pour gérer des dépendances entre une queue principale et une liste de     */
/*   queues secondaires                                                      */
/*---------------------------------------------------------------------------*/
class MainAndSecondaryQueuesDependencies
{
 public:
  MainAndSecondaryQueuesDependencies() {}
  ~MainAndSecondaryQueuesDependencies() {}

  /*---------------------------------------------------------------------------------*/
  /* Création d'une dépendance des queues secondaires sur une queue principale 

      m   s0   s1   s2 ..
      |   :    :    :
      |__ :___ :___ :
      :  \:   \:   \: 
      :   |    |    |
      :   |    |    |
      :   |    |    |
   */
  void secQueuesDependOnMainQueue(Runner& runner, ArrayView<Ref<RunQueue>> sec_queues, Ref<RunQueue> main_queue)
  {
    if (!m_are_events_initialized) {
      // Pré-construction des événements
      m_main_event = ax::makeEventRef(runner);
      m_sec_event.resize(sec_queues.size());
      for (Integer isec = 0; isec < m_sec_event.size(); ++isec) {
        m_sec_event[isec] = ax::makeEventRef(runner);
      }

      m_are_events_initialized = true;
    }
    if (sec_queues.size() != m_sec_event.size()) {
      throw FatalErrorException(A_FUNCINFO, "sec_queues.size() != m_sec_event.size()");
    }

    // On exprime la dépendance des queues secondaires à l'événement "racine" m_main_event
    main_queue->recordEvent(m_main_event); // événement "racine"

    for(auto sec_queue : sec_queues)
    {
      // la queue secondaire devra attendre l'occurence de l'événement m_main_queue
      sec_queue->waitEvent(m_main_event);  
    }
  }

  /*---------------------------------------------------------------------------------*/
  /* Création d'une dépendance de la queue principale sur les queues secondaires

      m   s0   s1   s2 ..
      :   |    |    |
      :   |    |    |
      :   |    |    |
      :__/:___/:___/:
      |   :    :    :
      |   :    :    :
   */
  void mainQueueDependsOnSecQueues(Ref<RunQueue> main_queue, ArrayView<Ref<RunQueue>> sec_queues)
  {
    if (!m_are_events_initialized) {
      throw FatalErrorException(
          A_FUNCINFO, "mainQueueDependsOnSecQueues ne peut pas être appelée avant secQueuesDependOnMainQueue");
    }
    if (sec_queues.size() != m_sec_event.size()) {
      throw FatalErrorException(A_FUNCINFO, "sec_queues.size() != m_sec_event.size()");
    }

    // On exprime la dépendance de main_queue aux événements de "fin de calcul pour isec"
    for (Integer isec = 0; isec < m_sec_event.size(); ++isec) 
    {
      // Enregistrement de l'événment "fin de calcul pour isec"
      sec_queues[isec]->recordEvent(m_sec_event[isec]);
      // La queue principale attendra l'occurence de l'événement "fin de calcul pour isec"
      main_queue->waitEvent(m_sec_event[isec]);
    }
  }

 private:
  // Evenements pré-construits pour synchroniser la queue principale et les queues secondaires
  bool m_are_events_initialized = false; //! Passera à vrai quand les événements ci-dessous seront créés
  Ref<ax::RunQueueEvent> m_main_event;
  UniqueArray< Ref<ax::RunQueueEvent> > m_sec_event;
};

/*---------------------------------------------------------------------------*/
/* Collection d'un ensemble de MainAndSecondaryQueuesDependencies            */
/*---------------------------------------------------------------------------*/
class MainAndSecondaryQueuesDependenciesPool
{
 public:
  MainAndSecondaryQueuesDependenciesPool() {}
  ~MainAndSecondaryQueuesDependenciesPool() {}

  //! A appeler avant une utilisation du pool
  void init()
  {
    m_cur_dep=0;
  }

  //! Récupère une référence sur un MainAndSecondaryQueuesDependencies
  // (à appeler entre init() et finalize())
  Ref<MainAndSecondaryQueuesDependencies> popDependenciesRef()
  {
    if (m_cur_dep == m_all_dependencies.size())
    {
      m_all_dependencies.add(makeRef(new MainAndSecondaryQueuesDependencies()));
    }
    auto ref_dependencies = m_all_dependencies[m_cur_dep];
    m_cur_dep++;
    return ref_dependencies;
  }

  //! A appeler après une utilisation du pool
  void finalize()
  {
    // Pour l'instant, pas trop utile mais pourrait servir
    m_cur_dep = 0;
  }

  //! Détruit toutes les références de MainAndSecondaryQueuesDependencies du pool
  void clear()
  {
    m_all_dependencies.clear();
  }

 private:
  Integer m_cur_dep = 0;
  UniqueArray< Ref<MainAndSecondaryQueuesDependencies> > m_all_dependencies;
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
