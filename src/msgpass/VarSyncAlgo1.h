#ifndef MSG_PASS_VAR_SYNC_ALGO1_H
#define MSG_PASS_VAR_SYNC_ALGO1_H

#include "msgpass/IAlgo1SyncData.h"

#include <arcane/IParallelMng.h>

/*---------------------------------------------------------------------------*/
/* \class VarSyncAlgo1                                                       */
/* \brief Algorithm (v1) to synchronize mesh variables                       */
/*---------------------------------------------------------------------------*/
class VarSyncAlgo1 {
 public:
  VarSyncAlgo1(IParallelMng* pm, Int32ConstArrayView neigh_ranks);
  virtual ~VarSyncAlgo1();

  //! Synchronize variables encapsulated into sync_data
  void synchronize(IAlgo1SyncData* sync_data);
 protected:
  IParallelMng* m_pm=nullptr;  //! To perform send/recv
  Int32ConstArrayView m_neigh_ranks;  //! List of neighbour ranks
  Integer m_nb_nei;  //! Number of neighbours (m_neigh_ranks.size())
};

#endif

