// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#ifndef MSG_PASS_ALGO1_SYNC_DATA_D_H
#define MSG_PASS_ALGO1_SYNC_DATA_D_H

#include "msgpass/IAlgo1SyncData.h"
#include "msgpass/MeshVariableSynchronizerList.h"
#include "msgpass/SyncItems.h"
#include "msgpass/SyncBuffers.h"

/*---------------------------------------------------------------------------*/
/* \class Algo1SyncDataD                                                     */
/* \brief Implementation of IAlgo1SyncData for all variable types            */
/*   packing/unpacking on Device                                             */
/*   communicating (MPI) on Host                                             */
/*---------------------------------------------------------------------------*/
class Algo1SyncDataD : public IAlgo1SyncData {
 public:
  /*!
   * \brief Persistent information which don't depend on the variables
   */
  class PersistentInfo {
    friend class Algo1SyncDataD;
   public:
    PersistentInfo(bool is_device_aware,
	Integer nb_nei,
        Runner& runner,
        SyncBuffers* sync_buffers);
    virtual ~PersistentInfo();
   protected:
    SyncBuffers* m_sync_buffers=nullptr;
    Integer m_nb_nei=0;
    bool m_is_device_aware=false;

    UniqueArray<Ref<ax::RunQueueEvent>> m_pack_events;  //! Les evenements pour le packing des données
  };

 public:
  Algo1SyncDataD(MeshVariableSynchronizerList& vars,
      Ref<RunQueue> ref_queue,
      PersistentInfo& pi);

  virtual ~Algo1SyncDataD();

  //! True if there is no data to synchronize
  bool isEmpty() override;

  //! Initialization step before beginning communications
  void initComm() override;

  //! Get the receive buffer for the neighbour inei
  ArrayView<Byte> recvBuf(Integer inei) override;

  //! Get the send buffer for the neighbour inei
  ArrayView<Byte> sendBuf(Integer inei) override;

  //! Prepare the sendings for every neighbour
  void initSendings() override;

  //! Prepare the sending for one neighbour
  void finalizePackBeforeSend(Integer inei) override;

  //! Finalize the sendings for every neighbour
  void finalizeSendings() override;

  //! Data treatement after receiving the message from the neighbour inei
  void unpackAfterRecv(Integer inei) override;

  //! Finalize the receipts for every neighbour
  void finalizeReceipts() override;

  //! Finalize if no communication needed
  void finalizeWoComm() override;

 protected:
  MeshVariableSynchronizerList& m_vars;
  Ref<RunQueue> m_ref_queue;
  PersistentInfo& m_pi;

  MultiBufView2 m_buf_snd_d;  //! Buffers on Device (_d) to send
  MultiBufView2 m_buf_rcv_d;  //! Buffers on Device (_d) to recv
};

#endif

