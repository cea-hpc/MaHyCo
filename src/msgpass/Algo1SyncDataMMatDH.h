#ifndef MSG_PASS_ALGO1_SYNC_DATA_MMAT_DH_H
#define MSG_PASS_ALGO1_SYNC_DATA_MMAT_DH_H

#include "msgpass/IAlgo1SyncData.h"
#include "msgpass/MeshVariableSynchronizerList.h"
#include "msgpass/SyncEnvIndexes.h"
#include "msgpass/SyncBuffers.h"

/*---------------------------------------------------------------------------*/
/* \class Algo1SyncDataMMatDH                                                */
/* \brief Implementation of IAlgo1SyncData for multi-material variables      */
/*   packing/unpacking on Device                                             */
/*   communicating (MPI) on Host                                             */
/*---------------------------------------------------------------------------*/
class Algo1SyncDataMMatDH : public IAlgo1SyncData {
 public:
  /*!
   * \brief Persistent information which don't depend on the variables
   */
  class PersistentInfo {
    friend class Algo1SyncDataMMatDH;
   public:
    PersistentInfo(Integer nb_nei,
        Runner& runner,
        SyncEnvIndexes* sync_evi,
        SyncBuffers* sync_buffers);
    virtual ~PersistentInfo();
   protected:
    SyncEnvIndexes* m_sync_evi=nullptr;
    SyncBuffers* m_sync_buffers=nullptr;
    Integer m_nb_nei=0;

    Ref<ax::RunQueue> m_ref_queue_data;  //! Référence sur une queue prioritaire pour le transfert des données
    UniqueArray<Ref<ax::RunQueueEvent>> m_pack_events;  //! Les evenements pour le packing des données
    UniqueArray<Ref<ax::RunQueueEvent>> m_transfer_events;  //! Les evenements pour le transfert des données
  };

 public:
  Algo1SyncDataMMatDH(MeshVariableSynchronizerList& vars,
      Ref<RunQueue> ref_queue,
      PersistentInfo& pi);

  virtual ~Algo1SyncDataMMatDH();

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

 protected:
  MeshVariableSynchronizerList& m_vars;
  Ref<RunQueue> m_ref_queue;
  PersistentInfo& m_pi;

  ConstMultiArray2View<EnvVarIndex> m_ghost_evi_pn;  //! Ghost EnvIndex(es) to update

  MultiBufView2 m_buf_snd_h;  //! Buffers on Host (_h) to send
  MultiBufView2 m_buf_rcv_h;  //! Buffers on Host (_h) to recv
  MultiBufView2 m_buf_snd_d;  //! Buffers on Device (_d) to send
  MultiBufView2 m_buf_rcv_d;  //! Buffers on Device (_d) to recv
};

#endif

