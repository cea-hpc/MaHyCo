#ifndef MSG_PASS_I_ALGO1_SYNC_DATA_H
#define MSG_PASS_I_ALGO1_SYNC_DATA_H

#include <arcane/utils/ArrayView.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Class interface to manage data used by VarSyncAlgo1::synchronize          */
/*---------------------------------------------------------------------------*/
class IAlgo1SyncData {
 public:

  //! Destructor to override
  virtual ~IAlgo1SyncData() {}

  //! True if there is no data to synchronize
  virtual bool isEmpty() = 0;

  //! Initialization step before beginning communications
  virtual void initComm() = 0;

  //! Get the receive buffer for the neighbour inei
  virtual ArrayView<Byte> recvBuf(Integer inei) = 0;

  //! Get the send buffer for the neighbour inei
  virtual ArrayView<Byte> sendBuf(Integer inei) = 0;

  //! Prepare the sendings for every neighbour
  virtual void initSendings() = 0;

  //! Prepare the sending for one neighbour
  virtual void finalizePackBeforeSend(Integer inei) = 0;

  //! Finalize the sendings for every neighbour
  virtual void finalizeSendings() = 0;

  //! Data treatement after receiving the message from the neighbour inei
  virtual void unpackAfterRecv(Integer inei) = 0;

  //! Finalize the receipts for every neighbour
  virtual void finalizeReceipts() = 0;

  //! Finalize if no communication needed 
  virtual void finalizeWoComm() = 0;
};

#endif
