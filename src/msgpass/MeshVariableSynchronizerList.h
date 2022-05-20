#ifndef MSG_PASS_MESH_VARIABLE_SYNCHRONIZER_LIST_H
#define MSG_PASS_MESH_VARIABLE_SYNCHRONIZER_LIST_H

#include "accenv/MultiEnvUtils.h"
#include "msgpass/SyncItems.h"

#include <arcane/materials/MeshMaterialVariable.h>

class VarSyncMng;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Interface of a mesh variable to synchronize
 */
class IMeshVarSync {
 public:
  struct SizeInfos {
     size_t alignOf;
     size_t sizeOf;
     size_t sizeOfItem;
  };

  //! Type of items for synchronization (shared or ghost)
  enum eItemSync {
    IS_owned = 0,
    IS_ghost
  };

 public:
  IMeshVarSync() {}
  virtual ~IMeshVarSync() {}

  //! Different sizes/properties depending on unit type
  virtual SizeInfos sizeInfos() const = 0;

  //! Pointer to MeshMaterialVariable if it exists (nullptr otherwise)
  virtual MeshMaterialVariable* materialVariable() = 0;

  //! Pointer to IVariable if it exists (nullptr otherwise)
  virtual IVariable* variable() = 0;

  //! Estimate an upper bound of the buffer size to pack <item_sizes> values
  // TODO : à supprimer
  virtual Int64 estimatedMaxBufSz(IntegerConstArrayView item_sizes) const = 0;

  //! Estimate an upper bound of the buffer size to pack/unpack the variable values
  virtual Int64 estimatedMaxBufSz() const = 0;

  //! Space in bytes to store the variable values on the item_sync items for the neighbour inei
  virtual size_t sizeInBytes(eItemSync item_sync, Integer inei) const = 0;

  //! Asynchronously pack "shared" cell (levis) into the buffer (buf)
  // TODO : à supprimer
  virtual void asyncPackIntoBuf(ConstArrayView<EnvVarIndex> levis,
      ArrayView<Byte> buf, RunQueue& queue) = 0;

  //! Asynchronously unpack buffer (buf) into "ghost" cell (levis)
  // TODO : à supprimer
  virtual void asyncUnpackFromBuf(ConstArrayView<EnvVarIndex> levis,
      ArrayView<Byte> buf, RunQueue& queue) = 0;

  //! Asynchronously pack "shared" items with neighbour <inei> into the buffer (buf)
  virtual void asyncPackOwnedIntoBuf(Integer inei, ArrayView<Byte> buf, 
      RunQueue& queue) = 0;

  //! Asynchronously unpack "ghost" items with neighbour <inei> from the buffer (buf)
  virtual void asyncUnpackGhostFromBuf(Integer inei, ArrayView<Byte> buf, 
      RunQueue& queue) = 0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief a Cell multi-mat variable to synchronize
 */
template<typename DataType>
class CellMatVarScalSync : public IMeshVarSync {
 public:
  CellMatVarScalSync(CellMaterialVariableScalarRef<DataType> var,
      BufAddrMng* bam);

  virtual ~CellMatVarScalSync();

  //! Different sizes/properties depending on unit type
  SizeInfos sizeInfos() const override;

  //! Pointer to MeshMaterialVariable if it exists (nullptr otherwise)
  MeshMaterialVariable* materialVariable() override;

  //! Pointer to IVariable if it exists (nullptr otherwise)
  IVariable* variable() override;

  //! Estimate an upper bound of the buffer size to pack <item_sizes> values
  // TODO : a supprimer
  Int64 estimatedMaxBufSz(IntegerConstArrayView item_sizes) const override;

  //! Estimate an upper bound of the buffer size to pack/unpack the variable values
  // TODO : à implémenter
  Int64 estimatedMaxBufSz() const override;

  //! Space in bytes to store the variable values on the item_sync items for the neighbour inei
  // TODO : à implémenter
  size_t sizeInBytes(eItemSync item_sync, Integer inei) const override;

  //! Asynchronously pack "shared" cell (levis) into the buffer (buf)
  // TODO : à supprimer
  void asyncPackIntoBuf(ConstArrayView<EnvVarIndex> levis,
      ArrayView<Byte> buf, RunQueue& queue) override;

  //! Asynchronously unpack buffer (buf) into "ghost" cell (levis)
  // TODO : à supprimer
  void asyncUnpackFromBuf(ConstArrayView<EnvVarIndex> levis,
      ArrayView<Byte> buf, RunQueue& queue) override;

  //! Asynchronously pack "shared" items with neighbour <inei> into the buffer (buf)
  // TODO : à implémenter
  void asyncPackOwnedIntoBuf(Integer inei, ArrayView<Byte> buf, 
      RunQueue& queue) override;

  //! Asynchronously unpack "ghost" items with neighbour <inei> from the buffer (buf)
  // TODO : à implémenter
  void asyncUnpackGhostFromBuf(Integer inei, ArrayView<Byte> buf, 
      RunQueue& queue) override;

 protected:
  CellMaterialVariableScalarRef<DataType> m_var;  //! Variable to synchronize
  MultiEnvVarHD<DataType> m_menv_var;  //! View memories on multi-mat data in HOST/DEVICE
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief a global variable to synchronize
 */
template<typename MeshVariableRefT>
class GlobVarSync : public IMeshVarSync {
  using DataType = typename MeshVariableRefT::DataType;
  using ItemType = typename MeshVariableRefT::ItemType;
 public:
  GlobVarSync(MeshVariableRefT var, 
      SyncItems<typename MeshVariableRefT::ItemType>* sync_items);

  virtual ~GlobVarSync();

  //! Different sizes/properties depending on unit type
  SizeInfos sizeInfos() const override;

  //! Pointer to MeshMaterialVariable if it exists (nullptr otherwise)
  MeshMaterialVariable* materialVariable() override;

  //! Pointer to IVariable if it exists (nullptr otherwise)
  IVariable* variable() override;

  //! Estimate an upper bound of the buffer size to pack <item_sizes> values
  // TODO : à supprimer
  Int64 estimatedMaxBufSz(IntegerConstArrayView item_sizes) const override;

  //! Estimate an upper bound of the buffer size to pack/unpack the variable values
  Int64 estimatedMaxBufSz() const override;

  //! Space in bytes to store the variable values on the item_sync items for the neighbour inei
  size_t sizeInBytes(eItemSync item_sync, Integer inei) const override;

  //! Asynchronously pack "shared" cell (levis) into the buffer (buf)
  // TODO : à supprimer
  void asyncPackIntoBuf(ConstArrayView<EnvVarIndex> levis,
      ArrayView<Byte> buf, RunQueue& queue) override;

  //! Asynchronously unpack buffer (buf) into "ghost" cell (levis)
  // TODO : à supprimer
  void asyncUnpackFromBuf(ConstArrayView<EnvVarIndex> levis,
      ArrayView<Byte> buf, RunQueue& queue) override;

  //! Asynchronously pack "shared" items with neighbour <inei> into the buffer (buf)
  void asyncPackOwnedIntoBuf(Integer inei, ArrayView<Byte> buf, 
      RunQueue& queue) override;

  //! Asynchronously unpack "ghost" items with neighbour <inei> from the buffer (buf)
  void asyncUnpackGhostFromBuf(Integer inei, ArrayView<Byte> buf, 
      RunQueue& queue) override;

 protected:
  MeshVariableRefT m_var;  //! Variable to synchronize
  SyncItems<ItemType>* m_sync_items;  //! Items to synchronize
  Integer m_degree;  //! For a given ItemType, how many are DataType used ?
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief List of mesh variables to synchronize
 */
class MeshVariableSynchronizerList {
 public:
  MeshVariableSynchronizerList(VarSyncMng* vsync_mng);

  virtual ~MeshVariableSynchronizerList();

  //! Add a multi-mat variable into the list of variables to synchronize
  template<typename DataType>
  void add(CellMaterialVariableScalarRef<DataType> var_menv);

  //! Add a global variable into the list of variables to synchronize
  template<typename MeshVariableRefT>
  void add(MeshVariableRefT var);

  //! Return the list of variables to synchronize
  ConstArrayView<IMeshVarSync*> varsList() const;

  //! Asynchronous pointers tranfer onto device
  void asyncHToD(RunQueue& queue);

 protected:
  VarSyncMng* m_vsync_mng=nullptr;
  BufAddrMng* m_buf_addr_mng=nullptr;
  UniqueArray<IMeshVarSync*> m_vars;  //! List of variables to synchronize
};

#endif

