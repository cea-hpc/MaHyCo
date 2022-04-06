#ifndef MSG_PASS_MESH_VARIABLE_SYNCHRONIZER_LIST_H
#define MSG_PASS_MESH_VARIABLE_SYNCHRONIZER_LIST_H

#include "accenv/MultiEnvUtils.h"

#include <arcane/materials/MeshMaterialVariable.h>

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

 public:
  IMeshVarSync() {}
  virtual ~IMeshVarSync() {}

  //! Different sizes/properties depending on unit type
  virtual SizeInfos sizeInfos() const = 0;

  //! Pointer to MeshMaterialVariable if it exists (nullptr otherwise)
  virtual MeshMaterialVariable* materialVariable() = 0;

  //! Estimate an upper bound of the buffer size to pack <item_sizes> values
  virtual Int64 estimatedMaxBufSz(IntegerConstArrayView item_sizes) const = 0;

  //! Asynchronously pack "shared" cell (levis) into the buffer (buf)
  virtual void asyncPackIntoBuf(ConstArrayView<EnvVarIndex> levis,
      ArrayView<Byte> buf, RunQueue& queue) = 0;

  //! Asynchronously unpack buffer (buf) into "ghost" cell (levis)
  virtual void asyncUnpackFromBuf(ConstArrayView<EnvVarIndex> levis,
      ArrayView<Byte> buf, RunQueue& queue) = 0;
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

  //! Estimate an upper bound of the buffer size to pack <item_sizes> values
  Int64 estimatedMaxBufSz(IntegerConstArrayView item_sizes) const override;

  //! Asynchronously pack "shared" cell (levis) into the buffer (buf)
  void asyncPackIntoBuf(ConstArrayView<EnvVarIndex> levis,
      ArrayView<Byte> buf, RunQueue& queue) override;

  //! Asynchronously unpack buffer (buf) into "ghost" cell (levis)
  void asyncUnpackFromBuf(ConstArrayView<EnvVarIndex> levis,
      ArrayView<Byte> buf, RunQueue& queue) override;

 protected:
  CellMaterialVariableScalarRef<DataType> m_var;  //! Variable to synchronize
  MultiEnvVarHD<DataType> m_menv_var;  //! View memories on multi-mat data in HOST/DEVICE
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief List of mesh variables to synchronize
 */
class MeshVariableSynchronizerList {
 public:
  MeshVariableSynchronizerList(BufAddrMng* bam);

  virtual ~MeshVariableSynchronizerList();

  //! Add a multi-mat variable into the list of variables to synchronize
  template<typename DataType>
  void add(CellMaterialVariableScalarRef<DataType> var_menv);

  //! Return the list of variables to synchronize
  ConstArrayView<IMeshVarSync*> varsList() const;

  //! Asynchronous pointers tranfer onto device
  void asyncHToD(RunQueue& queue);

 protected:
  BufAddrMng* m_buf_addr_mng=nullptr;
  UniqueArray<IMeshVarSync*> m_vars;  //! List of variables to synchronize
};

#endif

