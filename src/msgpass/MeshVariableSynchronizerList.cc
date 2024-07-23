// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#include "msgpass/MeshVariableSynchronizerList.h"
#include "msgpass/SyncBuffers.h"
#include "msgpass/PackTransfer.h"
#include "msgpass/VarSyncMng.h"

/*---------------------------------------------------------------------------*/
/* CellMatVarScalSync<DataType> : a Cell multi-mat variable to synchronize   */
/*---------------------------------------------------------------------------*/
template<typename DataType>
CellMatVarScalSync<DataType>::CellMatVarScalSync(
    CellMaterialVariableScalarRef<DataType> var,
    SyncEnvIndexes* sync_evi) : 
  IMeshVarSync(),
  m_var (var),
  m_sync_evi (sync_evi)
{
}

template<typename DataType>
CellMatVarScalSync<DataType>::~CellMatVarScalSync() {
}

//! Different sizes/properties depending on unit type
template<typename DataType>
IMeshVarSync::SizeInfos CellMatVarScalSync<DataType>::sizeInfos() const {
  return {alignof(DataType), sizeof(DataType), sizeof(DataType)*1};
}

//! Pointer to MeshMaterialVariable if it exists (nullptr otherwise)
template<typename DataType>
MeshMaterialVariable* CellMatVarScalSync<DataType>::materialVariable() {
  IMeshMaterialVariable* mvar = m_var.materialVariable();
  return dynamic_cast<MeshMaterialVariable*>(mvar);
}

//! Pointer to IVariable if it exists (nullptr otherwise)
template<typename DataType>
IVariable* CellMatVarScalSync<DataType>::variable() {
  return nullptr;
}

//! Estimate an upper bound of the buffer size to pack/unpack the variable values
template<typename DataType>
Int64 CellMatVarScalSync<DataType>::estimatedMaxBufSz() const
{
  auto nb_owned_mvi_pn = m_sync_evi->nbOwnedMviPn();
  auto nb_ghost_mvi_pn = m_sync_evi->nbGhostMviPn();

  Int64 buf_estim_sz=0;
  buf_estim_sz += SyncBuffers::estimatedMaxBufSz<DataType>(nb_owned_mvi_pn, /*degree=*/1);
  buf_estim_sz += SyncBuffers::estimatedMaxBufSz<DataType>(nb_ghost_mvi_pn, /*degree=*/1);

  return buf_estim_sz;
}

//! Space in bytes to store the variable values on the item_sync items for the neighbour inei
template<typename DataType>
size_t CellMatVarScalSync<DataType>::sizeInBytes(eItemSync item_sync, Integer inei) const
{
  auto item_sizes = (item_sync == IS_owned ? 
      m_sync_evi->nbOwnedMviPn() :
      m_sync_evi->nbGhostMviPn());

  size_t sizeof_item = sizeof(DataType)*1;  // 1 = degree
  size_t sz_nei_in_bytes = item_sizes[inei]*sizeof_item;

  return sz_nei_in_bytes;
}

//! Asynchronously pack "shared" items with neighbour <inei> into the buffer (buf)
template<typename DataType>
void CellMatVarScalSync<DataType>::asyncPackOwnedIntoBuf(
    Integer inei,
    ArrayView<Byte> buf, RunQueue& queue)  
{
  ConstArrayView<MatVarIndex> lmvis = m_sync_evi->ownedMviPn()[inei];  // Owned
  auto command = makeCommand(queue);

  auto in_var = ax::viewIn(command, m_var);
  Span<const MatVarIndex> in_lmvis(lmvis);
  Span<DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_mvis = lmvis.size();

  command << RUNCOMMAND_LOOP1(iter, nb_mvis) {
    auto [i] = iter();
    buf_vals[i] = in_var[ ComponentItemLocalId(in_lmvis[i]) ];
  }; // asynchronous
}

//! Asynchronously unpack "ghost" items with neighbour <inei> from the buffer (buf)
template<typename DataType>
void CellMatVarScalSync<DataType>::asyncUnpackGhostFromBuf(
    Integer inei,
    ArrayView<Byte> buf, RunQueue& queue)  
{
  ConstArrayView<MatVarIndex> lmvis = m_sync_evi->ghostMviPn()[inei];  // Ghost
  auto command = makeCommand(queue);

  auto out_var = ax::viewOut(command, m_var);
  Span<const MatVarIndex> in_lmvis(lmvis);
  Span<const DataType>    buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_mvis = lmvis.size();

  command << RUNCOMMAND_LOOP1(iter, nb_mvis) {
    auto [i] = iter();
    out_var[ComponentItemLocalId(in_lmvis[i])] = buf_vals[i];
  }; // asynchrone
}

/*---------------------------------------------------------------------------*/
/* GlobVarSync<MeshVariableRefT> : a global variable to synchronize          */
/*---------------------------------------------------------------------------*/
template<typename MeshVariableRefT>
GlobVarSync<MeshVariableRefT>::GlobVarSync(
    MeshVariableRefT var, 
    SyncItems<typename MeshVariableRefT::ItemType>* sync_items) : 
  IMeshVarSync(),
  m_var (var),
  m_sync_items (sync_items)
{
  // For a given ItemType, how many are DataType used ? => m_degree
  m_degree = get_var_degree(m_var);
}

template<typename MeshVariableRefT>
GlobVarSync<MeshVariableRefT>::~GlobVarSync() {
}

//! Different sizes/properties depending on unit type
template<typename MeshVariableRefT>
IMeshVarSync::SizeInfos GlobVarSync<MeshVariableRefT>::sizeInfos() const {
  return {alignof(DataType), sizeof(DataType), sizeof(DataType)*m_degree};
}

//! Pointer to MeshMaterialVariable if it exists (nullptr otherwise)
template<typename MeshVariableRefT>
MeshMaterialVariable* GlobVarSync<MeshVariableRefT>::materialVariable() {
  return nullptr;
}

//! Pointer to IVariable if it exists (nullptr otherwise)
template<typename MeshVariableRefT>
IVariable* GlobVarSync<MeshVariableRefT>::variable() {
  return m_var.variable();
}

//! Estimate an upper bound of the buffer size to pack/unpack the variable values
template<typename MeshVariableRefT>
Int64 GlobVarSync<MeshVariableRefT>::estimatedMaxBufSz() const
{
  auto nb_owned_item_idx_pn = m_sync_items->nbOwnedItemIdxPn();
  auto nb_ghost_item_idx_pn = m_sync_items->nbGhostItemIdxPn();

  Int64 buf_estim_sz=0;
  buf_estim_sz += SyncBuffers::estimatedMaxBufSz<DataType>(nb_owned_item_idx_pn, m_degree);
  buf_estim_sz += SyncBuffers::estimatedMaxBufSz<DataType>(nb_ghost_item_idx_pn, m_degree);

  return buf_estim_sz;
}

//! Space in bytes to store the variable values on the item_sync items for the neighbour inei
template<typename MeshVariableRefT>
size_t GlobVarSync<MeshVariableRefT>::sizeInBytes(eItemSync item_sync, Integer inei) const
{
  auto item_sizes = (item_sync == IS_owned ? 
      m_sync_items->nbOwnedItemIdxPn() :
      m_sync_items->nbGhostItemIdxPn());

  size_t sizeof_item = sizeof(DataType)*m_degree;
  size_t sz_nei_in_bytes = item_sizes[inei]*sizeof_item;

  return sz_nei_in_bytes;
}

//! Asynchronously pack "shared" items with neighbour <inei> into the buffer (buf)
template<typename MeshVariableRefT>
void GlobVarSync<MeshVariableRefT>::asyncPackOwnedIntoBuf(
    Integer inei,
    ArrayView<Byte> buf, RunQueue& queue)  
{
  auto owned_item_idx = m_sync_items->ownedItemIdxPn()[inei];
  async_pack_var2buf(owned_item_idx, m_var, buf, queue);
}

//! Asynchronously unpack "ghost" items with neighbour <inei> from the buffer (buf)
template<typename MeshVariableRefT>
void GlobVarSync<MeshVariableRefT>::asyncUnpackGhostFromBuf(
    Integer inei,
    ArrayView<Byte> buf, RunQueue& queue)  
{
  auto ghost_item_idx = m_sync_items->ghostItemIdxPn()[inei];
  async_unpack_buf2var(ghost_item_idx, buf, m_var, queue);
}

/*---------------------------------------------------------------------------*/
/* MeshVariableSynchronizerList : List of mesh variables to synchronize      */
/*---------------------------------------------------------------------------*/
MeshVariableSynchronizerList::MeshVariableSynchronizerList(VarSyncMng* vsync_mng) :
  m_vsync_mng (vsync_mng)
{
  m_buf_addr_mng = m_vsync_mng->bufAddrMng();
}

MeshVariableSynchronizerList::~MeshVariableSynchronizerList() {
  for(auto v : m_vars) {
    delete v;
  }
  m_buf_addr_mng->reset();
}

//! Add a multi-mat variable into the list of variables to synchronize
template<typename DataType>
void MeshVariableSynchronizerList::add(CellMaterialVariableScalarRef<DataType> var_menv) {
  auto sync_evi = m_vsync_mng->syncEnvIndexes();
  m_vars.add(new CellMatVarScalSync<DataType>(var_menv, sync_evi));
}

//! Add a global variable into the list of variables to synchronize
template<typename MeshVariableRefT>
void MeshVariableSynchronizerList::add(MeshVariableRefT var) {
  using ItemType = typename MeshVariableRefT::ItemType;
  SyncItems<ItemType>* sync_items = m_vsync_mng->getSyncItems<ItemType>();

  m_vars.add(new GlobVarSync<MeshVariableRefT>(var, sync_items));
}

//! Return the list of variables to synchronize
ConstArrayView<IMeshVarSync*> MeshVariableSynchronizerList::varsList() const {
  return m_vars;
}

//! Asynchronous pointers tranfer onto device
void MeshVariableSynchronizerList::asyncHToD(RunQueue& queue) {
  m_buf_addr_mng->asyncCpyHToD(queue);
}

/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/
#include <arcane/utils/Real3x3.h>

#define INST_MESH_VAR_SYNC_LIST_ADD(__DataType__) \
  template void MeshVariableSynchronizerList::add(CellMaterialVariableScalarRef<__DataType__> var_menv)

INST_MESH_VAR_SYNC_LIST_ADD(Byte);
INST_MESH_VAR_SYNC_LIST_ADD(Int16);
INST_MESH_VAR_SYNC_LIST_ADD(Int32);
INST_MESH_VAR_SYNC_LIST_ADD(Int64);
INST_MESH_VAR_SYNC_LIST_ADD(Real);
INST_MESH_VAR_SYNC_LIST_ADD(Real2);
INST_MESH_VAR_SYNC_LIST_ADD(Real2x2);
INST_MESH_VAR_SYNC_LIST_ADD(Real3);
INST_MESH_VAR_SYNC_LIST_ADD(Real3x3);

#define INST_MESH_VAR_SYNC_LIST_ADDG(__MeshVariableRefT__) \
  template void MeshVariableSynchronizerList::add(__MeshVariableRefT__ var)

/*------------------------- Cell --------------------------------------------*/

// Scalar
INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellByte);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellInteger);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellInt64);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellReal);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellReal3);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellReal3x3);

// Array
INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellArrayByte);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellArrayInteger);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellArrayInt64);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellArrayReal);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellArrayReal3);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableCellArrayReal3x3);

/*------------------------- Node --------------------------------------------*/

// Scalar
INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeByte);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeInteger);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeInt64);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeReal);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeReal3);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeReal3x3);

// Array
INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeArrayByte);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeArrayInteger);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeArrayInt64);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeArrayReal);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeArrayReal3);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableNodeArrayReal3x3);

/*------------------------- Face --------------------------------------------*/

// Scalar
INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceByte);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceInteger);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceInt64);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceReal);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceReal3);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceReal3x3);

// Array
INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceArrayByte);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceArrayInteger);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceArrayInt64);

INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceArrayReal);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceArrayReal3);
INST_MESH_VAR_SYNC_LIST_ADDG(VariableFaceArrayReal3x3);
