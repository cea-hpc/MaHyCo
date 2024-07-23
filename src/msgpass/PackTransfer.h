// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#ifndef MSG_PASS_PACK_TRANSFER_H
#define MSG_PASS_PACK_TRANSFER_H

#include "accenv/AcceleratorUtils.h"
#include "accenv/MultiEnvUtils.h"
#include "msgpass/SyncBuffers.h"

using namespace Arcane;

void async_transfer(Span<Byte> dst_buf, Span<const Byte> src_buf, RunQueue& queue);

void async_transfer(MultiBufView out_buf, MultiBufView in_buf, RunQueue& queue);

/*---------------------------------------------------------------------------*/
/* async_pack_var2buf */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void async_pack_var2buf(IntegerConstArrayView item_idx, 
    const MeshVarRefT<ItemType, DataType>& var,
    ArrayView<Byte> buf, RunQueue& queue) 
{
  // Ne devrait jamais être appelé
  throw NotSupportedException(A_FUNCINFO, String("async_pack_var2buf à spécialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void async_pack_var2buf(IntegerConstArrayView item_idx, 
    const MeshVariableScalarRefT<ItemType, DataType>& var,
    ArrayView<Byte> buf, RunQueue& queue) 
{
  using ItemIdType = typename ItemType::LocalIdType;

  auto command = makeCommand(queue);

  auto in_var = ax::viewIn(command, var);

  Span<const Integer> in_item_idx(item_idx);
  Integer nb_item_idx = item_idx.size();

  Span<DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  command << RUNCOMMAND_LOOP1(iter, nb_item_idx) {
    auto [i] = iter();
    ItemIdType lid(in_item_idx[i]); 
    buf_vals[i] = in_var[lid];
  }; // asynchrone
}

// Spécialisation pour MeshVariable***Array***RefT
template<typename ItemType, typename DataType>
void async_pack_var2buf(IntegerConstArrayView item_idx, 
    const MeshVariableArrayRefT<ItemType, DataType>& var,
    ArrayView<Byte> buf, RunQueue& queue) 
{
  using ItemIdType = typename ItemType::LocalIdType;

  auto command = makeCommand(queue);

  auto in_var = ax::viewIn(command, var);

  Span<const Integer> in_item_idx(item_idx);
  Integer nb_item_idx = item_idx.size();

  Integer degree = var.arraySize();
  // Vue sur tableau 2D [nb_item][degree]
  Span2<DataType> buf_vals(MultiBufView::valBuf2<DataType>(buf, degree));

  command << RUNCOMMAND_LOOP1(iter, nb_item_idx) {
    auto [i] = iter();
    ItemIdType lid(in_item_idx[i]); 

    Span<const DataType> in_var_lid  (in_var[lid]);
    Span<DataType>       buf_vals_i  (buf_vals[i]);
    for(Integer j=0 ; j<degree ; ++j) {
      buf_vals_i[j] = in_var_lid[j];
    }
  }; // asynchrone
}

/*---------------------------------------------------------------------------*/
/* async_unpack_buf2var */
/*---------------------------------------------------------------------------*/
template<typename ItemType, typename DataType, template<typename, typename> class MeshVarRefT>
void async_unpack_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVarRefT<ItemType, DataType> &var, RunQueue& queue) 
{
  // Ne devrait jamais être appelé
  throw NotSupportedException(A_FUNCINFO, String("async_unpack_buf2var à spécialiser"));
}

// Spécialisation pour MeshVariable***Scalar***RefT
template<typename ItemType, typename DataType>
void async_unpack_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVariableScalarRefT<ItemType, DataType> &var, RunQueue& queue)
{
  using ItemIdType = typename ItemType::LocalIdType;

  auto command = makeCommand(queue);

  Span<const Integer>  in_item_idx(item_idx);
  Integer nb_item_idx = item_idx.size();

  Span<const DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  auto out_var = ax::viewOut(command, var);

  command << RUNCOMMAND_LOOP1(iter, nb_item_idx) {
    auto [i] = iter();
    ItemIdType lid(in_item_idx[i]);
    out_var[lid] = buf_vals[i];
  }; // asynchrone
}

// Spécialisation pour MeshVariable***Array***RefT
template<typename ItemType, typename DataType>
void async_unpack_buf2var(IntegerConstArrayView item_idx, 
    ArrayView<Byte> buf,
    MeshVariableArrayRefT<ItemType, DataType> &var, RunQueue& queue)
{
  using ItemIdType = typename ItemType::LocalIdType;

  auto command = makeCommand(queue);

  Span<const Integer>  in_item_idx(item_idx);
  Integer nb_item_idx = item_idx.size();

  Integer degree = var.arraySize();
  // Vue sur tableau 2D [nb_item][degree]
  Span2<const DataType> buf_vals(MultiBufView::valBuf2<DataType>(buf, degree));

  auto out_var = ax::viewOut(command, var);

  command << RUNCOMMAND_LOOP1(iter, nb_item_idx) {
    auto [i] = iter();
    ItemIdType lid(in_item_idx[i]);

    out_var[lid].copy(buf_vals[i].smallView());
  }; // asynchrone
}


/*---------------------------------------------------------------------------*/
/* async_pack_varmenv2buf */
/*---------------------------------------------------------------------------*/
template< typename DataType>
void async_pack_varmenv2buf(ConstArrayView<EnvVarIndex> levis, 
    MultiEnvData<DataType>& in_var_menv,
    ArrayView<Byte> buf, RunQueue& queue) 
{
  auto command = makeCommand(queue);

  Span<const EnvVarIndex> in_levis(levis);
  Span<DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_evis = levis.size();

  command << RUNCOMMAND_LOOP1(iter, nb_evis) {
    auto [i] = iter();
    buf_vals[i] = in_var_menv[ in_levis[i] ];
  }; // asynchrone
}


/*---------------------------------------------------------------------------*/
/* async_unpack_buf2varmenv */
/*---------------------------------------------------------------------------*/
template< typename DataType>
void async_unpack_buf2varmenv(ConstArrayView<EnvVarIndex> levis,
    ArrayView<Byte> buf, 
    MultiEnvData<DataType>& out_var_menv, RunQueue& queue) 
{
  auto command = makeCommand(queue);

  Span<const EnvVarIndex> in_levis(levis);
  Span<const DataType>    buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_evis = levis.size();

  command << RUNCOMMAND_LOOP1(iter, nb_evis) {
    auto [i] = iter();
    out_var_menv.setValue(in_levis[i], buf_vals[i]);
  }; // asynchrone
}
#endif

