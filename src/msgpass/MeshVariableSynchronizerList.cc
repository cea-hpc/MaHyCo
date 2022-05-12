#include "msgpass/MeshVariableSynchronizerList.h"
#include "msgpass/SyncBuffers.h"

/*---------------------------------------------------------------------------*/
/* CellMatVarScalSync<DataType> : a Cell multi-mat variable to synchronize   */
/*---------------------------------------------------------------------------*/
template<typename DataType>
CellMatVarScalSync<DataType>::CellMatVarScalSync(
    CellMaterialVariableScalarRef<DataType> var,
    BufAddrMng* bam) : 
  IMeshVarSync(),
  m_var (var),
  m_menv_var(var, bam)
{
  // m_menv_var on DEVICE will be update by bam
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

//! Estimate an upper bound of the buffer size to pack <item_sizes> values
template<typename DataType>
Int64 CellMatVarScalSync<DataType>::estimatedMaxBufSz(IntegerConstArrayView item_sizes) const
{
  return SyncBuffers::estimatedMaxBufSz<DataType>(item_sizes, /*degree=*/1);
}

//! Asynchronously pack "shared" cell (levis) into the buffer (buf)
template<typename DataType>
void CellMatVarScalSync<DataType>::asyncPackIntoBuf(
    ConstArrayView<EnvVarIndex> levis,
    ArrayView<Byte> buf, RunQueue& queue)  
{
  auto command = makeCommand(queue);

  auto in_var_menv = m_menv_var.spanD();
  Span<const EnvVarIndex> in_levis(levis);
  Span<DataType> buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_evis = levis.size();

  command << RUNCOMMAND_LOOP1(iter, nb_evis) {
    auto [i] = iter();
    buf_vals[i] = in_var_menv[ in_levis[i] ];
  }; // asynchronous
}

//! Asynchronously unpack buffer (buf) into "ghost" cell (levis)
template<typename DataType>
void CellMatVarScalSync<DataType>::asyncUnpackFromBuf(
    ConstArrayView<EnvVarIndex> levis,
    ArrayView<Byte> buf, RunQueue& queue) 
{
  auto command = makeCommand(queue);

  auto out_var_menv = m_menv_var.spanD();
  Span<const EnvVarIndex> in_levis(levis);
  Span<const DataType>    buf_vals(MultiBufView::valBuf<DataType>(buf));

  Integer nb_evis = levis.size();

  command << RUNCOMMAND_LOOP1(iter, nb_evis) {
    auto [i] = iter();
    out_var_menv.setValue(in_levis[i], buf_vals[i]);
  }; // asynchrone
}

/*---------------------------------------------------------------------------*/
/* MeshVariableSynchronizerList : List of mesh variables to synchronize      */
/*---------------------------------------------------------------------------*/
MeshVariableSynchronizerList::MeshVariableSynchronizerList(BufAddrMng* bam) :
  m_buf_addr_mng (bam)
{
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
  m_vars.add(new CellMatVarScalSync<DataType>(var_menv, m_buf_addr_mng));
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

INST_MESH_VAR_SYNC_LIST_ADD(Integer);
INST_MESH_VAR_SYNC_LIST_ADD(Real);
INST_MESH_VAR_SYNC_LIST_ADD(Real3);
INST_MESH_VAR_SYNC_LIST_ADD(Real3x3);

