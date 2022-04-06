#ifndef ACC_ENV_BUF_ADDR_MNG_H
#define ACC_ENV_BUF_ADDR_MNG_H

#include "accenv/AcceleratorUtils.h"

#include <arcane/accelerator/core/RunQueueEvent.h>
#include <arcane/utils/IMemoryRessourceMng.h>
#include <arcane/materials/IMeshMaterialMng.h>

using namespace Arcane;
using namespace Arcane::Materials;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Manager of buffers of Int64
 */
class BufAddrMng {
 public:
  BufAddrMng(ax::Runner& runner, IMeshMaterialMng* mm) :
    m_mesh_mat_mng (mm)
  {
    bool is_acc_avl = AcceleratorUtils::isAvailable(runner);
    eMemoryRessource mem_h = (is_acc_avl ? eMemoryRessource::HostPinned : eMemoryRessource::Host);
    eMemoryRessource mem_d = (is_acc_avl ? eMemoryRessource::Device : eMemoryRessource::Host);

    IMemoryAllocator* alloc_h = platform::getDataMemoryRessourceMng()->getAllocator(mem_h);
    IMemoryAllocator* alloc_d = platform::getDataMemoryRessourceMng()->getAllocator(mem_d);

    m_nb_addr_per_buf = mm->environments().size()+1;
    Integer sz = 10*m_nb_addr_per_buf; 

    m_buf_addr_h = new UniqueArray<Int64>(alloc_h, sz);
    m_buf_addr_d = new UniqueArray<Int64>(alloc_d, sz);

    reset();

    // Event created once
    m_end_cpy_evt = ax::makeEventRef(runner);
  }

  virtual ~BufAddrMng() {
    delete m_buf_addr_h;
    delete m_buf_addr_d;
  }

  void reset() {
    m_cur_h=0;
    m_cur_d=0;
  }

  IMeshMaterialMng* materialMng() {
    return m_mesh_mat_mng;
  }

  //! Number of Int64 per buffer
  Integer nbAddr() const {
    return m_nb_addr_per_buf;
  }

  //! View for nbAddr() Int64 in HOST memory
  ArrayView<Int64> nextHostView() {
    return _nextView(m_buf_addr_h, m_cur_h);
  }

  //! View for nbAddr() Int64 in DEVICE memory
  ArrayView<Int64> nextDeviceView() {
    return _nextView(m_buf_addr_d, m_cur_d);
  }

  //! Asynchronous copy for all current buffers from Host to Device
  // Return an event which will occur after the copy
  Ref<ax::RunQueueEvent> asyncCpyHToD(RunQueue& queue) {
    if (m_cur_h != m_cur_d) {
      throw FatalErrorException(A_FUNCINFO,
          String::format(
            "Different buffer cursors between host m_cur_h={0} and device m_cur_d={1}",
            m_cur_h, m_cur_d));
    }
    Integer sz = m_cur_h*m_nb_addr_per_buf;

    // Asynchronous copy from Host to Device
    queue.copyMemory(ax::MemoryCopyArgs(
          m_buf_addr_d->subView(0, sz).data(), 
          m_buf_addr_h->subView(0, sz).data(), 
          sz*sizeof(Int64)).addAsync());

    queue.recordEvent(m_end_cpy_evt);

    reset(); // prepare the next stage

    return m_end_cpy_evt;
  }

 protected:

  ArrayView<Int64> _nextView(UniqueArray<Int64>* buf_addr, Integer& cur) {
    if (cur == (buf_addr->size()/m_nb_addr_per_buf)) {
      Integer sz = buf_addr->size() + 10*m_nb_addr_per_buf;
      // FIXME : resize plante si eMemoryRessource::Device, utiliser un
      // NumArray à la place ?
      buf_addr->resize(sz);
    }
    Integer n = m_nb_addr_per_buf;
    ArrayView<Int64> view = buf_addr->subView(cur*n, n);
    cur++;
    return view;
  }

 protected:
  IMeshMaterialMng* m_mesh_mat_mng=nullptr;
  Ref<ax::RunQueueEvent> m_end_cpy_evt;  //! Event occurs after copy
  UniqueArray<Int64>* m_buf_addr_h=nullptr;
  UniqueArray<Int64>* m_buf_addr_d=nullptr;
  Integer m_nb_addr_per_buf=0;
  Integer m_cur_h=0;
  Integer m_cur_d=0;
};

#endif

