#include "msgpass/PackTransfer.h"

/*---------------------------------------------------------------------------*/
/* async_transfer */
/*---------------------------------------------------------------------------*/
void async_transfer(Span<Byte> dst_buf, Span<const Byte> src_buf, RunQueue& queue) {

  ARCANE_ASSERT(src_buf.size()==dst_buf.size(), ("Les buffers src et dst n'ont pas la meme taille"));

  queue.copyMemory(ax::MemoryCopyArgs(dst_buf.data(), src_buf.data(), src_buf.size()).addAsync());
}

/*---------------------------------------------------------------------------*/
/* async_transfer */
/*---------------------------------------------------------------------------*/
void async_transfer(MultiBufView out_buf, MultiBufView in_buf, RunQueue& queue) {
  auto dst = out_buf.rangeSpan();
  auto src = in_buf.rangeSpan();
  queue.copyMemory(ax::MemoryCopyArgs(dst.data(), src.data(), src.size()).addAsync());
}

