#include "msgpass/SyncBuffers.h"

#include <arcane/utils/IMemoryRessourceMng.h>
#include <arcane/utils/IndexOutOfRangeException.h>

/*---------------------------------------------------------------------------*/
/* MultiBufView                                                              */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Encapsule des vues sur plusieurs buffers de communication                 */
/*---------------------------------------------------------------------------*/
MultiBufView::MultiBufView()
{ }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MultiBufView::MultiBufView(ArrayView<Byte*> ptrs, Int64ConstArrayView sizes) :
  m_ptrs    (ptrs),
  m_sizes   (sizes)
{
  ARCANE_ASSERT(ptrs.size()==sizes.size(), ("ptrs.size()!=sizes.size()"));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MultiBufView::MultiBufView(const MultiBufView& rhs) :
  m_ptrs    (rhs.m_ptrs),
  m_sizes   (rhs.m_sizes)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MultiBufView& MultiBufView::operator=(const MultiBufView& rhs)
{
  m_ptrs    = rhs.m_ptrs;
  m_sizes   = rhs.m_sizes;
  return *this;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Convertit un buffer d'octets en buffer de DataType
template<typename DataType>
ArrayView<DataType> MultiBufView::valBuf(ArrayView<Byte> buf) {
  // buf.data() doit être aligné sur alignof(DataType)
  ARCANE_ASSERT(reinterpret_cast<size_t>(buf.data())%alignof(DataType)==0, 
      ("L'adresse buf.data() n'est pas aligne sur alignof(DataType)"));

  // buf.size() doit être un multiple de sizeof(DataType)
  ARCANE_ASSERT(buf.size()%sizeof(DataType)==0, 
      ("buf.size() n'est pas un multiple de sizeof(DataType)"));

  return ArrayView<DataType>(
      static_cast<Integer>(buf.size()/sizeof(DataType)), 
      reinterpret_cast<DataType*>(buf.data()));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Convertit un buffer d'octets en buffer 2D de DataType dont la taille dans la deuxième dimension est dim2_size
template<typename DataType>
Array2View<DataType> MultiBufView::valBuf2(ArrayView<Byte> buf, Integer dim2_size) {
  // buf.data() doit être aligné sur alignof(DataType)
  ARCANE_ASSERT(reinterpret_cast<size_t>(buf.data())%alignof(DataType)==0, 
      ("L'adresse buf.data() n'est pas aligne sur alignof(DataType)"));

  // buf.size() doit être un multiple de dim2_size*sizeof(DataType)
  ARCANE_ASSERT(buf.size()%(dim2_size*sizeof(DataType))==0, 
      ("buf.size() n'est pas un multiple de dim2_size*sizeof(DataType)"));

  return Array2View<DataType>(reinterpret_cast<DataType*>(buf.data()),
      static_cast<Integer>(buf.size()/(dim2_size*sizeof(DataType))), 
      dim2_size);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Accès en lecture/écriture au i-ème buffer d'octets
ArrayView<Byte> MultiBufView::byteBuf(Integer i) {
  return ArrayView<Byte>(m_sizes[i], m_ptrs[i]);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Retourne [beg_ptr, end_ptr[ qui contient tous les buffers (peut-être espacés de trous)
Span<Byte> MultiBufView::rangeSpan() {
  if (m_ptrs.size()==0) {
    return Span<Byte>();
  } else {
    Byte* beg_ptr=m_ptrs[0];
    Integer last = m_ptrs.size()-1;
    Byte* end_ptr=m_ptrs[last]+m_sizes[last];
    Int64 sz = end_ptr-beg_ptr;
    return Span<Byte>(beg_ptr, sz);
  }
}

ArrayView<Byte> MultiBufView::rangeView() {
  auto sp = rangeSpan();
  return ArrayView<Byte>(sp.size(), sp.data());
}

/*---------------------------------------------------------------------------*/
/* MultiBufView2                                                             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Encapsule des vues sur plusieurs buffers de communication en 2 dimensions */
/*---------------------------------------------------------------------------*/
MultiBufView2::MultiBufView2()
{ }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MultiBufView2::MultiBufView2(ArrayView<Byte*> ptrs, Int64ConstArrayView sizes,
    Integer dim1_sz, Integer dim2_sz) :
  m_dim1_sz (dim1_sz),
  m_dim2_sz (dim2_sz),
  m_ptrs    (ptrs),
  m_sizes   (sizes)
{
  ARCANE_ASSERT(ptrs.size()==sizes.size(), ("ptrs.size()!=sizes.size()"));
  ARCANE_ASSERT(ptrs.size()==(m_dim1_sz*m_dim2_sz), 
      ("ptrs.size()!=(dim1_sz*dim2_sz)"));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MultiBufView2::MultiBufView2(const MultiBufView2& rhs) :
  m_dim1_sz (rhs.m_dim1_sz),
  m_dim2_sz (rhs.m_dim2_sz),
  m_ptrs    (rhs.m_ptrs),
  m_sizes   (rhs.m_sizes)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
MultiBufView2& MultiBufView2::operator=(const MultiBufView2& rhs)
{
  m_dim1_sz = rhs.m_dim1_sz;
  m_dim2_sz = rhs.m_dim2_sz;
  m_ptrs    = rhs.m_ptrs;
  m_sizes   = rhs.m_sizes;
  return *this;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Accès en lecture/écriture au i-ème buffer d'octets
MultiBufView MultiBufView2::multiView(Integer i1) {
  if (i1<0 || i1>=m_dim1_sz) {
    throw IndexOutOfRangeException(A_FUNCINFO,
        String::format("Invalid dim1 index value={0}",i1),
        i1, 0, m_dim1_sz);
  }
  return MultiBufView(
      m_ptrs.subView(i1*m_dim2_sz, m_dim2_sz),
      m_sizes.subView(i1*m_dim2_sz, m_dim2_sz));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Retourne [beg_ptr, end_ptr[ qui contient tous les buffers (peut-être espacés de trous)
Span<Byte> MultiBufView2::rangeSpan() {
  if (m_ptrs.size()==0) {
    return Span<Byte>();
  } else {
    Byte* beg_ptr=m_ptrs[0];
    Integer last = m_ptrs.size()-1;
    Byte* end_ptr=m_ptrs[last]+m_sizes[last];
    Int64 sz = end_ptr-beg_ptr;
    return Span<Byte>(beg_ptr, sz);
  }
}

/*---------------------------------------------------------------------------*/
/* SyncBuffers                                                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
SyncBuffers::SyncBuffers(bool is_acc_avl) :
  m_is_accelerator_available (is_acc_avl)
{
  for(Integer imem(0) ; imem<2 ; ++imem) {
    m_buf_mem[imem].m_buf=nullptr;
    m_buf_mem[imem].m_first_av_pos=0;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
SyncBuffers::~SyncBuffers() {
  for(Integer imem(0) ; imem<2 ; ++imem) {
    delete m_buf_mem[imem].m_buf;
  }
}

/*---------------------------------------------------------------------------*/
/* A partir des items à communiquer, estime une borne sup de la taille du    */ 
/* buffer en octets                                                          */
/*---------------------------------------------------------------------------*/
template<typename DataType>
Int64 SyncBuffers::estimatedMaxBufSz(IntegerConstArrayView item_sizes, 
    Integer degree) {
  // HYPOTHESE 1 : même valeur de sizeof(DataType) sur CPU et GPU
  // HYPOTHESE 2 : même valeur de alignof(DataType) sur CPU et GPU
  // TODO : comment le vérifier ?
  Integer nb_nei = item_sizes.size(); // nb de voisins
  Int64 estim_max_buf_sz = 0;
  Int64 sizeof_item = sizeof(DataType)*degree;
  for(Integer inei=0 ; inei<nb_nei ; ++inei) {
    // Dans le pire des cas, le décalage est d'au plus alignof(DataType)-1 octets
    estim_max_buf_sz += (item_sizes[inei]*sizeof_item + alignof(DataType)-1);
  }
  return estim_max_buf_sz;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void SyncBuffers::resetBuf() {
  m_buf_estim_sz = 0;
  for(Integer imem(0) ; imem<2 ; ++imem) {
    m_buf_mem[imem].m_first_av_pos=0;
  }
}

/*---------------------------------------------------------------------------*/
/* Ajout à partir du nb d'items par voisin                                   */
/*---------------------------------------------------------------------------*/
template<typename DataType>
void SyncBuffers::addEstimatedMaxSz(IntegerConstArrayView item_sizes,
    Integer degree) {
  m_buf_estim_sz += estimatedMaxBufSz<DataType>(item_sizes, degree);
}

/*---------------------------------------------------------------------------*/
/* Reallocation dans la mémoire hôte */
/*---------------------------------------------------------------------------*/
void SyncBuffers::BufMem::reallocIfNeededOnHost(Int64 wanted_size, bool is_acc_avl) {
  if (m_buf==nullptr) {
    eMemoryRessource mem_res = (is_acc_avl ? eMemoryRessource::HostPinned : eMemoryRessource::Host);
    IMemoryAllocator* allocator = platform::getDataMemoryRessourceMng()->getAllocator(mem_res);
    m_buf = new UniqueArray<Byte>(allocator, wanted_size);
  }
  m_buf->resize(wanted_size);
}

/*---------------------------------------------------------------------------*/
/* Reallocation dans la mémoire device */
/*---------------------------------------------------------------------------*/
void SyncBuffers::BufMem::reallocIfNeededOnDevice(Int64 wanted_size) {
  if (m_buf==nullptr) {
    IMemoryAllocator* allocator = 
      platform::getDataMemoryRessourceMng()->getAllocator(eMemoryRessource::Device);
    m_buf = new UniqueArray<Byte>(allocator, wanted_size);
  }
  m_buf->resize(wanted_size);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void SyncBuffers::allocIfNeeded(Int64 buf_estim_sz) {
  m_buf_estim_sz = buf_estim_sz;
  // D'abord l'hote
  m_buf_mem[0].reallocIfNeededOnHost(buf_estim_sz, m_is_accelerator_available);

  // Puis le device si celui-ci existe
  if (m_is_accelerator_available) {
    m_buf_mem[1].reallocIfNeededOnDevice(buf_estim_sz);
  }
  if (!m_is_accelerator_available) {
    // Pour débugger, le buffer "device" se trouve dans la mémoire hôte
    m_buf_mem[1].reallocIfNeededOnHost(buf_estim_sz, m_is_accelerator_available);
  }
}

void SyncBuffers::allocIfNeeded() {
  allocIfNeeded(m_buf_estim_sz);
}

/*---------------------------------------------------------------------------*/
/* A partir du nb d'items par voisin item_sizes et d'un buffer de données déjà
 * alloué buf_bytes,
 * retourne une vue par voisin des buffers
 */
/*---------------------------------------------------------------------------*/
template<typename DataType>
MultiBufView SyncBuffers::_multiBufView(
    IntegerConstArrayView item_sizes, Integer degree,
    Span<Byte> buf_bytes) {

  if (estimatedMaxBufSz<DataType>(item_sizes, degree)>buf_bytes.size()) {
    // Il y a un risque que le buffer déjà alloué ne soit pas assez grand
    return MultiBufView();
  }

  Integer nb_nei = item_sizes.size(); // nb de voisins
  UniqueArray<Byte*> ptrs(nb_nei); // le pointeur de base du buffer par voisin
  Int64UniqueArray sizes_in_bytes(nb_nei); // la taille en octets du buffer par voisin

  Byte* cur_ptr{buf_bytes.data()};
  size_t available_space = buf_bytes.size();
  size_t sizeof_item = sizeof(DataType)*degree;
  Integer inei;

  for(inei=0 ; available_space>0 && inei<nb_nei ; ++inei) {
    // Par voisin, le tableau de valeurs doit être aligné sur alignof(DataType)
    void* cur_ptr_v = static_cast<void*>(cur_ptr);
    if (std::align(alignof(DataType), sizeof(DataType), cur_ptr_v, available_space)) {

      cur_ptr = static_cast<Byte*>(cur_ptr_v); // cur_ptr_v a été potentiellement modifié

      // Ici, cur_ptr a été modifié et est aligné sur alignof(DataType)
      // available_space a été diminué du nb d'octets = cur_ptr(après appel) - cur_ptr(avant appel)

      // Calcul en octets de l'occupation des valeurs pour le voisin inei
      size_t sz_nei_in_bytes = item_sizes[inei]*sizeof_item;

      ptrs[inei] = cur_ptr;
      sizes_in_bytes[inei] = sz_nei_in_bytes;

      cur_ptr += sz_nei_in_bytes; // ici, cur_ptr n'est plus forcement aligné avec alignof(T)
      if (sz_nei_in_bytes <= available_space) {
        available_space -= sz_nei_in_bytes;
      } else {
        ARCANE_ASSERT(false, ("Espace insuffisant pour aligner les données dans le buffer, available_space va devenir négatif"));
        break; // available_space ne pourra jamais être négatif car size_t est non signé
      }
    } else {
      ARCANE_ASSERT(false, ("Espace insuffisant pour aligner les données dans le buffer d'après std::align"));
      break;
    }
  }

  if (inei==nb_nei) {
    MultiBufView mb(ptrs, sizes_in_bytes);
    return mb;
  } else {
    // On ne devait jamais arriver là
    ARCANE_ASSERT(false, ("On ne devrait pas etre la"));
    return MultiBufView();
  }
}

/*---------------------------------------------------------------------------*/
/* */
/*---------------------------------------------------------------------------*/
template<typename DataType>
MultiBufView SyncBuffers::multiBufView(
    IntegerConstArrayView item_sizes, Integer degree, Integer imem) {

  auto& buf_mem = m_buf_mem[imem];
  Byte* new_ptr = buf_mem.m_buf->data()+buf_mem.m_first_av_pos;
  Int64 av_space = buf_mem.m_buf->size()-buf_mem.m_first_av_pos;
  Span<Byte> buf_bytes(new_ptr, av_space);

  auto mb = _multiBufView<DataType>(item_sizes, degree, buf_bytes);

  auto rg{mb.rangeSpan()}; // Encapsule [beg_ptr, end_ptr[
  Byte* end_ptr = rg.data()+rg.size();
  buf_mem.m_first_av_pos = (end_ptr - buf_mem.m_buf->data());
  return mb;
}

/*---------------------------------------------------------------------------*/
/* TODO
 */
/*---------------------------------------------------------------------------*/
MultiBufView2 SyncBuffers::_multiBufViewVars(ConstArrayView<IMeshVarSync*> vars,
    IntegerConstArrayView item_sizes,
    Span<Byte> buf_bytes) {

  Integer nb_nei = item_sizes.size(); // nb de voisins
  Integer nb_var = vars.size();  // nb de variables à synchroniser
  UniqueArray<Byte*> ptrs(nb_nei*nb_var); // le pointeur de base du buffer par voisin et par variable
  Int64UniqueArray sizes_in_bytes(nb_nei*nb_var); // la taille en octets du buffer par voisin et par variable

  Byte* cur_ptr{buf_bytes.data()};
  size_t available_space = buf_bytes.size();
  Integer inei;

  for(inei=0 ; inei<nb_nei ; ++inei) {
    for(Integer ivar=0 ; available_space>0 && ivar<nb_var ; ++ivar) {

      // Par voisin et par variable, le tableau de valeurs doit être aligné sur size_infos.alignOf;
      auto size_infos = vars[ivar]->sizeInfos();

      void* cur_ptr_v = static_cast<void*>(cur_ptr);
      if (std::align(size_infos.alignOf, size_infos.sizeOf, cur_ptr_v, available_space)) {

        cur_ptr = static_cast<Byte*>(cur_ptr_v); // cur_ptr_v a été potentiellement modifié

        // Ici, cur_ptr a été modifié et est aligné sur size_infos.alignOf
        // available_space a été diminué du nb d'octets = cur_ptr(après appel) - cur_ptr(avant appel)

        // Calcul en octets de l'occupation des valeurs pour le voisin inei
        size_t sz_nei_in_bytes = item_sizes[inei]*size_infos.sizeOfItem;

        ptrs[inei*nb_var+ivar] = cur_ptr;
        sizes_in_bytes[inei*nb_var+ivar] = sz_nei_in_bytes;

        cur_ptr += sz_nei_in_bytes; // ici, cur_ptr n'est plus forcement aligné avec alignof(T)
        if (sz_nei_in_bytes <= available_space) {
          available_space -= sz_nei_in_bytes;
        } else {
          throw NotSupportedException(A_FUNCINFO, 
              String("Espace insuffisant pour aligner les données dans le buffer, available_space va devenir négatif"));
          break; // available_space ne pourra jamais être négatif car size_t est non signé
        }
      } else {
        throw NotSupportedException(A_FUNCINFO, 
            String("Espace insuffisant pour aligner les données dans le buffer d'après std::align"));
        break;
      }
    }
  }

  if (inei==nb_nei) {
    MultiBufView2 mb2(ptrs, sizes_in_bytes, nb_nei, nb_var);
    return mb2;
  } else {
    // On ne devait jamais arriver là
    throw NotSupportedException(A_FUNCINFO, String("On ne devrait pas etre la"));
    return MultiBufView2();
  }
}

/*---------------------------------------------------------------------------*/
/* */
/*---------------------------------------------------------------------------*/
MultiBufView2 SyncBuffers::multiBufViewVars(
    ConstArrayView<IMeshVarSync*> vars,
    IntegerConstArrayView item_sizes, Integer imem) {

  auto& buf_mem = m_buf_mem[imem];
  Byte* new_ptr = buf_mem.m_buf->data()+buf_mem.m_first_av_pos;
  Int64 av_space = buf_mem.m_buf->size()-buf_mem.m_first_av_pos;
  Span<Byte> buf_bytes(new_ptr, av_space);

  auto mb2 = _multiBufViewVars(vars, item_sizes, buf_bytes);

  auto rg{mb2.rangeSpan()}; // Encapsule [beg_ptr, end_ptr[
  Byte* end_ptr = rg.data()+rg.size();
  buf_mem.m_first_av_pos = (end_ptr - buf_mem.m_buf->data());
  return mb2;
}


/*---------------------------------------------------------------------------*/
/* INSTANCIATIONS STATIQUES                                                  */
/*---------------------------------------------------------------------------*/
#include <arcane/utils/Real3x3.h>

#define INST_SYNC_BUFFERS(__DataType__) \
  template ArrayView<__DataType__> MultiBufView::valBuf<__DataType__>(ArrayView<Byte> buf); \
  template Array2View<__DataType__> MultiBufView::valBuf2<__DataType__>(ArrayView<Byte> buf, Integer dim2_size); \
  template MultiBufView SyncBuffers::multiBufView<__DataType__>(IntegerConstArrayView item_sizes, Integer degree, Integer imem); \
  template void SyncBuffers::addEstimatedMaxSz<__DataType__>(IntegerConstArrayView item_sizes, Integer degree); \
  template Int64 SyncBuffers::estimatedMaxBufSz<__DataType__>(IntegerConstArrayView item_sizes, Integer degree)

INST_SYNC_BUFFERS(Integer);
INST_SYNC_BUFFERS(Real);
INST_SYNC_BUFFERS(Real3);
INST_SYNC_BUFFERS(Real3x3);

