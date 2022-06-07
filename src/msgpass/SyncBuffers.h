#ifndef MSG_PASS_SYNC_BUFFERS_H
#define MSG_PASS_SYNC_BUFFERS_H

#include "accenv/AcceleratorUtils.h"
#include "msgpass/MeshVariableSynchronizerList.h"

#include <arcane/utils/MultiArray2.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Encapsule les buffers des valeurs à envoyer/recevoir                      */
/*---------------------------------------------------------------------------*/
class MultiBufView {
 public:
  MultiBufView();
  MultiBufView(ArrayView<Byte*> ptrs, Int64ConstArrayView sizes);
  MultiBufView(const MultiBufView& rhs);

  MultiBufView& operator=(const MultiBufView& rhs);

  //! Convertit un buffer d'octets en buffer de DataType
  template<typename DataType>
  static ArrayView<DataType> valBuf(ArrayView<Byte> buf);

  //! Convertit un buffer d'octets en buffer 2D de DataType dont la taille dans la deuxième dimension est dim2_size
  template<typename DataType>
  static Array2View<DataType> valBuf2(ArrayView<Byte> buf, Integer dim2_size);

  //! Accès en lecture/écriture au i-ème buffer d'octets
  ArrayView<Byte> byteBuf(Integer i);

  //! Retourne [beg_ptr, end_ptr[ qui contient tous les buffers (peut-être espacés de trous)
  Span<Byte> rangeSpan();
  ArrayView<Byte> rangeView();

 protected:
  SharedArray<Byte*> m_ptrs;
  SharedArray<Int64> m_sizes;
};

/*---------------------------------------------------------------------------*/
/* Encapsule les buffers des valeurs à envoyer/recevoir en 2 dimensions      */
/*---------------------------------------------------------------------------*/
class MultiBufView2 {
 public:
  MultiBufView2();
  MultiBufView2(ArrayView<Byte*> ptrs, Int64ConstArrayView sizes, 
      Integer dim1_sz, Integer dim2_sz);
  MultiBufView2(const MultiBufView2& rhs);

  MultiBufView2& operator=(const MultiBufView2& rhs);

  //! Accès à la sous-vue MultiBufView pour l'indice i1 dans la première dimension
  MultiBufView multiView(Integer i1);

  //! Retourne [beg_ptr, end_ptr[ qui contient tous les buffers (peut-être espacés de trous)
  Span<Byte> rangeSpan();

 protected:
  Integer m_dim1_sz;
  Integer m_dim2_sz;
  SharedArray<Byte*> m_ptrs;
  SharedArray<Int64> m_sizes;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
class SyncBuffers {
 public:
  SyncBuffers(bool is_acc_avl);
  virtual ~SyncBuffers();

  //
  void resetBuf();

  // Ajout à partir du nb d'items par voisin
  template<typename DataType>
  void addEstimatedMaxSz(IntegerConstArrayView item_sizes, Integer degree);

  // 
  void allocIfNeeded(Int64 buf_estim_sz);
  void allocIfNeeded();

  /*!
   * \brief A partir des nb d'items à communiquer, estime une borne sup de la taille du buffer en octets
   */
  template<typename DataType>
  static Int64 estimatedMaxBufSz(IntegerConstArrayView item_sizes, Integer degree);

  /*!
   * \brief A partir de la vue sur m_buf_bytes, construit une vue par voisin des buffers
   */
  template<typename DataType>
  MultiBufView multiBufView(
    IntegerConstArrayView item_sizes, Integer degree, Integer imem);

  /*!
   * \brief Construit des vues par voisin et par variable
   */
  MultiBufView2 multiBufViewVars(
    ConstArrayView<IMeshVarSync*> vars,
    IntegerConstArrayView item_sizes, Integer imem);

  /*!
   * \brief TODO : à commenter
   */
  MultiBufView2 multiBufViewVars(
    ConstArrayView<IMeshVarSync*> vars,
    Integer nb_nei, IMeshVarSync::eItemSync item_sync, Integer imem);

 protected:
  /*!
   * \brief A partir de la vue sur un buffer déjà alloué, construit une vue par voisin des buffers
   */
  template<typename DataType>
  MultiBufView _multiBufView(
      IntegerConstArrayView item_sizes, Integer degree,
      Span<Byte> buf_bytes);
  
  /*!
   * \brief TODO : à commenter
   */
  MultiBufView2 _multiBufViewVars(ConstArrayView<IMeshVarSync*> vars,
      IntegerConstArrayView item_sizes,
      Span<Byte> buf_bytes);

  /*!
   * \brief TODO : à commenter
   */
  MultiBufView2 _multiBufViewVars(ConstArrayView<IMeshVarSync*> vars,
    Integer nb_nei, IMeshVarSync::eItemSync item_sync,
    Span<Byte> buf_bytes);

 protected:
  struct BufMem {
    UniqueArray<Byte> *m_buf=nullptr;
    Int64 m_first_av_pos=0;  //! Première position disponible
    void reallocIfNeededOnHost(Int64 wanted_size, bool is_acc_avl);
    void reallocIfNeededOnDevice(Int64 wanted_size);
  };
 protected:
  bool m_is_accelerator_available=false;  //! Vrai si un GPU est disponible pour les calculs
  Int64 m_buf_estim_sz=0;  //! Taille qui va servir à allouer
  // Pour gérer les buffers sur l'hote et le device
  BufMem m_buf_mem[2];
};

#endif

