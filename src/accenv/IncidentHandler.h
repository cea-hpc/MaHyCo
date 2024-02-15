#ifndef ERROR_HANDLER_H
#define ERROR_HANDLER_H

#include "UtilsAtomic.h"

using namespace Arcane;
/*!
 * \brief Opérations disponibles pour le choix de la valeur à garder après incident
 */
enum class IncidentHandlerOperation {Min,Max};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*!
 * \brief Vue sur un IncidentHandler
 */
class IncidentHandlerView
{
  public:
  /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur de la classe IncidentHandlerView
 *@param flag adresse du booléen stockant l'état d'incident
 */ 
/*---------------------------------------------------------------------------*/
  IncidentHandlerView(bool* flag)
  {
    m_out_flag = flag;
  }
  
  /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur par copie de la classe IncidentHandlerView
 * @param hand Instance à copier
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE IncidentHandlerView(const IncidentHandlerView& hand):
    m_out_flag(hand.m_out_flag)
  {
  }
  /*---------------------------------------------------------------------------*/
/*!
 * \brief Destructeur de la classe IncidentHandlerView
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE ~IncidentHandlerView()
  {
  }
  /*---------------------------------------------------------------------------*/
/*!
 * \brief Déclare qu'un incident à eu lieu pour l'IncidentHandler associé
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE void declare() const
  {
    *m_out_flag = true;
  }


  private: 
  /*! pointeur sur le booléen du handler associé*/
  mutable bool*m_out_flag;

};
/*!
 * \brief Vue générique sur un ValuedIncidentHandler
 */
template <typename T, IncidentHandlerOperation OP> class ValuedIncidentHandlerView
{
  public:
  /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur de la classe ValuedIncidentHandlerView
 *@param stored adresse de la valeur associée à l'incident
 *@param utils adresse de l'objet pour les opérations atomiques
 */
/*------------- --------------------------------------------------------------*/
  ValuedIncidentHandlerView(int*stored, UtilsAtomic*utils)
  {
    m_out_stored_value = stored;
    m_utils_ptr = utils;
  }
  /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur par copie de la classe ValuedIncidentHandlerView
 * @param hand Instance à copier
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE ValuedIncidentHandlerView(const ValuedIncidentHandlerView& hand):
    m_out_stored_value(hand.m_out_stored_value),
    m_utils_ptr(hand.m_utils_ptr)
  {
  }
   /*---------------------------------------------------------------------------*/
/*!
 * \brief Destructeur de la classe ValuedIncidentHandlerView
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE ~ValuedIncidentHandlerView()
  {
  }
   /*---------------------------------------------------------------------------*/
/*!
 * \brief Déclare qu'un incident à eu lieu pour le ValuedIncidentHandler associé
 * @param item valeur à associer au handler
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE void declare(int item) const{
    ARCANE_ASSERT(false, ("ValuedIncidentHandlerView<T,OP>::declare(...) doit être spécialisé"));
  }


  private:
  /*!pointeur sur la valeur*/ 
  mutable T*m_out_stored_value;
  /*pointeur sur l'utilitaire pour les fonctions atomiques*/
  UtilsAtomic*m_utils_ptr;

};
/*!
 * \brief Spécialisation du ValuedIncidentHandlerView pour l'operation Max
 */
template <typename T> class ValuedIncidentHandlerView<T,IncidentHandlerOperation::Max> 
{
  public:
   /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur de la classe ValuedIncidentHandlerView
 *@param flag adresse du booléen stockant l'état d'incident
 */
/*---------------------------------------------------------------------------*/
  ValuedIncidentHandlerView(T*stored, UtilsAtomic*utils)
  {
    m_out_stored_value = stored;
    m_utils_ptr = utils;
  }
   /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur par copie de la classe ValuedIncidentHandlerView
 * @param hand Instance à copier
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE ValuedIncidentHandlerView(const ValuedIncidentHandlerView& hand):
    m_out_stored_value(hand.m_out_stored_value),
    m_utils_ptr(hand.m_utils_ptr)
  {
  }
   /*---------------------------------------------------------------------------*/
/*!
 * \brief Destructeur de la classe ValuedIncidentHandlerView
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE ~ValuedIncidentHandlerView()
  {
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Déclare qu'un incident à eu lieu pour le ValuedIncidentHandler associé
 * @param item valeur à associer au handler
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE void declare(T item) const
  {
    m_utils_ptr->implem_atomicMax(m_out_stored_value,item);
  }

  private:
  /*!pointeur sur la valeur*/ 
  mutable T*m_out_stored_value;
  /*pointeur sur l'utilitaire pour les fonctions atomiques*/
  UtilsAtomic*m_utils_ptr;
};
/*!
 * \brief Spécialisation du ValuedIncidentHandlerView pour l'operation Min
 */
template <typename T> class ValuedIncidentHandlerView<T,IncidentHandlerOperation::Min> 
{
  public:
   /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur de la classe ValuedIncidentHandlerView
 *@param flag adresse du booléen stockant l'état d'incident
 */
/*---------------------------------------------------------------------------*/
   ValuedIncidentHandlerView(T*stored, UtilsAtomic*utils)
  {
    m_out_stored_value = stored;
    m_utils_ptr = utils;
  }
 /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur par copie de la classe ValuedIncidentHandlerView
 * @param hand Instance à copier
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE ValuedIncidentHandlerView(const ValuedIncidentHandlerView& hand):
    m_out_stored_value(hand.m_out_stored_value),
    m_utils_ptr(hand.m_utils_ptr)
  {
  }
   /*---------------------------------------------------------------------------*/
/*!
 * \brief Destructeur de la classe ValuedIncidentHandlerView
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE ~ValuedIncidentHandlerView()
  {
  }
  /*---------------------------------------------------------------------------*/
/*!
 * \brief Déclare qu'un incident à eu lieu pour le ValuedIncidentHandler associé
 * @param item valeur à associer au handler
 */
/*---------------------------------------------------------------------------*/
  ARCCORE_HOST_DEVICE void declare(T item) const
  {
   m_utils_ptr->implem_atomicMin(m_out_stored_value,item);
  }
  private:
  /*!pointeur sur la valeur*/ 
  mutable T*m_out_stored_value;
  /*pointeur sur l'utilitaire pour les fonctions atomiques*/
  UtilsAtomic*m_utils_ptr;
};

/*!
 * \brief Classe pour gerer les incidents dans un contexte multithread et GPU sans stockage de valeur.
 */
class IncidentHandler
{
 public:
 /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur de la classe IncidentHandler
 * Initialise le flag à false
 */
/*---------------------------------------------------------------------------*/
  IncidentHandler()
      : m_flag(MemoryAllocationOptions(platform::getAcceleratorHostMemoryAllocator(), eMemoryLocationHint::MainlyHost), 1)
  {
    init();
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur par copie de la classe IncidentHandlerView
 * @param rhs Instance à copier
 */
/*---------------------------------------------------------------------------*/
  IncidentHandler(const IncidentHandler& rhs) : m_flag(rhs.m_flag) 
  {
  }
  /*---------------------------------------------------------------------------*/
/*!
 * \brief Destructeur de la classe IncidentHandler
 */
/*---------------------------------------------------------------------------*/
  ~IncidentHandler()
  {
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialise le handler associé 
 */
/*---------------------------------------------------------------------------*/
  void init()
  {
    m_flag[0] = false;
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief crée une vue associée au Incident Handler
 */
/*---------------------------------------------------------------------------*/
  IncidentHandlerView createView()
  {
    IncidentHandlerView v(&m_flag[0]);
    return v;
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialise le Incident Handler et crée une vue associée 
 */
/*---------------------------------------------------------------------------*/
  IncidentHandlerView initAndCreateView()
  {
    init();
    return createView();
  }

  /*---------------------------------------------------------------------------*/
/*!
 * \brief Teste si un incident à eu lieu
 @return renvoie true si un incident à eu lieu, false sinon
 */
/*---------------------------------------------------------------------------*/
  bool
  happened() const
  {
    return m_flag[0];
  }

 private:
 /*!flag de stockage de l'état d'incident*/
  UniqueArray<bool> m_flag;
};

/*!
 * \brief Classe générique pour gerer les incidents dans un contexte multithread et GPU avec stockage d'une valeur relative à l'incident
 */
template<typename T, IncidentHandlerOperation OP> 
class ValuedIncidentHandler
{
  public:
  /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur de la classe ValuedIncidentHandler
 */
/*---------------------------------------------------------------------------*/
  ValuedIncidentHandler():
  m_stored_value(MemoryAllocationOptions(platform::getAcceleratorHostMemoryAllocator(),eMemoryLocationHint::MainlyHost),1)
  {
    ARCANE_ASSERT(false, ("ValuedIncidentHandlerView<T,OP>::declare(...) doit être spécialisé"));
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur par copie de la classe ValuedIncidentHandler
 */
/*---------------------------------------------------------------------------*/
    ValuedIncidentHandler(const ValuedIncidentHandler& rhs):
  m_stored_value(rhs.m_stored_value),
    m_utils(rhs.m_utils)
  {
  }


/*---------------------------------------------------------------------------*/
/*!
 * \brief Destructeur de la classe ValuedIncidentHandler
 */
/*---------------------------------------------------------------------------*/
  ~ValuedIncidentHandler()
  {
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief crée une vue associée au Valued Incident Handler
 */
/*---------------------------------------------------------------------------*/
  ValuedIncidentHandlerView<T,OP> createView()
  {
     ValuedIncidentHandlerView<T,OP> v(&m_stored_value[0],&m_utils);
     return v;
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief accesseur de la valeur associée au handler
 */
/*---------------------------------------------------------------------------*/
  T value() const
  {
    return m_stored_value[0];
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Teste si un incident à eu lieu
 @return renvoie true si un incident à eu lieu, false sinon
 */
/*---------------------------------------------------------------------------*/
  bool happened() const
  {
    return m_stored_value[0] != std::numeric_limits<T>::max();
  }
 private:
 /*!variable stockant la valeur associée au handler*/
  UniqueArray<T> m_stored_value;
  /*pointeur sur l'utilitaire pour les fonctions atomiques*/
  UtilsAtomic m_utils;



};
 
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Spécialisation du ValuedIncidentHandler pour l'operation Min
 */
template<typename T>
class ValuedIncidentHandler<T,IncidentHandlerOperation::Min>
{
  public:
   /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur de la classe ValuedIncidentHandler
 */
/*---------------------------------------------------------------------------*/
  ValuedIncidentHandler():
    m_stored_value(MemoryAllocationOptions(platform::getAcceleratorHostMemoryAllocator(),eMemoryLocationHint::MainlyHost),1)
  {
    init();
  }
 /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur par copie de la classe ValuedIncidentHandler
 */
/*---------------------------------------------------------------------------*/
  ValuedIncidentHandler(const ValuedIncidentHandler& rhs):
  m_stored_value(rhs.m_stored_value),
    m_utils(rhs.m_utils)
  {
  }
 /*---------------------------------------------------------------------------*/
/*!
 * \brief Destructeur de la classe ValuedIncidentHandler
 */
/*---------------------------------------------------------------------------*/
  ~ValuedIncidentHandler() 
  {
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialise le handler associé 
 */
/*---------------------------------------------------------------------------*/
  void init()
  {
    m_stored_value[0] = std::numeric_limits<T>::max();
  }

  ValuedIncidentHandlerView<T,IncidentHandlerOperation::Min> createView()
  {
     ValuedIncidentHandlerView<T,IncidentHandlerOperation::Min> v(&m_stored_value[0],&m_utils);
     return v;
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialise le Incident Handler et crée une vue associée 
 */
/*---------------------------------------------------------------------------*/
  ValuedIncidentHandlerView<T,IncidentHandlerOperation::Min> initAndCreateView()
  {
    init();
    ValuedIncidentHandlerView<T,IncidentHandlerOperation::Min> v(&m_stored_value[0],&m_utils);
    return v;
  }

/*---------------------------------------------------------------------------*/
/*!
 * \brief accesseur de la valeur associée au handler
 */
/*---------------------------------------------------------------------------*/
  T value() const
  {
    return m_stored_value[0];
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Teste si un incident à eu lieu
 @return renvoie true si un incident à eu lieu, false sinon
 */
/*---------------------------------------------------------------------------*/
  bool happened() const
  {
    return m_stored_value[0] != std::numeric_limits<T>::max();
  }

 private:
 /*!variable stockant la valeur associée au handler*/
  UniqueArray<T> m_stored_value;
  /*pointeur sur l'utilitaire pour les fonctions atomiques*/
  UtilsAtomic m_utils;



};
/*!
 * \brief Spécialisation du ValuedIncidentHandler pour l'operation Max
 */
template<typename T>
class ValuedIncidentHandler<T,IncidentHandlerOperation::Max>
{
  public:
   /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur de la classe ValuedIncidentHandler
 */
/*---------------------------------------------------------------------------*/
  ValuedIncidentHandler():
    m_stored_value(MemoryAllocationOptions(platform::getAcceleratorHostMemoryAllocator(),eMemoryLocationHint::MainlyHost),1)
  {
    init();
  }
 /*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur par copie de la classe ValuedIncidentHandler
 */
/*---------------------------------------------------------------------------*/
  ValuedIncidentHandler(const ValuedIncidentHandler& rhs):
  m_stored_value(rhs.m_stored_value),
    m_utils(rhs.m_utils)
  {
  }
 /*---------------------------------------------------------------------------*/
/*!
 * \brief Destructeur de la classe ValuedIncidentHandler
 */
/*---------------------------------------------------------------------------*/
  ~ValuedIncidentHandler()
  {
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialise le handler associé 
 */
/*---------------------------------------------------------------------------*/
  void init()
  {
    m_stored_value[0] = std::numeric_limits<T>::lowest();
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief crée une vue associée au Valued Incident Handler
 */
/*---------------------------------------------------------------------------*/
  ValuedIncidentHandlerView<T,IncidentHandlerOperation::Max> createView()
  {
     ValuedIncidentHandlerView<T,IncidentHandlerOperation::Max> v(&m_stored_value[0],&m_utils);
     return v;
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Initialise le Incident Handler et crée une vue associée 
 */
/*---------------------------------------------------------------------------*/
  ValuedIncidentHandlerView<T,IncidentHandlerOperation::Max> initAndCreateView()
  {
    init();
    ValuedIncidentHandlerView<T,IncidentHandlerOperation::Max> v(&m_stored_value[0],&m_utils);
    return v;
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief accesseur de la valeur associée au handler
 */
/*---------------------------------------------------------------------------*/
 T value() const
  {
    return m_stored_value[0];
  }
/*---------------------------------------------------------------------------*/
/*!
 * \brief Teste si un incident à eu lieu
 @return renvoie true si un incident à eu lieu, false sinon
 */
/*---------------------------------------------------------------------------*/
  bool happened() const
  {
    return m_stored_value[0] != std::numeric_limits<T>::lowest();
  }
 private:
 /*!variable stockant la valeur associée au handler*/
  UniqueArray<T> m_stored_value;
  /*pointeur sur l'utilitaire pour les fonctions atomiques*/
  UtilsAtomic m_utils;

};

#endif
