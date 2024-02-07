#ifndef UTILS_ATOMIC_H
#define UTILS_ATOMIC_H

#include <mutex>

#ifdef ARCCORE_DEVICE_CODE
#include <cuda/atomic>
#endif
using namespace Arcane;

/*!
 * \brief Classe utilitaire d'opérations valides en multithread et en CUDA
 */
class UtilsAtomic
{
 public:
/*---------------------------------------------------------------------------*/
/*-------------------------------------------------------s--------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur de la classe UtilsAtomic
 */ 
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE UtilsAtomic():
  val_mutex()
{
}
/*---------------------------------------------------------------------------*/
/*!
 * \brief Constructeur par copie
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE UtilsAtomic([[maybe_unused]] const UtilsAtomic& rhs):
  val_mutex()
{
}
/*---------------------------------------------------------------------------*/
/*!
 * \brief Destructeur de la classe UtilsAtomic

 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE ~UtilsAtomic(){};

//*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Calcul atomique du minimum (int) entre val et old et écriture dans old, compatible Multithread et CUDA
 *
 * @param old : adresse du minimum actuel
 * @param val : valeur à tester
 */
/*---------------------------------------------------------------------------*/

 ARCCORE_HOST_DEVICE void implem_atomicMin(int*old, int val) const{
#ifdef ARCCORE_DEVICE_CODE
   atomicMin(old,val);
#else
   val_mutex.lock();
   *old = (val < *old)?val:*old;
   val_mutex.unlock();
#endif
}
/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Calcul atomique du maximum (int) entre val et old et écriture dans old, compatible Multithread et CUDA
 *
 * @param old : adresse du maximum actuel
 * @param val : valeur à tester
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE void implem_atomicMax(int*old, int val) const{
#ifdef ARCCORE_DEVICE_CODE
   atomicMax(old,val);
#else
   val_mutex.lock();
   *old = (val > *old)?val:*old;
   val_mutex.unlock();
#endif
}
/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Calcul atomique du minimum (Real) entre val et old et écriture dans old, compatible Multithread et CUDA
 *
 * @param old : adresse du minimum actuel
 * @param val : valeur à tester
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE void implem_atomicMin(Real*old, Real val) const{
#ifdef ARCCORE_DEVICE_CODE
   unsigned long long oldval, newval, readback;
   oldval = __double_as_longlong(*old);
   if(val > __longlong_as_double(oldval))
  {
   return;
  }
   newval = __double_as_longlong(val);
   while((readback = atomicCAS((unsigned long long*)old,oldval,newval)) != oldval)
   {
      oldval = readback;

      newval = newval = __double_as_longlong(val < __longlong_as_double(oldval)?val:__longlong_as_double(oldval));
   }
#else
   val_mutex.lock();
   *old = (val < *old)?val:*old;
   val_mutex.unlock();
#endif
}
/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Calcul atomique du maximum (Real) entre val et old et écriture dans old, compatible Multithread et CUDA
 *
 * @param old : adresse du maximum actuel
 * @param val : valeur à tester
 */
/*---------------------------------------------------------------------------*/
ARCCORE_HOST_DEVICE void implem_atomicMax(Real*old, Real val) const{
#ifdef ARCCORE_DEVICE_CODE

  unsigned long long oldval, newval, readback;
  oldval = __double_as_longlong(*old);
  if(val < __longlong_as_double(oldval))
  {
   return;
  }
   newval = __double_as_longlong(val);
   while((readback = atomicCAS((unsigned long long*)old,oldval,newval)) != oldval)
   {
      oldval = readback;

      newval = newval = __double_as_longlong(val > __longlong_as_double(oldval)?val:__longlong_as_double(oldval));
   }
#else
   val_mutex.lock();
   *old = (val > *old)?val:*old;
   val_mutex.unlock();
#endif
}

 private:
 //! mutex utilisé uniquement lorsque les opérations atomiques ne sont pas disponibles
   mutable std::mutex val_mutex;
};


#endif
