// Copyright 2000-2024 CEA (www.cea.fr) 
// See the top-level COPYRIGHT file for details. 
// SPDX-License-Identifier: Apache-2.0
#ifndef _EOS_STDPERFECTGASACC1_MALLOCATOR_H
#define _EOS_STDPERFECTGASACC1_MALLOCATOR_H

#ifdef MAHYCO_HAS_BOOST
#include <boost/core/noinit_adaptor.hpp>
#endif

#include <arcane/utils/PlatformUtils.h>
#include <arccore/collections/IMemoryAllocator.h>
#include <arcane/utils/IMemoryRessourceMng.h>
#include <arccore/base/FatalErrorException.h>
#include "accenv/ProfAcc.h"

namespace Stdperfectgasacc1 {

template <class T>
struct ArcaneMemAllocator {
  typedef T value_type;

  ArcaneMemAllocator() = default;

  template <class U>
    constexpr ArcaneMemAllocator(const ArcaneMemAllocator<U>&) noexcept {}

  /*[[nodiscard]]*/ T* allocate(std::size_t n) 
  {
    if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
      throw std::bad_array_new_length();

    Arccore::MemoryAllocationArgs param_alloc;
    param_alloc.setMemoryLocationHint(Arccore::eMemoryLocationHint::MainlyDevice);

    auto allocated_mem_info = getAllocator()->allocate(param_alloc, n * sizeof(T));

    if (auto p = static_cast<T*>(allocated_mem_info.baseAddress())) {
      report(p, allocated_mem_info.capacity(), true);
      return p;
    }

    throw std::bad_alloc();
  }

  void deallocate(T* p, std::size_t n) {
    report(p, n, false);
    getAllocator()->deallocate(Arccore::MemoryAllocationArgs(), Arccore::AllocatedMemoryInfo(p, n, n));
  }

 private:
  void report([[maybe_unused]] T* p, [[maybe_unused]] std::size_t n, bool alloc) const 
  {
    /*
    std::cout << (alloc ? "Alloc: " : "Dealloc: ") << sizeof(T) * n
      << " bytes at " << std::hex << std::showbase
      << reinterpret_cast<void*>(p) << std::dec << '\n';
      */
    const char* str_msg = (alloc ? "MemAlloc" : "MemFree");
    prof_acc_mark(str_msg);
  }

  Arcane::IMemoryAllocator* getAllocator() {
    Arcane::IMemoryAllocator* allocator = Arcane::platform::getDefaultDataAllocator();
    //Arcane::IMemoryAllocator* allocator = Arcane::platform::getDataMemoryRessourceMng()->getAllocator(Arcane::eMemoryRessource::Device);
    if (!allocator)
      throw Arcane::FatalErrorException(A_FUNCINFO, "getDefaultDataAllocator() inexistant");
    return allocator;
  }
};


template<class T, class U>
bool operator==(const ArcaneMemAllocator <T>&, const ArcaneMemAllocator <U>&) { return true; }

template<class T, class U>
bool operator!=(const ArcaneMemAllocator <T>&, const ArcaneMemAllocator <U>&) { return false; }

#ifdef MALLOCATOR_IS_STD_ALLOC
template<typename T> using Mallocator = std::allocator<T>;
#else

#ifdef MAHYCO_HAS_BOOST
//#warning "MAHYCO_HAS_BOOST : boost::noinit_adaptor"
template<typename T> using Mallocator = boost::noinit_adaptor< ArcaneMemAllocator<T> >;
#else
//#warning "MAHYCO_HAS_BOOST NON DEFINI"
template<typename T> using Mallocator = ArcaneMemAllocator<T>;
#endif

#endif

} // Stdperfectgasacc1

#endif

