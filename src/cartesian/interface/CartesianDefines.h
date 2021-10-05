#ifndef CARTESIAN_INTERFACE_CARTESIAN_DEFINES_H
#define CARTESIAN_INTERFACE_CARTESIAN_DEFINES_H

/*
 * Si -DIMPL_CART_CARTESIAN : utilisation de l'implementation cartésienne pure Cartesian:: (cartesian/)
 * Si -DIMPL_CART_ARCANE    : utilisation de l'implémentation Arcane:: dans arcane/cea/
 * Si rien n'est défini, Cartesian:: par défaut
 */
#if defined(IMPL_CART_CARTESIAN) && defined(IMPL_CART_ARCANE)
  #error "IMPL_CART_CARTESIAN et IMPL_CART_ARCANE ne peuvent pas etre definies en même temps"
#endif

#if defined(IMPL_CART_CARTESIAN)
//  #warning "USE_CARTESIAN_IMPL defini car IMPL_CART_CARTESIAN est defini"
  #ifndef USE_CARTESIAN_IMPL
  #define USE_CARTESIAN_IMPL
  #endif

#elif defined(IMPL_CART_ARCANE)
//  #warning "USE_CARTESIAN_IMPL n'est pas defini car IMPL_CART_ARCANE est defini"
  #undef USE_CARTESIAN_IMPL

#else
//  #warning "USE_CARTESIAN_IMPL defini par defaut"
  #ifndef USE_CARTESIAN_IMPL
  #define USE_CARTESIAN_IMPL
  #endif

#endif

#endif

