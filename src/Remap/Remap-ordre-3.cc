// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "RemapADIService.h"

// fonctions pour l'ordre 3
// ----------------------------------
// fonction pour evaluer le gradient
Real RemapADIService::evaluate_grad(Real hm, Real h0, Real hp, Real ym,
                            Real y0, Real yp) {
  Real grad;
  grad = h0 / (hm + h0 + hp) *
         ((2. * hm + h0) / (h0 + hp) * (yp - y0) +
          (h0 + 2. * hp) / (hm + h0) * (y0 - ym));
  return grad;
}
// ----------------------------------
// fonction pour évaluer ystar
Real RemapADIService::evaluate_ystar(Real hmm, Real hm, Real hp, Real hpp,
                             Real ymm, Real ym, Real yp, Real ypp,
                             Real gradm, Real gradp) {
  Real ystar, tmp1, tmp2;
  tmp1 = (2. * hp * hm) / (hm + hp) *
         ((hmm + hm) / (2. * hm + hp) - (hpp + hp) / (2. * hp + hm)) *
         (yp - ym);
  tmp2 = -hm * (hmm + hm) / (2. * hm + hp) * gradp +
         hp * (hp + hpp) / (hm + 2. * hp) * gradm;
  ystar = ym + hm / (hm + hp) * (yp - ym) +
          1. / (hmm + hm + hp + hpp) * (tmp1 + tmp2);
  return ystar;
}
// ----------------------------------
// fonction pour évaluer fm
Real RemapADIService::evaluate_fm(Real x, Real dx, Real up, Real du,
                          Real u6) {
  Real fm;
  fm = up - 0.5 * x / dx * (du - (1. - 2. / 3. * x / dx) * u6);
  return fm;
}
// ----------------------------------
// fonction pour évaluer fr
Real RemapADIService::evaluate_fp(Real x, Real dx, Real um, Real du,
                          Real u6) {
  Real fp;
  fp = um + 0.5 * x / dx * (du - (1. - 2. / 3. * x / dx) * u6);
  return fp;
}
// ----------------------------------
// fonction pour initialiser la structure interval
Real2 RemapADIService::define_interval(Real a, Real b) {
  Real2 I;
  I[0] = math::min(a, b);
  I[1] = math::max(a, b);
  return I;
}
// ----------------------------------
// fonction pour calculer l'intersection entre deux intervals
Real2 RemapADIService::intersection(Real2 I1, Real2 I2) {
  Real2 I;
  if ((I1[1] < I2[0]) || (I2[1] < I1[0])) {
    I[0] = 0.;
    I[1] = 0.;
  } else {
    I[0] = math::max(I1[0], I2[0]);
    I[1] = math::min(I1[1], I2[1]);
  }
  return I;
}
// ----------------------------------
// fonction pour calculer le flux
Real RemapADIService::ComputeFluxOrdre3(Real ymmm, Real ymm, Real ym, Real yp,
                                Real ypp, Real yppp, Real hmmm,
                                Real hmm, Real hm, Real hp, Real hpp,
                                Real hppp, Real vdt) {

  Real flux;
  Real gradmm, gradm, gradp, gradpp;
  Real ystarm, ystar, ystarp;
  Real ym_m, ym_p, yp_m, yp_p;
  Real grad_m, grad_p, ym6, yp6;
  //
  gradmm = evaluate_grad(hmmm, hmm, hm, ymmm, ymm, ym);
  gradm = evaluate_grad(hmm, hm, hp, ymm, ym, yp);
  gradp = evaluate_grad(hm, hp, hpp, ym, yp, ypp);
  gradpp = evaluate_grad(hp, hpp, hppp, yp, ypp, yppp);
  //
  ystarm = evaluate_ystar(hmmm, hmm, hm, hp, ymmm, ymm, ym, yp, gradmm, gradm);
  ystar = evaluate_ystar(hmm, hm, hp, hpp, ymm, ym, yp, ypp, gradm, gradp);
  ystarp = evaluate_ystar(hm, hp, hpp, hppp, ym, yp, ypp, yppp, gradp, gradpp);
  //
  ym_m = ystarm;
  ym_p = ystar;
  yp_m = ystar;
  yp_p = ystarp;
  //
  grad_m = ym_p - ym_m;
  grad_p = yp_p - yp_m;
  //
  ym6 = 6. * (ym - 0.5 * (ym_m + ym_p));
  yp6 = 6. * (yp - 0.5 * (yp_m + yp_p));
  //
  if (vdt > 0.) {
    flux = evaluate_fm(vdt, hm, ym_p, grad_m, ym6);
  } else if (vdt < 0.) {
    flux = evaluate_fp(-vdt, hp, yp_m, grad_p, yp6);
  } else if (vdt == 0.) {
    return 0.;
  }
  // Limitation TVD
  Real num, nup, ym_ym, yp_ym;
  Real2 I1, I2, limiteur;
  num = vdt / hm;
  nup = vdt / hp;
  ym_ym = ym + (1. - num) / num * (ym - ymm);
  yp_ym = yp - (1. + nup) / nup * (yp - ypp);
  
  if (vdt >= 0.) {
    I1 = define_interval(ym, yp);
    I2 = define_interval(ym, ym_ym);
  } else {
    I1 = define_interval(ym, yp);
    I2 = define_interval(yp, yp_ym);
  }
  limiteur = intersection(I1, I2);
  if (flux < limiteur[0]) {
    flux = limiteur[0];
  }
  if (flux > limiteur[1]) {
    flux = limiteur[1];
  }
  //
  return flux;
}
