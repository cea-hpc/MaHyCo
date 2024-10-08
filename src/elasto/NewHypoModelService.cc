// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "NewHypoModelService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/* Mise à zero de tous les tenseurs ou vecteurs liés à l'elasto-plasticité   */
/*---------------------------------------------------------------------------*/

void NewHypoModelService::initElasto(IMeshEnvironment* env)
{
     Real3x3 Identity = Real3x3(Real3(1.0, 0.0, 0.0), Real3(0.0, 1.0, 0.0), Real3(0.0, 0.0, 1.0));
     m_velocity_gradient.fill(Real3x3::zero());
     m_spin_rate.fill(Real3::zero());
     m_deformation_rate.fill(Real3x3::zero());
     m_strain_tensor.fill(Real3x3::zero());
     m_strain_tensor_n.fill(Real3x3::zero());
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    m_tensorF[cell] = Identity;
  }
    
}

/*---------------------------------------------------------------------------*/
/* Calcul des gradients de vitesses à la cell */
/*---------------------------------------------------------------------------*/

void NewHypoModelService::ComputeVelocityGradient(Real delta_t)
{
  Real3x3 velocity_gradient;
  Real3x3 displacement_gradient;
  Real volume;
  Real3x3 Identity = Real3x3(Real3(1.0, 0.0, 0.0), Real3(0.0, 1.0, 0.0), Real3(0.0, 0.0, 1.0));
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    velocity_gradient = Real3x3::zero();
    displacement_gradient = Real3x3::zero();
    for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
      const Real3 vi = m_velocity[cell.node(inode)];
      const Real3 ui = m_displacement[cell.node(inode)];
      Real3 cqsndemi = .5*(m_cell_cqs_n[icell][inode]+m_cell_cqs[icell][inode]);
      velocity_gradient +=math::prodTens(vi, cqsndemi);
    }
    
    volume = .5*(m_cell_volume_n[cell]+m_cell_volume[cell]);
    m_velocity_gradient[cell] = velocity_gradient / volume;
    m_displacement_gradient[cell] =  Identity  + velocity_gradient / volume * delta_t;
  }
}
/*---------------------------------------------------------------------------*/
/* Calcul du taux de déformations et du vecteur de rotation                  */
/*---------------------------------------------------------------------------*/

void NewHypoModelService::ComputeDeformationAndRotation()
{
  
 Real KoneOverTwo = 0.5;
    
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    m_deformation_rate[cell].x.x = m_velocity_gradient[cell].x.x;
    m_deformation_rate[cell].y.y = m_velocity_gradient[cell].y.y;
    m_deformation_rate[cell].z.z = m_velocity_gradient[cell].z.z;
    
    m_deformation_rate[cell].x.y = KoneOverTwo * (m_velocity_gradient[cell].x.y + m_velocity_gradient[cell].y.x);
    m_deformation_rate[cell].x.z = KoneOverTwo * (m_velocity_gradient[cell].x.z + m_velocity_gradient[cell].z.x);
    m_deformation_rate[cell].y.z = KoneOverTwo * (m_velocity_gradient[cell].y.z + m_velocity_gradient[cell].z.y);
    
    m_deformation_rate[cell].y.x = m_deformation_rate[cell].x.y;
    m_deformation_rate[cell].z.x = m_deformation_rate[cell].x.z;
    m_deformation_rate[cell].z.y = m_deformation_rate[cell].y.z;
    
    m_spin_rate[cell].x = KoneOverTwo * (m_velocity_gradient[cell].y.z - m_velocity_gradient[cell].z.y);
    m_spin_rate[cell].y = KoneOverTwo * (m_velocity_gradient[cell].z.x - m_velocity_gradient[cell].x.z);
    m_spin_rate[cell].z = KoneOverTwo * (m_velocity_gradient[cell].x.y - m_velocity_gradient[cell].y.x);
  
    m_tensorF[cell] = math::matrixProduct(m_displacement_gradient[cell],m_tensorF[cell]);
    Real3x3  transpose_tensorF = math::transpose(m_tensorF[cell]);
    m_gauchy_green_tensor[cell] = math::matrixProduct(m_tensorF[cell], transpose_tensorF);
    
  }
}
 /*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NewHypoModelService::ComputeElasticity(IMeshEnvironment* env, Real delta_t, Integer dim)
{
 Integer order(options()->ordreRotation);
 Real mu = getElasticCst(env);
 Real trace;
 Real3x3 Identity = Real3x3(Real3(1.0, 0.0, 0.0), Real3(0.0, 1.0, 0.0), Real3(0.0, 0.0, 1.0));
 Real3x3 strain_tensor_point(Real3x3::zero());
 Real FiveoverThree = 5./3.; 
 Real KoneOverThree = 1./3.;
 Real TwoOverThree = 2./3.;
 ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Cell cell = ev.globalCell();
    // sauvegarde du tenseur de l'iteration précédente
    m_strain_tensor_n[ev] = m_strain_tensor[ev];
    
    Real J = m_density_0[ev] / m_density[ev] ;
    // Comparaison avec F
    Real Jprime = determinant(m_tensorF[cell]);
    // calcul du nouveau tenseur
    // pour retrouver les résultats hypo s = s + 2 * mu (D-trace(D)/3) 
    // on pose 
    // J = 1;
    // et 
    // m_gauchy_green_tensor[cell] = Identity;
    Real mu_over_J = mu / math::pow(J, FiveoverThree);
    Real3x3 DB = math::matrixProduct(m_deformation_rate[cell], m_gauchy_green_tensor[cell]);
    Real3x3 BD = math::matrixProduct(m_gauchy_green_tensor[cell], m_deformation_rate[cell]);
    Real traceD = m_deformation_rate[cell].x.x + m_deformation_rate[cell].y.y + m_deformation_rate[cell].z.z ;
    Real traceB = m_gauchy_green_tensor[cell].x.x + m_gauchy_green_tensor[cell].y.y + m_gauchy_green_tensor[cell].z.z;
    Real3x3 devB = m_gauchy_green_tensor[cell] - KoneOverThree * traceB * Identity;
    
    
    strain_tensor_point = mu_over_J * ( DB + BD - TwoOverThree * math::doubleContraction(m_gauchy_green_tensor[cell], m_deformation_rate[cell]) * Identity
                                                              - FiveoverThree * traceD * devB) * delta_t ;
    // simplifié : m_strain_tensor[ev] = m_strain_tensor_n[ev] + mu * ( DB + BD - TwoOverThree * math::doubleContraction(m_gauchy_green_tensor[cell], m_deformation_rate[cell]) * Identity) * delta_t ;

   
    /* rotation sur le résultat */
    if (dim == 2) {
        Real sxx = m_strain_tensor_n[ev].x.x;
        Real sxy = m_strain_tensor_n[ev].x.y;
        Real syy = m_strain_tensor_n[ev].y.y;
        Real syx = m_strain_tensor_n[ev].y.x;

        if (order ==1 ) {
          // utilisation poduit tensoriel ?
          strain_tensor_point.x.x += 2 * delta_t * sxy * m_spin_rate[cell].z;
          strain_tensor_point.y.y -= 2 * delta_t * sxy * m_spin_rate[cell].z;
          strain_tensor_point.x.y -= delta_t *(sxx - syy) * m_spin_rate[cell].z;
          // symetrie
          strain_tensor_point.y.x =  strain_tensor_point.x.y;
        } else {
          // prise en compte consistante de la rotation
          //Q  = (Id +0.5dt*rotz) /  (Id -0.5dt*rotz)
          // les dexu matrices commutent et rotz etant antisymétrique, matrice de rotation
          // en posant A = - 0.5dt*rotz :  **** Report de la correction ordre 1 TESTER OK
          // Q  = I + (2/det(I+A))*(A+A*A)
          //      ( c -s 0 )
          // Q =  ( s  c 0 ) avec c = (1-a*a)/(1+a*a) et s = (2*a) / (1+a*a)
          //      ( 0  0 1 )
          Real a = -0.5 * delta_t * m_spin_rate[cell].z;
          Real qxx = (1.-a*a) / (1.+a*a);
          Real qyy = qxx; 
          Real qxy = -2.*a / (1.+a*a);
          Real qyx = -qxy;
          Real blxx = sxx*qxx+ sxy*qxy;
          Real blxy = sxx*qyx+ sxy*qyy;
          Real blyx = syx*qxx+ syy*qxy;
          Real blyy = syx*qyx+ syy*qyy;
          strain_tensor_point.x.x = qxx*blxx+qxy*blyx;
          strain_tensor_point.y.y = qyx*blxy+qyy*blyy;
          strain_tensor_point.x.y = qxx*blxy+qxy*blyy;
          // symetrie
          strain_tensor_point.y.x = strain_tensor_point.x.y;
        }
    } else {
        Real sxx = m_strain_tensor_n[ev].x.x;
        Real sxy = m_strain_tensor_n[ev].x.y;
        Real syy = m_strain_tensor_n[ev].y.y;
        Real szz = m_strain_tensor_n[ev].z.z;
        Real syz = m_strain_tensor_n[ev].y.z;
        Real szx = m_strain_tensor_n[ev].z.x;
        // cela reste à trace nulle
        // verifier le signe --> correction faite comme en 2D
        strain_tensor_point.x.x -= 2 * delta_t * (szx * m_spin_rate[cell].y - sxy * m_spin_rate[cell].z);
        strain_tensor_point.y.y -= 2 * delta_t * (sxy * m_spin_rate[cell].z - syz * m_spin_rate[cell].x);
        strain_tensor_point.z.z -= 2 * delta_t * (syz * m_spin_rate[cell].x - szx * m_spin_rate[cell].y);
        strain_tensor_point.x.y -= delta_t * (sxx - syy) * m_spin_rate[cell].z 
        - delta_t * syz * m_spin_rate[cell].y
        - delta_t * szx * m_spin_rate[cell].x;
        strain_tensor_point.y.z -= delta_t * (syy - szz) * m_spin_rate[cell].x 
        - delta_t * szx * m_spin_rate[cell].z
        - delta_t * sxy * m_spin_rate[cell].y;
        strain_tensor_point.z.x -= delta_t * (szz - sxx) * m_spin_rate[cell].y 
        - delta_t * sxy * m_spin_rate[cell].x
        - delta_t * syz * m_spin_rate[cell].z;
        // symetrie
        strain_tensor_point.y.x = strain_tensor_point.x.y;
        strain_tensor_point.z.y = strain_tensor_point.y.z;
        strain_tensor_point.x.z = strain_tensor_point.z.x;
    }
    m_strain_tensor[ev] =  m_strain_tensor_n[ev] + strain_tensor_point;

  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NewHypoModelService::ComputePlasticity(IMeshEnvironment* env, Real delta_t, Integer dim)
{
 // yield strength
 Real mu = getElasticCst(env);
 Real yield_strength = getLimitElasticCst(env);
 ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Cell cell = ev.globalCell();
    // coeff retour radial
    Real coeff(1.);
    Real intensite_deviateur(0.);
    if (dim == 2) {
      // invariant ici 0.5 (sij::sij) que l'on compare à (y**2)/3
      // ce qui revient à comparer sij::sij à 2 (y**2)/3
        intensite_deviateur  = math::pow(m_strain_tensor[ev].x.x,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].y.y,2.);
        intensite_deviateur += m_strain_tensor[ev].x.x * m_strain_tensor[ev].y.y;
        intensite_deviateur += math::pow(m_strain_tensor[ev].x.y,2.);
        // ou plus simple
        // intensite_deviateur  = 0.5 * math::doubleContraction(m_strain_tensor[ev],m_strain_tensor[ev]);
        
    } else {
      // invariant ici 0.5 (sij::sij)
        intensite_deviateur  = math::pow(m_strain_tensor[ev].x.x,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].y.y,2.);
        intensite_deviateur += m_strain_tensor[ev].x.x * m_strain_tensor[ev].y.y;
        intensite_deviateur += math::pow(m_strain_tensor[ev].x.y,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].y.z,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].z.x,2.);
        // ou plus simple
        // intensite_deviateur  = 0.5 * math::doubleContraction(m_strain_tensor[ev],m_strain_tensor[ev]);
        
    }
    if ( intensite_deviateur > math::pow(yield_strength,2.)/3.) {
        coeff = yield_strength/math::sqrt(3.*intensite_deviateur); 
        // retour radial
        m_strain_tensor[ev] *= coeff;
        // vitesse de déformation plastique
        // mu non nul car EPP
        m_plastic_deformation_velocity[ev] = (1.-coeff)*yield_strength/(3.*mu*coeff*delta_t);
        // deformation plastique
        m_plastic_deformation[ev] += m_plastic_deformation_velocity[ev];
        
        intensite_deviateur = math::pow(yield_strength,2.)/3.;
    }
    // pour les sorties
    m_strain_tensor_xx[cell] = - m_strain_tensor[ev].x.x;
    m_strain_tensor_yy[cell] = - m_strain_tensor[ev].y.y;
    m_strain_tensor_xy[cell] = - m_strain_tensor[ev].x.y;
    m_strain_tensor_xz[cell] = - m_strain_tensor[ev].x.z;
    m_strain_tensor_yz[cell] = - m_strain_tensor[ev].y.z;
    m_strain_norm[cell] = math::sqrt(intensite_deviateur);
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NewHypoModelService::ComputeElastoEnergie(IMeshEnvironment* env, Real delta_t)
{
 Real KoneOverFour = 0.25;
 ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;  
    Cell cell = ev.globalCell();
    
    m_internal_energy[ev] += KoneOverFour * delta_t * (m_cell_volume[ev] + m_cell_volume_n[ev]) * 
        math::doubleContraction(m_strain_tensor[ev], m_deformation_rate[cell])
        / m_cell_mass[ev];
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real NewHypoModelService::getElasticCst(IMeshEnvironment* env) { return options()->elasticCst();}
Real NewHypoModelService::getLimitElasticCst(IMeshEnvironment* env) { return options()->limitElastic();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_NEWHYPOMODEL(NewHypoModel, NewHypoModelService);
