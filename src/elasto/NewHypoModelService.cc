// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "NewHypoModelService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/* Mise à zero de tous les tenseurs ou vecteurs liés à l'elasto-plasticité   */
/*---------------------------------------------------------------------------*/

void NewHypoModelService::initElasto(IMeshEnvironment* env)
{
     m_velocity_gradient.fill(Real3x3::zero());
     m_spin_rate.fill(Real3::zero());
     m_deformation_rate.fill(Real3x3::zero());
     m_strain_tensor.fill(Real3x3::zero());
     m_strain_tensor_n.fill(Real3x3::zero());
    
}

/*---------------------------------------------------------------------------*/
/* Calcul des gradients de vitesses à la cell */
/*---------------------------------------------------------------------------*/

void NewHypoModelService::ComputeVelocityGradient()
{
  Real3x3 velocity_gradient;
  Real3x3 displacement_gradient;
  
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    velocity_gradient = Real3x3::zero();
    displacement_gradient = Real3x3::zero();
    for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
      const Real3 vi = m_velocity[cell.node(inode)];
      const Real3 ui = m_displacement[cell.node(inode)];
      velocity_gradient +=math::prodTens(vi, m_cell_cqs[icell] [inode]);
      displacement_gradient +=math::prodTens(ui, m_cell_cqs[icell] [inode]);
    }
    m_velocity_gradient[cell] = velocity_gradient / m_cell_volume[cell];
    m_displacement_gradient[cell] = displacement_gradient / m_cell_volume[cell];
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
    
    Real3x3  tensorF = Real3x3(Real3(1.0, 0.0, 0.0), Real3(0.0, 1.0, 0.0), Real3(0.0, 0.0, 1.0));
    tensorF +=  m_displacement_gradient[cell];
    Real3x3  transpose_tensorF = math::transpose(tensorF);
    m_gauchy_green_tensor[cell] = math::matrixProduct(tensorF, transpose_tensorF);
    
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NewHypoModelService::ComputeElasticity(IMeshEnvironment* env, Real delta_t, Integer dim)
{
 Real mu = getElasticCst(env);
 Real3x3 Identity = Real3x3(Real3(1.0, 0.0, 0.0), Real3(0.0, 1.0, 0.0), Real3(0.0, 0.0, 1.0));
 Real FiveoverThree = 5./3.; 
 Real KoneOverThree = 1./3.;
 Real TwoOverThree = 2./3.;
 ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Cell cell = ev.globalCell();
    // sauvegarde du tenseur de l'iteration précédente
    m_strain_tensor_n[ev] = m_strain_tensor[ev];
    
    Real J = m_density[ev] / m_density_0[ev] ;
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
    
    
    m_strain_tensor[ev] = m_strain_tensor_n[ev] + mu_over_J * ( DB + BD - TwoOverThree * math::doubleContraction(m_gauchy_green_tensor[cell], m_deformation_rate[cell]) * Identity
        - FiveoverThree * traceD * devB) * delta_t ;
    
    // rotation sur le résultat
    if (dim == 2) {
        Real sxx = m_strain_tensor[ev].x.x;
        Real sxy = m_strain_tensor[ev].x.y;
        Real syy = m_strain_tensor[ev].y.y;
        
        // utilisation poduit tensoriel ?
        m_strain_tensor[ev].x.x -= 2 * delta_t * sxy * m_spin_rate[cell].z;
        m_strain_tensor[ev].y.y += 2 * delta_t * sxy * m_spin_rate[cell].z;
        m_strain_tensor[ev].x.y += delta_t *(sxx - syy) * m_spin_rate[cell].z;
        // symetrie
        m_strain_tensor[ev].y.x = m_strain_tensor[ev].x.y;
    } else {
        Real sxx = m_strain_tensor[ev].x.x;
        Real sxy = m_strain_tensor[ev].x.y;
        Real syy = m_strain_tensor[ev].y.y;
        Real szz = m_strain_tensor[ev].z.z;
        Real syz = m_strain_tensor[ev].y.z;
        Real szx = m_strain_tensor[ev].z.x;
        // cela reste à trace nulle 
        m_strain_tensor[ev].x.x += 2 * delta_t * (szx * m_spin_rate[cell].y - sxy * m_spin_rate[cell].z);
        m_strain_tensor[ev].y.y += 2 * delta_t * (sxy * m_spin_rate[cell].z - syz * m_spin_rate[cell].x);
        m_strain_tensor[ev].z.z += 2 * delta_t * (syz * m_spin_rate[cell].x - szx * m_spin_rate[cell].y);
        m_strain_tensor[ev].x.y += delta_t * (sxx - syy) * m_spin_rate[cell].z 
        - delta_t * syz * m_spin_rate[cell].y
        - delta_t * szx * m_spin_rate[cell].x;
        m_strain_tensor[ev].y.z += delta_t * (syy - szz) * m_spin_rate[cell].x 
        - delta_t * szx * m_spin_rate[cell].z
        - delta_t * sxy * m_spin_rate[cell].y;
        m_strain_tensor[ev].z.x += delta_t * (szz - sxx) * m_spin_rate[cell].y 
        - delta_t * sxy * m_spin_rate[cell].x
        - delta_t * syz * m_spin_rate[cell].z;
        // symetrie
        m_strain_tensor[ev].y.x = m_strain_tensor[ev].x.y;
        m_strain_tensor[ev].z.y = m_strain_tensor[ev].y.z;
        m_strain_tensor[ev].x.z = m_strain_tensor[ev].z.x;
    }
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
        intensite_deviateur  = math::pow(m_strain_tensor[ev].x.x,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].y.y,2.);
        intensite_deviateur += m_strain_tensor[ev].x.x * m_strain_tensor[ev].y.y;
        intensite_deviateur += math::pow(m_strain_tensor[ev].x.y,2.);
    } else {
        intensite_deviateur  = math::pow(m_strain_tensor[ev].x.x,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].y.y,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].z.z,2.);
        intensite_deviateur +=  m_strain_tensor[ev].x.x * m_strain_tensor[ev].y.y;
        intensite_deviateur += math::pow(m_strain_tensor[ev].x.y,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].y.z,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].z.x,2.);
    }
    if ( intensite_deviateur > 2.*math::pow(yield_strength,2.)/3.) {
        coeff = yield_strength/math::sqrt(3.*intensite_deviateur/2.); 
        // retour radial
        m_strain_tensor[ev] *= coeff;
        // vitesse de déformation plastique
        // mu non nul car EPP
        m_plastic_deformation_velocity[ev] = (1.-coeff)*yield_strength/(3.*mu*coeff*delta_t);
        // deformation plastique
        m_plastic_deformation[ev] += m_plastic_deformation_velocity[ev];
        
        intensite_deviateur = 2.*math::pow(yield_strength,2.)/3.;
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
