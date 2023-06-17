// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "DefaultModelService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/* Mise à zero de tous les tenseurs ou vecteurs liés à l'elasto-plasticité   */
/*---------------------------------------------------------------------------*/

void DefaultModelService::initElasto(IMeshEnvironment* env)
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

void DefaultModelService::ComputeVelocityGradient()
{
  Real3x3 velocity_gradient;
  
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    velocity_gradient = Real3x3::zero();
    for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
      const Real3 vi = m_velocity[cell.node(inode)];
      velocity_gradient +=math::prodTens(vi, m_cell_cqs[icell] [inode]);
    }
    m_velocity_gradient[cell] = velocity_gradient / m_cell_volume[cell];
  }
}
/*---------------------------------------------------------------------------*/
/* Calcul du taux de déformations et du vecteur de rotation                  */
/*---------------------------------------------------------------------------*/

void DefaultModelService::ComputeDeformationAndRotation()
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
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DefaultModelService::ComputeElasticity(IMeshEnvironment* env, Real delta_t, Integer dim)
{
 Real mu = getElasticCst(env);
 Real trace;
 Real3x3 identity = Real3x3(Real3(1.0, 0.0, 0.0), Real3(0.0, 1.0, 0.0), Real3(0.0, 0.0, 1.0));
 ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    Cell cell = ev.globalCell();
    // sauvegarde du tenseur de l'iteration précédente
    m_strain_tensor_n[ev] = m_strain_tensor[ev];
    // calcul du nouveau tenseur
    m_strain_tensor[ev] += 2*mu*m_deformation_rate[cell] * delta_t;
    // trace
    trace = m_strain_tensor[ev].x.x + m_strain_tensor[ev].y.y + m_strain_tensor[ev].z.z;
    // on retire la trace sur les élément diagonaux
    m_strain_tensor[ev] -= trace * identity;
    // rotation
    if (dim == 2) {
        // utilisation poduit tensoriel ?
        m_strain_tensor[ev].x.x -= 2 * delta_t * m_strain_tensor_n[ev].x.y * m_spin_rate[cell].z;
        m_strain_tensor[ev].y.y += 2 * delta_t * m_strain_tensor_n[ev].x.y * m_spin_rate[cell].z;
        m_strain_tensor[ev].x.y += delta_t *(m_strain_tensor_n[ev].y.y - m_strain_tensor_n[ev].y.y) * m_spin_rate[cell].z;
        // symetrie
        m_strain_tensor[ev].y.x = m_strain_tensor[ev].x.y;
    } else {
        // cela reste à trace nulle 
        m_strain_tensor[ev].x.x += 2 * delta_t * (m_strain_tensor_n[ev].z.x * m_spin_rate[cell].y - m_strain_tensor_n[ev].x.y * m_spin_rate[cell].z);
        m_strain_tensor[ev].y.y += 2 * delta_t * (m_strain_tensor_n[ev].x.y * m_spin_rate[cell].z - m_strain_tensor_n[ev].y.z * m_spin_rate[cell].x);
        m_strain_tensor[ev].z.z += 2 * delta_t * (m_strain_tensor_n[ev].y.z * m_spin_rate[cell].x - m_strain_tensor_n[ev].z.x * m_spin_rate[cell].y);
        m_strain_tensor[ev].x.y += delta_t * (m_strain_tensor_n[ev].x.x - m_strain_tensor_n[ev].y.y ) * m_spin_rate[cell].z 
        - delta_t * m_strain_tensor_n[ev].y.z * m_spin_rate[cell].y
        - delta_t * m_strain_tensor_n[ev].z.x * m_spin_rate[cell].x;
        m_strain_tensor[ev].y.z += delta_t * (m_strain_tensor_n[ev].y.y - m_strain_tensor_n[ev].z.z ) * m_spin_rate[cell].x 
        - delta_t * m_strain_tensor_n[ev].z.x * m_spin_rate[cell].z
        - delta_t * m_strain_tensor_n[ev].x.y * m_spin_rate[cell].y;
        m_strain_tensor[ev].z.x += delta_t * (m_strain_tensor_n[ev].z.z - m_strain_tensor_n[ev].x.x ) * m_spin_rate[cell].y 
        - delta_t * m_strain_tensor_n[ev].x.y * m_spin_rate[cell].x
        - delta_t * m_strain_tensor_n[ev].y.z * m_spin_rate[cell].z;
        // symetrie
        m_strain_tensor[ev].y.x = m_strain_tensor[ev].x.y;
        m_strain_tensor[ev].z.y = m_strain_tensor[ev].y.z;
        m_strain_tensor[ev].x.z = m_strain_tensor[ev].z.x;
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DefaultModelService::ComputePlasticity(IMeshEnvironment* env, Real delta_t, Integer dim)
{
 // yield strength
 Real mu = getElasticCst(env);
 Real yield_strength = math::pow(getLimitElasticCst(env),2.)/3.;
 ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;   
    // coeff retour radial
    Real coeff(1.);
    Real intensite_deviateur(0.);
    if (dim == 2) {
        intensite_deviateur  = math::pow(m_strain_tensor[ev].x.x,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].y.y,2.);
        intensite_deviateur += math::pow(m_strain_tensor[ev].z.z,2.);
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
    if ( intensite_deviateur > yield_strength) coeff = math::sqrt(yield_strength / intensite_deviateur); 
    // retour radial
    m_strain_tensor[ev] *= coeff;
    // vitesse de déformation plastique
    // mu non nul car EPP
    m_plastic_deformation_velocity[ev] = (1.-coeff)/mu*math::sqrt(intensite_deviateur/3.)/delta_t;
    // deformation plastique
    m_plastic_deformation[ev] += m_plastic_deformation_velocity[ev];
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void DefaultModelService::ComputeElastoEnergie(IMeshEnvironment* env, Real delta_t)
{
 Real KoneOverFour = 0.25;
 ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;  
    Cell cell = ev.globalCell();
    m_internal_energy[ev] += KoneOverFour * delta_t * (m_cell_volume[ev] + m_cell_volume_n[ev]) *(
        (m_strain_tensor[ev].x.x + m_strain_tensor_n[ev].x.x) * m_deformation_rate[cell].x.x 
       +(m_strain_tensor[ev].x.y + m_strain_tensor_n[ev].x.y) * m_deformation_rate[cell].x.y 
       +(m_strain_tensor[ev].x.z + m_strain_tensor_n[ev].x.z) * m_deformation_rate[cell].x.z 
       +(m_strain_tensor[ev].y.x + m_strain_tensor_n[ev].y.x) * m_deformation_rate[cell].y.x 
       +(m_strain_tensor[ev].y.y + m_strain_tensor_n[ev].y.y) * m_deformation_rate[cell].y.y 
       +(m_strain_tensor[ev].y.z + m_strain_tensor_n[ev].y.z) * m_deformation_rate[cell].y.z 
       +(m_strain_tensor[ev].z.x + m_strain_tensor_n[ev].z.x) * m_deformation_rate[cell].z.x 
       +(m_strain_tensor[ev].z.y + m_strain_tensor_n[ev].z.y) * m_deformation_rate[cell].z.y 
       +(m_strain_tensor[ev].z.z + m_strain_tensor_n[ev].z.z) * m_deformation_rate[cell].z.z 
    )/ m_cell_mass[ev];
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real DefaultModelService::getElasticCst(IMeshEnvironment* env) { return options()->elasticCst();}
Real DefaultModelService::getLimitElasticCst(IMeshEnvironment* env) { return options()->limitElastic();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_DEFAULTMODEL(DefaultModel, DefaultModelService);
