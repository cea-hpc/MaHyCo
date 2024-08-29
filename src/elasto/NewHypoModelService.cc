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
  Real volume;
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    velocity_gradient = Real3x3::zero();
    displacement_gradient = Real3x3::zero();
    for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
      const Real3 vi = m_velocity[cell.node(inode)];
      const Real3 ui = m_displacement[cell.node(inode)];
      Real3 cqsndemi = .5*(m_cell_cqs_n[icell][inode]+m_cell_cqs[icell][inode]);
      velocity_gradient +=math::prodTens(vi, cqsndemi);
      displacement_gradient +=math::prodTens(ui, cqsndemi);
    }
    
    volume = .5*(m_cell_volume_n[cell]+m_cell_volume[cell]);
    // m_velocity_gradient[cell] = velocity_gradient / m_cell_volume[cell];
    // m_displacement_gradient[cell] = displacement_gradient / m_cell_volume[cell];
    m_velocity_gradient[cell] = velocity_gradient / volume;
    m_displacement_gradient[cell] = displacement_gradient / volume;
    pinfo() << " volumes (new et old) " << m_cell_volume[cell] << " " << m_cell_volume_n[cell];
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
    
    pinfo() << " m_velocity_gradient[cell].x.x " << m_velocity_gradient[cell].x.x;
    pinfo() << " m_velocity_gradient[cell].y.y " << m_velocity_gradient[cell].y.y;
    pinfo() << " m_velocity_gradient[cell].x.y " << m_velocity_gradient[cell].x.y;
    pinfo() << " m_velocity_gradient[cell].y.x " << m_velocity_gradient[cell].y.x;

    pinfo() << " m_deformation_rate[cell].x.x " << m_deformation_rate[cell].x.x;
    pinfo() << " m_deformation_rate[cell].y.y " << m_deformation_rate[cell].y.y;
    pinfo() << " m_deformation_rate[cell].x.y " << m_deformation_rate[cell].x.y;
    pinfo() << " m_deformation_rate[cell].y.x " << m_deformation_rate[cell].y.x;

    pinfo() << " m_spin_rate[cell].z " << m_spin_rate[cell].z;
    pinfo() << "  m_deformation_rate[cell].z.z " <<m_deformation_rate[cell].z.z;

    pinfo() << " Gradient de déplacement " ;
    pinfo() << " m_displacement_gradient[cell].x.x " << m_displacement_gradient[cell].x.x;
    pinfo() << " m_displacement_gradient[cell].y.y " << m_displacement_gradient[cell].y.y;
    pinfo() << " m_displacement_gradient[cell].x.y " << m_displacement_gradient[cell].x.y;
    pinfo() << " m_displacement_gradient[cell].y.x " << m_displacement_gradient[cell].y.x;
    
    Real3x3  tensorF = Real3x3(Real3(1.0, 0.0, 0.0), Real3(0.0, 1.0, 0.0), Real3(0.0, 0.0, 1.0));
    tensorF +=  m_displacement_gradient[cell];
    Real3x3  transpose_tensorF = math::transpose(tensorF);
    m_gauchy_green_tensor[cell] = math::matrixProduct(tensorF, transpose_tensorF);
    pinfo() << " Tensor F " ;
    pinfo() << "tensorF.x.x " << tensorF.x.x;
    pinfo() << "tensorF.y.y " << tensorF.y.y;
    pinfo() << "tensorF.x.y " << tensorF.x.y;
    pinfo() << "tensorF.y.x " << tensorF.y.x;
    pinfo() << " transpose_Tensor F " ;
    pinfo() << " transpose_tensorF.x.x " <<  transpose_tensorF.x.x;
    pinfo() << " transpose_tensorF.y.y " <<  transpose_tensorF.y.y;
    pinfo() << " transpose_tensorF.x.y " <<  transpose_tensorF.x.y;
    pinfo() << " transpose_tensorF.y.x " <<  transpose_tensorF.y.x;
    
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
    pinfo() << " traceD " << traceD;
    pinfo() << " Matrice B ";
    pinfo() << "m_gauchy_green_tensor[cell].x.x " << m_gauchy_green_tensor[cell].x.x;
    pinfo() << "m_gauchy_green_tensor[cell].y.y " << m_gauchy_green_tensor[cell].y.y;
    pinfo() << "m_gauchy_green_tensor[cell].x.y " << m_gauchy_green_tensor[cell].x.y;
    pinfo() << "m_gauchy_green_tensor[cell].y.x " << m_gauchy_green_tensor[cell].y.x;
    pinfo() << " traceB " << traceB;
    
    trace = (m_strain_tensor[ev].x.x + m_strain_tensor[ev].y.y + m_strain_tensor[ev].z.z);
    pinfo() << " avant trace ";
    pinfo() << " m_strain_tensor[ev].x.x " << m_strain_tensor[ev].x.x;
    pinfo() << " m_strain_tensor[ev].y.y " << m_strain_tensor[ev].y.y;
    pinfo() << " m_strain_tensor[ev].x.y " << m_strain_tensor[ev].x.y;
    pinfo() << " m_strain_tensor[ev].y.x " << m_strain_tensor[ev].y.x;
    pinfo() << " trace " << trace/3.;

    // doit-on enlever la trace ?
    pinfo() << " avant rotation ";
    pinfo() << " m_strain_tensor[ev].x.x " << m_strain_tensor[ev].x.x;
    pinfo() << " m_strain_tensor[ev].y.y " << m_strain_tensor[ev].y.y;
    pinfo() << " m_strain_tensor[ev].x.y " << m_strain_tensor[ev].x.y;
    pinfo() << " m_strain_tensor[ev].y.x " << m_strain_tensor[ev].y.x;
    pinfo() << " trace " << trace/3.;
    
    /* rotation sur le résultat */
    if (dim == 2) {
        Real sxx = m_strain_tensor[ev].x.x;
        Real sxy = m_strain_tensor[ev].x.y;
        Real syy = m_strain_tensor[ev].y.y;
        Real syx = m_strain_tensor[ev].y.x;

        if (order ==1 ) {
        
          // utilisation poduit tensoriel ?
          m_strain_tensor[ev].x.x += 2 * delta_t * sxy * m_spin_rate[cell].z;
          m_strain_tensor[ev].y.y -= 2 * delta_t * sxy * m_spin_rate[cell].z;
          m_strain_tensor[ev].x.y -= delta_t *(sxx - syy) * m_spin_rate[cell].z;
          // symetrie
          m_strain_tensor[ev].y.x = m_strain_tensor[ev].x.y;
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
          m_strain_tensor[ev].x.x = qxx*blxx+qxy*blyx;
          m_strain_tensor[ev].y.y = qyx*blxy+qyy*blyy;
          m_strain_tensor[ev].x.y = qxx*blxy+qxy*blyy;
          // symetrie
          m_strain_tensor[ev].y.x = m_strain_tensor[ev].x.y;
        }

        /*
        trace = (m_strain_tensor[ev].x.x + m_strain_tensor[ev].y.y + m_strain_tensor[ev].z.z);
        pinfo() << " apres rotation ";
        pinfo() << " Effet sur xx et yy " << 2 * delta_t * m_strain_tensor[ev].x.y * m_spin_rate[cell].z;
        pinfo() << " m_strain_tensor[ev].x.x " << m_strain_tensor[ev].x.x;
        pinfo() << " m_strain_tensor[ev].y.y " << m_strain_tensor[ev].y.y;
        pinfo() << " m_strain_tensor[ev].x.y " << m_strain_tensor[ev].x.y;
        pinfo() << " m_strain_tensor[ev].y.x " << m_strain_tensor[ev].y.x;
        pinfo() << " trace " << trace/3.;
        pinfo() << " -------------------------------------------";
        pinfo() << " Solution analytique ";
        Real A(0.),l(0.), val(0.), val0(0.);
        Real Epsilon(1.e-8);
        Real Temps = (m_global_time() - m_global_deltat() * (1 - Epsilon));
        Real h(1.);
        Real Sf(0.6),Uf(0.5);
        Real Sb=Sf/h;
        Real Ub=Uf/h; 
        if (m_global_time() >= 0. &&  Temps <= 1.) {
          val0 =  math::log(1.+Uf/h);
          pinfo() << " FIN ETAPE 1";
          pinfo() << " Tensor.x.x analytique " << -(2./3.)*val0;
          pinfo() << " Tensor.y.y analytique " << (4./3.)*val0;
          val = 0.;
          pinfo() << " Tensor.x.y analytique " << val;
        }  else if ( Temps > 1.  &&  Temps <= 2.) {   
          val0 = math::log(1.+Uf/h);
          A = mu*(1.+math::log(1.+Uf/h));
          val = A*(1.-cos(Sf/(h+Uf)));
          pinfo() << " FIN ETAPE 2";
          pinfo() << " Tensor.x.x analytique " << -(2./3.)*(val0)+val;
          pinfo() << " Tensor.y.y analytique " << (4./3.)*(val0)-val;
          val = A * sin(Sf/(h+Uf));
          pinfo() << " Tensor.x.y analytique " << val;
        } else if ( Temps > 2. &&  Temps <= 3.) {
          A = mu*(1+math::log(1.+Uf/h));
          val = A*(1-cos(Sf/(h+Uf)));
          pinfo() << " FIN ETAPE 3";
          pinfo() << " Tensor.x.x analytique " << val;
          pinfo() << " Tensor.y.y analytique " << -val;
          val = A * sin(Sf/(h+Uf));
          pinfo() << " Tensor.x.y analytique " << val;
        } else if ( Temps > 3.) {
          val = 1-cos(Sb) + (1+math::log(1+Ub))*((1-cos((Sb)/(1+Ub)))*cos(Sb) - sin((Sb)/(1+Ub))*sin(Sb));
          pinfo() << " FIN ETAPE 4";
          pinfo() << " Tensor.x.x analytique " << val;
          pinfo() << " Tensor.y.y analytique " << -val;
          val = -sin(Sb) + (1.+math::log(1.+ Ub))*(sin(Sb/(1.+Ub))*cos(Sb)+(1-cos(Sb/(1+Ub)))*sin(Sb));
          pinfo() << " Tensor.x.y analytique " << val;
        }
        */
    } else {
        Real sxx = m_strain_tensor[ev].x.x;
        Real sxy = m_strain_tensor[ev].x.y;
        Real syy = m_strain_tensor[ev].y.y;
        Real szz = m_strain_tensor[ev].z.z;
        Real syz = m_strain_tensor[ev].y.z;
        Real szx = m_strain_tensor[ev].z.x;
        // cela reste à trace nulle
        // verifier le signe --> correction faite comme en 2D
        m_strain_tensor[ev].x.x -= 2 * delta_t * (szx * m_spin_rate[cell].y - sxy * m_spin_rate[cell].z);
        m_strain_tensor[ev].y.y -= 2 * delta_t * (sxy * m_spin_rate[cell].z - syz * m_spin_rate[cell].x);
        m_strain_tensor[ev].z.z -= 2 * delta_t * (syz * m_spin_rate[cell].x - szx * m_spin_rate[cell].y);
        m_strain_tensor[ev].x.y -= delta_t * (sxx - syy) * m_spin_rate[cell].z 
        - delta_t * syz * m_spin_rate[cell].y
        - delta_t * szx * m_spin_rate[cell].x;
        m_strain_tensor[ev].y.z -= delta_t * (syy - szz) * m_spin_rate[cell].x 
        - delta_t * szx * m_spin_rate[cell].z
        - delta_t * sxy * m_spin_rate[cell].y;
        m_strain_tensor[ev].z.x -= delta_t * (szz - sxx) * m_spin_rate[cell].y 
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
