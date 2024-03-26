#include "ADVECTIONService.h"


void ADVECTIONService::initMatMono(Integer dim)  {
    
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0;
  }
}
void ADVECTIONService::initMat(Integer dim)  {
    
  if (options()->casTest >= MonoAdvectionTx && options()->casTest <= MonoAdvectionRotation)  {
    initMatMono(dim);
    return;
  } 
  Real3 Xb;
  if (options()->casTest < AdvectionRotation) 
    Xb = {0.50, 0.50, 0.};
  else
    Xb = {0.50, 0.75, 0.};

  // rayon interne et externe
  double rb(0.25);
  
  
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    Real rmin(10.), rmax(0.);
    Real rminx(10.), rmaxx(0.);
    Real rminy(10.), rmaxy(0.);
    ENUMERATE_NODE(inode, cell.nodes()) {
      Real rnode = std::sqrt((m_node_coord[inode][0] - Xb[0]) *
                                (m_node_coord[inode][0] - Xb[0]) +
                                   (m_node_coord[inode][1] - Xb[1]) *
                                       (m_node_coord[inode][1] - Xb[1]));
      Real rnodex = m_node_coord[inode][0] - Xb[0];
      Real rnodey = m_node_coord[inode][1] - Xb[1];

      rmin = std::min(rmin, rnode);
      rmax = std::max(rmax, rnode);
      
      rminx = std::min(rminx, rnodex);
      rmaxx = std::max(rmaxx, rnodex);
      
      rminy = std::min(rminy, rnodey);
      rmaxy = std::max(rmaxy, rnodey);

    }
    
    // Air partout
    m_materiau[cell] = 0;
    // bulle surchargera l'aire
    // centre de la bulle
    double r = sqrt((m_cell_coord[cell][0] - Xb[0]) *
                        (m_cell_coord[cell][0] - Xb[0]) +
                    (m_cell_coord[cell][1] - Xb[1]) *
                        (m_cell_coord[cell][1] - Xb[1]));
    double rx = math::abs(m_cell_coord[cell][0] - Xb[0]);
    double ry = math::abs(m_cell_coord[cell][1] - Xb[1]);
    
    if ((rx < rb) && (ry < rb)) {
      // maille pure de carré
      m_materiau[cell] = 1;
    } 
  }
}
void ADVECTIONService::initVarMono(Integer dim, double* densite_initiale, double* energie_initiale, 
                                   double* pression_initiale, double* temperature_initiale, Real3x3 vitesse_initiale )  {
    
  Real3 Xb;
  if (options()->casTest < MonoAdvectionRotation) 
            Xb = {0.50, 0.50, 0.};
      else
            Xb = {0.50, 0.75, 0.};
  Real3 cc = {0.5, 0.5, 0.};
  // rayon interne et externe
  double rb(options()->rayonBulle);      
  // info() << " boucle sur les mailles";
  
  
  Real croissance_densite(1.);
  if (options()->densiteInterieureInitialeLineaire()) croissance_densite = 1./rb;
  
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    // parametres maille
    Real rmin(10.), rmax(0.);
    Real rminx(10.), rmaxx(0.);
    Real rminy(10.), rmaxy(0.);
    ENUMERATE_NODE(inode, cell.nodes()) {
      Real rnode = std::sqrt((m_node_coord[inode][0] - Xb[0]) *
                                       (m_node_coord[inode][0] - Xb[0]) +
                                   (m_node_coord[inode][1] - Xb[1]) *
                                       (m_node_coord[inode][1] - Xb[1]));
      Real rnodex = m_node_coord[inode][0] - Xb[0];
      Real rnodey = m_node_coord[inode][1] - Xb[1];
      
      rmin = std::min(rmin, rnode);
      rmax = std::max(rmax, rnode);
      
      rminx = std::min(rminx, rnodex);
      rmaxx = std::max(rmaxx, rnodex);
      
      rminy = std::min(rminy, rnodey);
      rmaxy = std::max(rmaxy, rnodey);
    }
    // Air partout
    m_density[cell] = densite_initiale[0];
    m_pressure[cell] = pression_initiale[0];
    // bulle surchargera l'aire
    // centre de la bulle
    double r = sqrt((m_cell_coord[cell][0] - Xb[0]) *
                        (m_cell_coord[cell][0] - Xb[0]) +
                    (m_cell_coord[cell][1] - Xb[1]) *
                        (m_cell_coord[cell][1] - Xb[1]));
    double rx = math::abs(m_cell_coord[cell][0] - Xb[0]);
    double ry = math::abs(m_cell_coord[cell][1] - Xb[1]);

    if ((rx < rb) && (ry < rb)) {
        // maille pure de carré
        m_density[cell] = densite_initiale[1];
        if (options()->densiteInterieureInitialeLineaire()) m_density[cell] = densite_initiale[1]*croissance_densite*r;
        m_pressure[cell] = pression_initiale[1];
    } 
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};    
    if ( options()->casTest == MonoAdvectionTx) m_velocity[inode].x = 1.;
    if ( options()->casTest == MonoAdvectionTy) m_velocity[inode].y = 1.;
    if ( options()->casTest == MonoAdvectionT45) m_velocity[inode] = {1.0, 1.0, 0.0};  
    if ( options()->casTest == MonoAdvectionRotation) {
      Real3 dd  = m_node_coord[inode] - cc;
      double theta = std::atan2(dd[1], dd[0]);
      double r = std::sqrt(dd[0] * dd[0] + dd[1] * dd[1]);
      double omega = 4. * Pi;
      m_velocity[inode].x = -r * omega * std::sin(omega * 0. + theta);
      m_velocity[inode].y = r * omega * std::cos(omega * 0. + theta);
    }
    // sauvegarde des valeurs initiales mises dans m_velocity_n
    m_velocity_n[inode] = m_velocity[inode];
  }
}
void ADVECTIONService::initVar(Integer dim, double* densite_initiale, double* energie_initiale, 
                               double* pression_initiale, double* temperature_initiale, Real3x3 vitesse_initiale)  {

  if (options()->casTest >= MonoAdvectionTx && options()->casTest <= MonoAdvectionRotation)  {
    initVarMono(dim, densite_initiale, energie_initiale, pression_initiale, temperature_initiale, vitesse_initiale);
    return;
  } 
  Real3 Xb;
  if (options()->casTest < AdvectionRotation) 
          Xb = {0.50, 0.50, 0.};
      else
          Xb = {0.50, 0.75, 0.};
  
  
  CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
  Real3 cc = {0.5, 0.5, 0.};
  // rayon interne et externe
  double rb(options()->rayonBulle);     
  Real croissance_densite(1.);
  if (options()->densiteInterieureInitialeLineaire()) croissance_densite = 1./rb;
  
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    // parametres maille
    Real rmin(10.), rmax(0.);
    Real rminx(10.), rmaxx(0.);
    Real rminy(10.), rmaxy(0.);

    ENUMERATE_NODE(inode, cell.nodes()) {
      Real rnode = std::sqrt((m_node_coord[inode][0] - Xb[0]) *
                                       (m_node_coord[inode][0] - Xb[0]) +
                                   (m_node_coord[inode][1] - Xb[1]) *
                                       (m_node_coord[inode][1] - Xb[1]));
      
      Real rnodex = m_node_coord[inode][0] - Xb[0];
      Real rnodey = m_node_coord[inode][1] - Xb[1];
      
      rmin = std::min(rmin, rnode);
      rmax = std::max(rmax, rnode);
      
      rminx = std::min(rminx, rnodex);
      rmaxx = std::max(rmaxx, rnodex);
      
      rminy = std::min(rminy, rnodey);
      rmaxy = std::max(rmaxy, rnodey);
    }
    // Air partout
    m_density[cell] = densite_initiale[0];
    m_pressure[cell] = pression_initiale[0];
    // bulle surchargera l'aire
    // centre de la bulle
    double r = sqrt((m_cell_coord[cell][0] - Xb[0]) *
                        (m_cell_coord[cell][0] - Xb[0]) +
                    (m_cell_coord[cell][1] - Xb[1]) *
                        (m_cell_coord[cell][1] - Xb[1]));
    double rx = math::abs(m_cell_coord[cell][0] - Xb[0]);
    double ry = math::abs(m_cell_coord[cell][1] - Xb[1]);
    if ((rx < rb) && (ry < rb)) {
      // maille pure de bulle
      m_density[cell] = densite_initiale[1];
      if (options()->densiteInterieureInitialeLineaire()) m_density[cell] = densite_initiale[1]*croissance_densite*r;
      m_pressure[cell] = pression_initiale[1];
    } /* else if (((rmaxx >= rb) && (rminx < rb)) || ((rmaxx >= rb) && (rminx < rb))) {
      // cas des cellules mailles mixtes
      double frac_b = (rb - rminx) / (rmaxx - rminx);
      AllEnvCell all_env_cell = all_env_cell_converter[cell];
      ENUMERATE_CELL_ENVCELL(ienvcell, all_env_cell) {
        EnvCell ev = *ienvcell;
        if (ev.environmentId() == 0) {
          m_density[ev] = 1.;
          m_fracvol[ev] = frac_b;
          m_mass_fraction[ev] = frac_b;
        }
        if (ev.environmentId() == 1) {
          m_density[ev] = 0.;
          m_fracvol[ev] = 1-frac_b;
          m_mass_fraction[ev] = 1-frac_b;
        }
        m_pressure[ev] = 0.;
      }
      m_density[cell] = frac_b;
      m_pressure[cell] = 0.;
    } */
  }
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};    
    if ( options()->casTest == AdvectionTx) m_velocity[inode].x = 1.;
    if ( options()->casTest == AdvectionTy) m_velocity[inode].y = 1.;
    if ( options()->casTest == AdvectionT45) m_velocity[inode] = {1.0, 1.0, 0.0};  
    if ( options()->casTest == AdvectionRotation) {
      Real3 dd  = m_node_coord[inode] - cc;
      double theta = std::atan2(dd[1], dd[0]);
      double r = std::sqrt(dd[0] * dd[0] + dd[1] * dd[1]);
      double omega = 4. * Pi;
      m_velocity[inode].x = -r * omega * std::sin(omega * 0. + theta);
      m_velocity[inode].y = r * omega * std::cos(omega * 0. + theta);
    }
    // sauvegarde des valeurs initiales mises dans m_velocity_n
    m_velocity_n[inode] = m_velocity[inode];
  }
}

/*---------------------------------------------------------------------------*/

bool ADVECTIONService::hasReverseOption() { return options()->reverseOption;}
Real ADVECTIONService::getReverseParameter() { return options()->parametre;}
bool ADVECTIONService::isInternalModel() { return true; }
void ADVECTIONService::initUtilisateur() { }

ARCANE_REGISTER_SERVICE_ADVECTION(ADVECTION, ADVECTIONService);
