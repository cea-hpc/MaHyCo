#include "RIDERService.h"

/*
 * Initialisation 1 matériau
 */
void RIDERService::initMatMono(Integer dim)
{
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0;
  }
}

/*
 * Positionnement du matériau
 */
void RIDERService::initMat(Integer dim)
{

  if (options()->casTest >= MonoRiderTx && options()->casTest <= MonoRiderDeformationTimeReverse) {
    // Cas mono matériau
    initMatMono(dim);
    return;
  }

  Real3 Xb;
  if (options()->casTest < RiderRotation)
    Xb = {0.20, 0.20, 0.};
  else
    Xb = {0.50, 0.75, 0.};
  // rayon interne et externe
  double rb(options()->rayonBulle);

  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    Real rmin(10.), rmax(0.);
    ENUMERATE_NODE(inode, cell.nodes()) {
      Real rnode = std::sqrt((m_node_coord[inode].x - Xb.x) *
                                 (m_node_coord[inode].x - Xb.x) +
                             (m_node_coord[inode].y - Xb.y) *
                                 (m_node_coord[inode].y - Xb.y));
      rmin = std::min(rmin, rnode);
      rmax = std::max(rmax, rnode);
    }

    // Air partout
    m_materiau[cell] = 0;
    // bulle surchargera l'air

    if (rmax < rb) {
      // maille pure de bulle
      m_materiau[cell] = 1;
    }
    else if ((rmax >= rb) && (rmin < rb)) {
      double frac_b = (rb - rmin) / (rmax - rmin);
      m_materiau[cell] = frac_b;
    }
  }
}

/*
 * Init variables avec un seul matériau
 */
void RIDERService::initVarMono(  Integer dim, 
  SharedArray<double> densite_initiale, 
  SharedArray<double> energie_initiale, 
  SharedArray<double> pression_initiale, 
  SharedArray<double> temperature_initiale, 
  SharedArray<Real3> vitesse_initiale)
{
  Real3 Xb;
  if (options()->casTest < MonoRiderRotation)
    Xb = {0.20, 0.20, 0.};
  else
    Xb = {0.50, 0.75, 0.};

  Real3 cc = {0.5, 0.5, 0.};
  // rayon interne et externe
  double rb(options()->rayonBulle);
  Real r_cell(0.);
  Real croissance_densite(1.);
  if (options()->densiteInterieureInitialeLineaire())
    croissance_densite = 1. / rb;

  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    // parametres maille
    Real rmin(10.), rmax(0.);
    ENUMERATE_NODE(inode, cell.nodes()) {
      Real rnode = std::sqrt((m_node_coord[inode].x - Xb.x) *
                                 (m_node_coord[inode].x - Xb.x) +
                             (m_node_coord[inode].y - Xb.y) *
                                 (m_node_coord[inode].y - Xb.y));
      rmin = std::min(rmin, rnode);
      rmax = std::max(rmax, rnode);
    }
    // Air partout
    m_density[cell] = densite_initiale[0];
    m_pressure[cell] = pression_initiale[0];
    m_internal_energy[cell] = energie_initiale[0];
    // bulle surchargera l'air
    // centre de la bulle
    double r = sqrt((m_cell_coord[cell].x - Xb.x) *
                        (m_cell_coord[cell].x - Xb.x) +
                    (m_cell_coord[cell].y - Xb.y) *
                        (m_cell_coord[cell].y - Xb.y));
    if (options()->densiteInterieureInitialeLineaire())
      r_cell = r;
    else
      r_cell = rmax;
    if (r_cell < rb) {
      // maille pure de bulle
      m_density[cell] = densite_initiale[1];
      if (options()->densiteInterieureInitialeLineaire())
        m_density[cell] = densite_initiale[1] * std::min(croissance_densite * 2 * (rb - r), 1.);
      m_pressure[cell] = pression_initiale[1];
    }
    else if ((rmax >= rb) && (rmin < rb)) {
      if (!options()->densiteInterieureInitialeLineaire()) {
        double frac_b = (rb - rmin) / (rmax - rmin);
        m_density[cell] = frac_b * densite_initiale[1];
        m_pressure[cell] = pression_initiale[1];
      }
    }
  }

  // Init vitesse
  ENUMERATE_NODE(inode, allNodes()) {
    m_velocity[inode] = {0.0, 0.0, 0.0};
    if (options()->casTest == MonoRiderTx)
      m_velocity[inode].x = 1.;
    if (options()->casTest == MonoRiderTy)
      m_velocity[inode].y = 1.;
    if (options()->casTest == MonoRiderT45)
      m_velocity[inode] = {1.0, 1.0, 0.0};
    if (options()->casTest == MonoRiderRotation) {
      Real3 dd = m_node_coord[inode] - cc;
      double theta = std::atan2(dd[1], dd[0]);
      double r = std::sqrt(dd[0] * dd[0] + dd[1] * dd[1]);
      double omega = 4. * Pi;
      m_velocity[inode].x = -r * omega * std::sin(omega * 0. + theta);
      m_velocity[inode].y = r * omega * std::cos(omega * 0. + theta);
    }
    if (options()->casTest == MonoRiderVortex ||
        options()->casTest == MonoRiderVortexTimeReverse) {
      Real3 dd = m_node_coord[inode];
      m_velocity[inode].x =
          -2. * std::cos(Pi * dd[1]) * std::sin(Pi * dd[1]) *
          std::sin(Pi * dd[0]) * std::sin(Pi * dd[0]);
      m_velocity[inode].y =
          2. * std::cos(Pi * dd[0]) * std::sin(Pi * dd[0]) *
          std::sin(Pi * dd[1]) * std::sin(Pi * dd[1]);
    }
    if (options()->casTest == MonoRiderDeformation ||
        options()->casTest == MonoRiderDeformationTimeReverse) {
      Real3 dd = m_node_coord[inode] + cc;
      m_velocity[inode].x =
          std::sin(4. * Pi * dd[0]) * std::sin(4. * Pi * dd[1]);
      m_velocity[inode].y =
          std::cos(4. * Pi * dd[0]) * std::cos(4. * Pi * dd[1]);
    }
    // sauvegarde des valeurs initiales mises dans m_velocity_n
    m_velocity_n[inode] = m_velocity[inode];
  }
}

/*
 * Initialisation des variables
 */
void RIDERService::initVar(  Integer dim, 
  SharedArray<double> densite_initiale, 
  SharedArray<double> energie_initiale, 
  SharedArray<double> pression_initiale, 
  SharedArray<double> temperature_initiale, 
  SharedArray<Real3> vitesse_initiale)
{
  if (options()->casTest >= MonoRiderTx && options()->casTest <= MonoRiderDeformationTimeReverse) {
    initVarMono(dim, densite_initiale, energie_initiale, pression_initiale, temperature_initiale, vitesse_initiale);
    return;
  }
  Real3 Xb;
  if (options()->casTest < RiderRotation)
    Xb = {0.20, 0.20, 0.};
  else
    Xb = {0.50, 0.75, 0.};

  CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
  Real3 cc = {0.5, 0.5, 0.};
  // rayon interne et externe
  double rb(options()->rayonBulle);
  Real croissance_densite(1.);
  Real r_cell(0.);
  if (options()->densiteInterieureInitialeLineaire())
    croissance_densite = 1. / rb;

  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    // parametres maille
    Real rmin(10.), rmax(0.);
    ENUMERATE_NODE(inode, cell.nodes()) {
      Real rnode = std::sqrt((m_node_coord[inode].x - Xb.x) *
                                 (m_node_coord[inode].x - Xb.x) +
                             (m_node_coord[inode].y - Xb.y) *
                                 (m_node_coord[inode].y- Xb.y));
      rmin = std::min(rmin, rnode);
      rmax = std::max(rmax, rnode);
    }
    // Air partout
    m_density[cell] = densite_initiale[0];
    m_pressure[cell] = pression_initiale[0];
    m_internal_energy[cell] = energie_initiale[0];
  
    // bulle surchargera l'air
    // centre de la bulle
    double r = sqrt((m_cell_coord[cell].x - Xb.x) *
                        (m_cell_coord[cell].x - Xb.x) +
                    (m_cell_coord[cell].y- Xb.y) *
                        (m_cell_coord[cell].y - Xb.y));
    if (options()->densiteInterieureInitialeLineaire())
      r_cell = r;
    else
      r_cell = rmax;
    if (r_cell < rb) {
      // maille pure de bulle
      m_density[cell] = densite_initiale[1];
      if (options()->densiteInterieureInitialeLineaire())
        m_density[cell] = densite_initiale[1] * std::min(croissance_densite * 2 * (rb - r), 1.);
      m_pressure[cell] = pression_initiale[1];
    }
    else if ((rmax >= rb) && (rmin < rb)) {

      if (!options()->densiteInterieureInitialeLineaire()) {
        // cas des cellules mailles mixtes
        double frac_b = (rb - rmin) / (rmax - rmin);
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        ENUMERATE_CELL_ENVCELL(ienvcell, all_env_cell) {
          EnvCell ev = *ienvcell;
          if (ev.environmentId() == 0) {
            m_density[ev] = densite_initiale[1];
            m_fracvol[ev] = frac_b;
            m_mass_fraction[ev] = frac_b;
          }
          if (ev.environmentId() == 1) {
            m_density[ev] = densite_initiale[0];
            m_fracvol[ev] = 1. - frac_b;
            m_mass_fraction[ev] = 1. - frac_b;
          }
          m_pressure[ev] = 0.;
        }
        m_density[cell] = frac_b * densite_initiale[1];
        m_pressure[cell] = pression_initiale[1];
      }
    }
  }
  ENUMERATE_NODE(inode, allNodes()) {
    m_velocity[inode] = {0.0, 0.0, 0.0};
    if (options()->casTest == RiderTx)
      m_velocity[inode].x = 1.;
    if (options()->casTest == RiderTy)
      m_velocity[inode].y = 1.;
    if (options()->casTest == RiderT45)
      m_velocity[inode] = {1.0, 1.0, 0.0};
    if (options()->casTest == RiderRotation) {
      Real3 dd = m_node_coord[inode] - cc;
      double theta = std::atan2(dd[1], dd[0]);
      double r = std::sqrt(dd[0] * dd[0] + dd[1] * dd[1]);
      double omega = 4. * Pi;
      m_velocity[inode].x = -r * omega * std::sin(omega * 0. + theta);
      m_velocity[inode].y = r * omega * std::cos(omega * 0. + theta);
    }
    if (options()->casTest == RiderVortex ||
        options()->casTest == RiderVortexTimeReverse) {
      Real3 dd = m_node_coord[inode];
      m_velocity[inode].x =
          -2. * std::cos(Pi * dd[1]) * std::sin(Pi * dd[1]) *
          std::sin(Pi * dd[0]) * std::sin(Pi * dd[0]);
      m_velocity[inode].y =
          2. * std::cos(Pi * dd[0]) * std::sin(Pi * dd[0]) *
          std::sin(Pi * dd[1]) * std::sin(Pi * dd[1]);
    }
    if (options()->casTest == RiderDeformation ||
        options()->casTest == RiderDeformationTimeReverse) {
      Real3 dd = m_node_coord[inode] + cc;
      m_velocity[inode].x =
          std::sin(4. * Pi * dd[0]) * std::sin(4. * Pi * dd[1]);
      m_velocity[inode].y =
          std::cos(4. * Pi * dd[0]) * std::cos(4. * Pi * dd[1]);
    }
    // sauvegarde des valeurs initiales mises dans m_velocity_n
    m_velocity_n[inode] = m_velocity[inode];
  }
}

/*---------------------------------------------------------------------------*/

bool RIDERService::hasReverseOption() { return options()->reverseOption; }
Real RIDERService::getReverseParameter() { return options()->parametre; }
bool RIDERService::isInternalModel() { return true; }
void RIDERService::initUtilisateur(Real3 vitesse_initiale) {}

ARCANE_REGISTER_SERVICE_RIDER(RIDER, RIDERService);
