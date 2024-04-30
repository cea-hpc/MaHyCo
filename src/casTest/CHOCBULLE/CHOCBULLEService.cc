#include "CHOCBULLEService.h"

void CHOCBULLEService::initMatMono(Integer dim)  {
    
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0.;
  }
}

void CHOCBULLEService::initVarMono(Integer dim, double* densite_initiale, double* energie_initiale, double* pression_initiale, 
                                    double* temperature_initiale, Real3x3 vitesse_initiale)  {
    
  /* Pas d'implémentation monoMat */
}
void CHOCBULLEService::initMat(Integer dim)  {
    
  String fichier_string = options()->fichier_bulles();
  const char* fichier= fichier_string.localstr();
  lectureFichier(fichier);
  Integer nb_bulles_fichier(m_posr.size());
  pinfo() << " Nombre de bulles " << nb_bulles_fichier;
  Integer nb_bulles(options()->nombreBullesX 
        * options()->nombreBullesY * options()->nombreBullesZ + options()->nombreBullesAleatoires + nb_bulles_fichier);
  Real3* Xb = (Real3 *)malloc(sizeof(Real3) * nb_bulles);
  Real* rb = (Real *)malloc(sizeof(Real) * nb_bulles);
  Real3 Position_1_bulle(options()->positionPremiereBulle);
  Real deltax(options()->deltaxBulle);
  Real deltay(options()->deltayBulle);
  Real deltaz(options()->deltazBulle);
  Real rayon(options()->rayonBulle);

  
  Real XposFictif(options()->xfinFictif);
  Real XposFictif2(options()->xdebutFictif2);
  Real YminFictif2haut(options()->ydebutFictif2);
  Real YmaxFictif2bas(options()->yfinFictif2);
  
  Real3 PosCreneauMinMin(options()->positionCrenauMinMin);
  Real3 PosCreneauMinMax(options()->positionCrenauMinMax);
  Real3 PosCreneauMaxMin(options()->positionCrenauMaxMin);
  Real3 PosCreneauMaxMax(options()->positionCrenauMaxMax); 
  
  Real r(0.),rn(0.);
  Real Rmin(10.), Rmax(0.);
  pinfo() << " DEBUT INITIALISATION DES " << nb_bulles << " BULLES"; 
  
  pinfo() << " DEBUT INITIALISATION DES " << options()->nombreBullesAleatoires << " BULLES"; 
  srand( (unsigned)time( NULL ) );
  m_rand1 = (float *) malloc(sizeof(float) * options()->nombreBullesAleatoires);
  m_rand2 = (float *) malloc(sizeof(float) * options()->nombreBullesAleatoires);
  m_rand3 = (float *) malloc(sizeof(float) * options()->nombreBullesAleatoires);
  m_rand4 = (float *) malloc(sizeof(float) * options()->nombreBullesAleatoires);
  for (Integer ibulle=0; ibulle < options()->nombreBullesAleatoires; ibulle++) {
    m_rand1[ibulle] = (float) std::rand() / RAND_MAX;
    m_rand2[ibulle] = (float) std::rand() / RAND_MAX;
    m_rand3[ibulle] = (float) std::rand() / RAND_MAX;
    m_rand4[ibulle] = (float) std::rand() / RAND_MAX;
  }
  
  ENUMERATE_CELL(icell, allCells()) {
    Cell cell = *icell;
    // Plaque par défaut
    m_materiau[cell] = 1.;
// -----------------------------------------------------------------------------
// ---- Bulles en réseau
// -----------------------------------------------------------------------------
    // Boucle sur le nombre de bulle
    for (Integer ibullex=0; ibullex <options()->nombreBullesX ; ibullex++) {
    for (Integer ibulley=0; ibulley <options()->nombreBullesY ; ibulley++) {
    for (Integer ibullez=0; ibullez <options()->nombreBullesZ ; ibullez++) {
            
        Integer ibulle= ibulley + ibullex * options()->nombreBullesX;
         // centre et rayon des bulles
        Xb[ibulle].x = Position_1_bulle.x + ibullex*deltax;
        Xb[ibulle].y = Position_1_bulle.y + ibulley*deltay;
        Xb[ibulle].z = Position_1_bulle.z + ibullez*deltaz;
        // pinfo() << " Bulle " << ibulle << " Position " << Xb[ibulle];
        // {0.35, 0.25, 0.0};
        rb[ibulle] = rayon;
        // 0.05;
        const Real x = m_cell_coord[cell].x - Xb[ibulle].x;
        const Real y = m_cell_coord[cell].y - Xb[ibulle].y;
        const Real z = m_cell_coord[cell].z - Xb[ibulle].z;
        r = std::sqrt(x*x + y*y);
        // bulle surchargera l'aire
        if (r <= rb[ibulle]) {
            m_materiau[cell] = 2.;
            ENUMERATE_NODE(inode, cell.nodes()) {
                rn = std::sqrt(m_node_coord[inode].x * m_node_coord[inode].x +
                    m_node_coord[inode].y * m_node_coord[inode].y);
                Rmin = std::min(rn, Rmin);
                Rmax = std::max(rn, Rmax);
            }
            if (Rmin < rb[ibulle] && Rmax > rb[ibulle]) {
                // maille bulle + matiere
                m_materiau[cell] = 1.5;
            }
        }
        
    }
    }
    }

// -----------------------------------------------------------------------------
// ---- Bulles aléatoires
// -----------------------------------------------------------------------------    
    nb_bulles = options()->nombreBullesAleatoires;
    Real3 position_min(options()->positionMin);
    Real3 position_max(options()->positionMax);
    Real rayon_min(options()->rayonMin);
    Real rayon_max(options()->rayonMax);
    for (Integer ibulle=0; ibulle < nb_bulles ; ibulle++) {
        Xb[ibulle].x = position_min.x + (position_max.x-position_min.x) *  m_rand1[ibulle];
        Xb[ibulle].y = position_min.y + (position_max.y-position_min.y) *  m_rand2[ibulle];
        Xb[ibulle].z = position_min.z + (position_max.z-position_min.z) *  m_rand3[ibulle];
        rb[ibulle] =  rayon_min + (rayon_max-rayon_min) *  m_rand4[ibulle];
        const Real x = m_cell_coord[cell].x - Xb[ibulle].x;
        const Real y = m_cell_coord[cell].y - Xb[ibulle].y;
        const Real z = m_cell_coord[cell].z - Xb[ibulle].z;
        r = std::sqrt(x*x + y*y);
        // pinfo() << " Position de la bulle " << ibulle << " : " << Xb[ibulle];
        // bulle surchargera l'aire
        if (r <= rb[ibulle]) {
            m_materiau[cell] = 2.;
            ENUMERATE_NODE(inode, cell.nodes()) {
                rn = std::sqrt(m_node_coord[inode].x * m_node_coord[inode].x +
                    m_node_coord[inode].y * m_node_coord[inode].y);
                Rmin = std::min(rn, Rmin);
                Rmax = std::max(rn, Rmax);
            }
            if (Rmin < rb[ibulle] && Rmax > rb[ibulle]) {
                // maille bulle + matiere
                m_materiau[cell] = 1.5;
            }
        }
        
    }
// -----------------------------------------------------------------------------
// ---- Bulles lues dans un fichier
// -----------------------------------------------------------------------------  
    nb_bulles = nb_bulles_fichier;
    for (Integer ibulle=0; ibulle < nb_bulles ; ibulle++) {
        const Real x = m_cell_coord[cell].x - m_posx[ibulle];
        const Real y = m_cell_coord[cell].y - m_posy[ibulle];
        r = std::sqrt(x*x + y*y);
        // pinfo() << " Position de la bulle " << ibulle << " : " << m_posx[ibulle] << " , " << m_posy[ibulle] ;
        // bulle surchargera l'aire
        if (r <= m_posr[ibulle]) {
            m_materiau[cell] = 2.;
            ENUMERATE_NODE(inode, cell.nodes()) {
                rn = std::sqrt(m_node_coord[inode].x * m_node_coord[inode].x +
                    m_node_coord[inode].y * m_node_coord[inode].y);
                Rmin = std::min(rn, Rmin);
                Rmax = std::max(rn, Rmax);
            }
            // pinfo() << " Cellule dans une bulle " << cell.localId();
            if (Rmin < m_posr[ibulle] && Rmax > m_posr[ibulle]) {
                // maille bulle + matiere
                m_materiau[cell] = 1.5;
            }
        }
    }
// -----------------------------------------------------------------------------
// ---- Fin de la locasition des Bulles
// -----------------------------------------------------------------------------   

     // Fictif en aval de la plaque :
    if (m_cell_coord[cell].x > XposFictif2) 
      m_materiau[cell] = 3.;
    
    // Fictif au dessus de la plaque :
    if (m_cell_coord[cell].y > YminFictif2haut) 
      m_materiau[cell] = 3.;
    
     // Fictif en dessous de la plaque :
    if (m_cell_coord[cell].y < YmaxFictif2bas) 
      m_materiau[cell] = 3.;
    
    if (m_cell_coord[cell].x > PosCreneauMinMin.x && m_cell_coord[cell].x < PosCreneauMinMax.x
            && m_cell_coord[cell].y > PosCreneauMinMin.y && m_cell_coord[cell].y < PosCreneauMaxMax.y) {          
      m_materiau[cell] = 3.;
    }
    // Fictif en amont de la plaque :
    if (m_cell_coord[cell].x < XposFictif) 
      m_materiau[cell] = 0.;
       
  }
}

void CHOCBULLEService::initVar(Integer dim, double* densite_initiale, double* energie_initiale, double* pression_initiale, 
                                    double* temperature_initiale, Real3x3 vitesse_initiale)  {

    Integer nb_bulles_fichier(m_posr.size());
    Integer nb_bulles(options()->nombreBullesX * 
        options()->nombreBullesY * options()->nombreBullesZ + options()->nombreBullesAleatoires + nb_bulles_fichier);
    Real3* Xb = (Real3 *)malloc(sizeof(Real3) * nb_bulles);
    Real* rb = (Real *)malloc(sizeof(Real) * nb_bulles);
    Real3 Position_1_bulle(options()->positionPremiereBulle);
    Real deltax(options()->deltaxBulle);
    Real deltay(options()->deltayBulle);
    Real deltaz(options()->deltazBulle);
    Real rayon(options()->rayonBulle);
    
    pinfo() << " Valeur nombre aleatoire " << m_rand1 << " " << m_rand2 << " " << m_rand3 << " " << m_rand4;
  
    Real XposFictif(options()->xfinFictif);
    Real XposFictif2(options()->xdebutFictif2);
    Real YminFictif2haut(options()->ydebutFictif2);
    Real YmaxFictif2bas(options()->yfinFictif2);
    
    
    Real3 PosCreneauMinMin(options()->positionCrenauMinMin);
    Real3 PosCreneauMinMax(options()->positionCrenauMinMax);
    Real3 PosCreneauMaxMin(options()->positionCrenauMaxMin);
    Real3 PosCreneauMaxMax(options()->positionCrenauMaxMax); 
  
    Real r(0.),rn(0.);
    Real Rmin(10.), Rmax(0.);
    Real fracvol(0.);
    // mise à zero puis initialisation des fractions de masses et volumes
    m_mass_fraction.fill(0.0);
    m_fracvol.fill(0.0);
    CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
    ENUMERATE_CELL(icell,allCells()) {
        Cell cell = *icell;
        AllEnvCell all_env_cell = all_env_cell_converter[cell]; 
        // Materiau par défaut
        m_density[cell] = densite_initiale[1];
        m_fracvol[cell] = 1.;
        m_mass_fraction[cell] = 1.;
        m_pressure[cell] = pression_initiale[1];
        m_internal_energy[cell] = energie_initiale[1];
        m_temperature[cell] = temperature_initiale[1];
// -----------------------------------------------------------------------------
// ---- Bulles en réseau
// -----------------------------------------------------------------------------  
        for (Integer ibullex=0; ibullex <options()->nombreBullesX ; ibullex++) {
        for (Integer ibulley=0; ibulley <options()->nombreBullesY ; ibulley++) {
        for (Integer ibullez=0; ibullez <options()->nombreBullesZ ; ibullez++) {
            
            Integer ibulle= ibulley + ibullex * options()->nombreBullesX;
            // centre et rayon des bulles
            Xb[ibulle].x = Position_1_bulle.x + ibullex*deltax;
            Xb[ibulle].y = Position_1_bulle.y + ibulley*deltay;
            Xb[ibulle].z = Position_1_bulle.z + ibullez*deltaz;
            // {0.35, 0.25, 0.0};
            rb[ibulle] = rayon;          
            // bulle surchargera le materiau en présence
            const Real x = m_cell_coord[cell].x - Xb[ibulle].x;
            const Real y = m_cell_coord[cell].y - Xb[ibulle].y;
            const Real z = m_cell_coord[cell].z - Xb[ibulle].z;
            r = std::sqrt(x*x + y*y);
            if (r <= rb[ibulle]) {
                m_density[cell] = densite_initiale[2];
                m_pressure[cell] = pression_initiale[2];
                m_internal_energy[cell] = energie_initiale[2];
                m_temperature[cell] = temperature_initiale[2];
                m_fracvol[cell] = 1.;
                m_mass_fraction[cell] = 1.;
                ENUMERATE_NODE(inode, cell.nodes()) {
                    rn = std::sqrt(m_node_coord[inode].x * m_node_coord[inode].x +
                        m_node_coord[inode].y * m_node_coord[inode].y);
                    Rmin = std::min(rn, Rmin);
                    Rmax = std::max(rn, Rmax);
                }
                if (Rmin < rb[ibulle] && Rmax > rb[ibulle]) {
                    if (all_env_cell.nbEnvironment() ==1) { 
                        pinfo() << "Maille limite de la bulle doit etre mixte"; 
                        exit(1);
                    }
                    fracvol = (rb[ibulle] - Rmin) / (Rmax - Rmin); // fraction de la bulle
                    m_density[cell] = fracvol * densite_initiale[2] + (1. - fracvol) * densite_initiale[1];
                    m_pressure[cell] = fracvol * pression_initiale[2] + (1. - fracvol) * pression_initiale[1];
                    m_internal_energy[cell] = fracvol * energie_initiale[2] + (1. - fracvol) * energie_initiale[1];
                    m_temperature[cell] = fracvol * temperature_initiale[2] + (1. - fracvol) * temperature_initiale[1];
                    
                    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
                        EnvCell ev = *ienvcell;  
                        Integer index_env = ev.environmentId();  
                        m_density[ev] = densite_initiale[index_env] ;
                        m_pressure[ev] = pression_initiale[index_env] ;
                        m_internal_energy[ev] = energie_initiale[index_env];
                        m_temperature[ev] = temperature_initiale[index_env];
                    }
                }
            }
        }
        }
        }
// -----------------------------------------------------------------------------
// ---- Bulles aléatoires
// -----------------------------------------------------------------------------    
        nb_bulles = options()->nombreBullesAleatoires;
        Real3 position_min(options()->positionMin);
        Real3 position_max(options()->positionMax);
        Real rayon_min(options()->rayonMin);
        Real rayon_max(options()->rayonMax);
        for (Integer ibulle=0; ibulle < nb_bulles ; ibulle++) {
            Xb[ibulle].x = position_min.x + (position_max.x-position_min.x) *  m_rand1[ibulle];
            Xb[ibulle].y = position_min.y + (position_max.y-position_min.y) *  m_rand2[ibulle];
            Xb[ibulle].z = position_min.z + (position_max.z-position_min.z) *  m_rand3[ibulle];
            rb[ibulle] =  rayon_min + (rayon_max-rayon_min) *  m_rand4[ibulle];
            const Real x = m_cell_coord[cell].x - Xb[ibulle].x;
            const Real y = m_cell_coord[cell].y - Xb[ibulle].y;
            const Real z = m_cell_coord[cell].z - Xb[ibulle].z;
            r = std::sqrt(x*x + y*y);
            if (r <= rb[ibulle]) {
                m_density[cell] = densite_initiale[2];
                m_pressure[cell] = pression_initiale[2];
                m_internal_energy[cell] = energie_initiale[2];
                m_temperature[cell] = temperature_initiale[2];
                m_fracvol[cell] = 1.;
                m_mass_fraction[cell] = 1.;
                ENUMERATE_NODE(inode, cell.nodes()) {
                    rn = std::sqrt(m_node_coord[inode].x * m_node_coord[inode].x +
                        m_node_coord[inode].y * m_node_coord[inode].y);
                    Rmin = std::min(rn, Rmin);
                    Rmax = std::max(rn, Rmax);
                }
                if (Rmin < rb[ibulle] && Rmax > rb[ibulle]) {
                    if (all_env_cell.nbEnvironment() ==1) { 
                        pinfo() << "Maille limite de la bulle doit etre mixte"; 
                        exit(1);
                    }
                    fracvol = (rb[ibulle] - Rmin) / (Rmax - Rmin); // fraction de la bulle
                    m_density[cell] = fracvol * densite_initiale[2] + (1. - fracvol) * densite_initiale[1];
                    m_pressure[cell] = fracvol * pression_initiale[2] + (1. - fracvol) * pression_initiale[1];
                    m_internal_energy[cell] = fracvol * energie_initiale[2] + (1. - fracvol) * energie_initiale[1];
                    m_temperature[cell] = fracvol * temperature_initiale[2] + (1. - fracvol) * temperature_initiale[1];
                    
                    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
                        EnvCell ev = *ienvcell;  
                        Integer index_env = ev.environmentId();  
                        m_density[ev] = densite_initiale[index_env] ;
                        m_pressure[ev] = pression_initiale[index_env] ;
                        m_internal_energy[ev] = energie_initiale[index_env];
                        m_temperature[ev] = temperature_initiale[index_env];
                    }
                }
            }
        }
// -----------------------------------------------------------------------------
// ---- Bulles lues dans un fichier
// -----------------------------------------------------------------------------  
        nb_bulles = nb_bulles_fichier;
        for (Integer ibulle=0; ibulle < nb_bulles ; ibulle++) {
            // les positions Xb et les rayons rb ont été lus dans un fichier
            const Real x = m_cell_coord[cell].x - m_posx[ibulle];
            const Real y = m_cell_coord[cell].y - m_posy[ibulle];
            r = std::sqrt(x*x + y*y);
            if (r <= m_posr[ibulle]) {
                m_density[cell] = densite_initiale[2];
                m_pressure[cell] = pression_initiale[2];
                m_internal_energy[cell] = energie_initiale[2];
                m_temperature[cell] = temperature_initiale[2];
                m_fracvol[cell] = 1.;
                m_mass_fraction[cell] = 1.;
                ENUMERATE_NODE(inode, cell.nodes()) {
                    rn = std::sqrt(m_node_coord[inode].x * m_node_coord[inode].x +
                        m_node_coord[inode].y * m_node_coord[inode].y);
                    Rmin = std::min(rn, Rmin);
                    Rmax = std::max(rn, Rmax);
                }
                if (Rmin < m_posr[ibulle] && Rmax > m_posr[ibulle]) {
                    if (all_env_cell.nbEnvironment() ==1) { 
                        pinfo() << "Maille limite de la bulle doit etre mixte"; 
                        exit(1);
                    }
                    fracvol = (rb[ibulle] - Rmin) / (Rmax - Rmin); // fraction de la bulle
                    m_density[cell] = fracvol * densite_initiale[2] + (1. - fracvol) * densite_initiale[1];
                    m_pressure[cell] = fracvol * pression_initiale[2] + (1. - fracvol) * pression_initiale[1];
                    m_internal_energy[cell] = fracvol * energie_initiale[2] + (1. - fracvol) * energie_initiale[1];
                    m_temperature[cell] = fracvol * temperature_initiale[2] + (1. - fracvol) * temperature_initiale[1];
                    
                    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
                        EnvCell ev = *ienvcell;  
                        Integer index_env = ev.environmentId();  
                        m_density[ev] = densite_initiale[index_env] ;
                        m_pressure[ev] = pression_initiale[index_env] ;
                        m_internal_energy[ev] = energie_initiale[index_env];
                        m_temperature[ev] = temperature_initiale[index_env];
                    }
                }
            }
        }
// -----------------------------------------------------------------------------
// ---- Fin de la locasition des Bulles
// -----------------------------------------------------------------------------   
        if (m_cell_coord[cell].x > XposFictif2) {
            m_density[cell] = densite_initiale[3];
            m_fracvol[cell] = 1.;
            m_mass_fraction[cell] = 1.;
            m_pressure[cell] = pression_initiale[3];
            m_internal_energy[cell] = energie_initiale[3];            
        }
        
        // Fictif au dessus de la plaque :
        if (m_cell_coord[cell].y > YminFictif2haut) {
            m_density[cell] = densite_initiale[3];
            m_fracvol[cell] = 1.;
            m_mass_fraction[cell] = 1.;
            m_pressure[cell] = pression_initiale[3];
            m_internal_energy[cell] = energie_initiale[3];   
        }
        // Fictif en dessous de la plaque :
        if (m_cell_coord[cell].y < YmaxFictif2bas) {
            m_density[cell] = densite_initiale[3];
            m_fracvol[cell] = 1.;
            m_mass_fraction[cell] = 1.;
            m_pressure[cell] = pression_initiale[3];
            m_internal_energy[cell] = energie_initiale[3]; 
        }
        if (m_cell_coord[cell].x > PosCreneauMinMin.x && m_cell_coord[cell].x < PosCreneauMinMax.x
                && m_cell_coord[cell].y > PosCreneauMinMin.y && m_cell_coord[cell].y < PosCreneauMaxMax.y) {             
            m_density[cell] = densite_initiale[3];
            m_fracvol[cell] = 1.;
            m_mass_fraction[cell] = 1.;
            m_pressure[cell] = pression_initiale[3];
            m_internal_energy[cell] = energie_initiale[3]; 
        }
        if (m_cell_coord[cell].x < XposFictif) {
            m_density[cell] = densite_initiale[0];
            m_fracvol[cell] = 1.;
            m_mass_fraction[cell] = 1.;
            m_pressure[cell] = pression_initiale[0];
            m_internal_energy[cell] = energie_initiale[0];
        }
    }
    
    ENUMERATE_NODE(inode, allNodes()){
        m_velocity[inode] = vitesse_initiale[0];
    } 
}
/**
 *******************************************************************************
 * \file lectureFichier()
 * \brief 
 *******************************************************************************
 */
void CHOCBULLEService::lectureFichier(const std::string& nomFichier)
{
    std::ifstream fichier(nomFichier); // Ouverture du fichier en lecture
    pinfo() << "Ouverture de " << nomFichier;
    
    if (fichier)
    {
        
        int compteur = 0;
        Real x(-1),y(-1),r(-1);
        String ligne; // premiere ligne à ne pas lire 
        fichier >> ligne; 
        while (fichier >> x >> y >> r)
        {
            compteur++;
            pinfo() << x << " " << y << " " << r ;
            m_posx.push_back(x);
            m_posy.push_back(y);
            m_posr.push_back(r);
        }
        fichier.close(); // Fermeture du fichier
    }
    else
    {
        std::cout << "Erreur lors de l'ouverture du fichier." << std::endl;
    }
}
/*---------------------------------------------------------------------------*/


bool CHOCBULLEService::hasReverseOption() { return options()->reverseOption;}
Real CHOCBULLEService::getReverseParameter() { return options()->parametre;}
bool CHOCBULLEService::isInternalModel() { return true; }
void CHOCBULLEService::initUtilisateur(Real3 vitesse_initiale) { }
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_CHOCBULLE(CHOCBULLE, CHOCBULLEService);
