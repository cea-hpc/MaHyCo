// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "TabuleeEOSService.h"


using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#define CONVERSION_DT 1.e3           // (de s à ms) 
#define CONVERSION_DENSITE 1.e3   // (de g/cm3 à kg/m3)
#define CONVERSION_PRESSION 1.e-1  // (de Erg/cm3 à Pa)
#define CONVERSION_ENERGIE  1.e-4  // (de J/Kg à Erg/g)
#define CONVERSION_TEMPERATURE 1.  // (de K à K)
#define CONVERSION_VITSON 1.e-6       // (de m2/s2 à m2.ms-2)

void TabuleeEOSService::initEOS(IMeshEnvironment* env)
{
    // rangement de la tabulation dans data
    initdata();
        
    std::string interp = "lin";
    double d0, t0, p0, e0;
    double epsilon ; 
    
    if (options()->initFormulation() == "VE") {
        
        info() << " Init avec densité et energie "; 
        ENUMERATE_ENVCELL(ienvcell,env) {
            EnvCell ev = *ienvcell;   
            Cell cell = ev.globalCell();
            
            m_density_0[ev] = m_density[ev];
            
            d0 = m_density[ev];  
            e0 = m_internal_energy[ev];
            
            calculPetT( d0,  e0,  t0,  p0,  interp );
            
            m_temperature[ev] = t0;
            m_pressure[ev] = p0;
            
            calculPetT( d0 + 1.e-7,  e0,  t0,  p0,  interp );
            
            if ((p0 -  m_pressure[ev]) > 0.) {
                m_sound_speed[ev] = math::sqrt((p0 - m_pressure[ev])/ 1.e-7);
            } else m_sound_speed[ev] = 1.e4;
            
            calculPetT( d0,  e0 + 1.e-7,  t0,  p0,  interp );
            
            m_dpde[ev] = (p0 -  m_pressure[ev])/1.e-7 ;
        }
    }
        
    if (options()->initFormulation() == "VT") {
        info() << " Init avec densité et temperature ";
        ENUMERATE_ENVCELL(ienvcell,env) {
            EnvCell ev = *ienvcell;   
            Cell cell = ev.globalCell();
            
            m_density_0[ev] = m_density[ev];
            
            d0 = m_density[ev] ;     
            t0 = m_temperature[ev] ;

            calculPetE( d0,  t0,  p0,  e0,  interp );
            
            m_pressure[ev] = p0 ;
            m_internal_energy[ev] = e0 ;
            
            calculPetE( d0 + 1.e-7,  t0,  p0,  e0,  interp );
            
            if ((p0 -  m_pressure[ev]) > 0.) {
                m_sound_speed[ev] = math::sqrt((p0 - m_pressure[ev])/ 1.e-7);
            } else m_sound_speed[ev] = 1.e4;
            
            m_dpde[ev] = 0.;
            
        }
    }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void TabuleeEOSService::ReinitEOS(IMeshEnvironment* env)
{
    // rangement de la tabulation dans data
    initdata();
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void TabuleeEOSService::applyEOS(IMeshEnvironment* env)
{
    std::string interp = "lin";
    double d0, t0, p0, e0;
    double epsilon ; 
    
    if (options()->formulation() == "VE") 
    // Init avec densité et energie 
    ENUMERATE_ENVCELL(ienvcell,env) {
        EnvCell ev = *ienvcell;   
        Cell cell = ev.globalCell();
        
        // std::cout << " cell " << cell.uniqueId() << " e = " << m_internal_energy[ev] << " emax " << emax << std::endl; 
        // pour éviter de sortir de la tabulation
        // d0 = std::max(std::min(m_density[ev], dmax/2), dmin);  
        // e0 = std::max(std::min(m_internal_energy[ev], emax/2), emin); 
        
       
        d0 = m_density[ev];  
        e0 = m_internal_energy[ev]; 
        
        calculPetT( d0,  e0,  t0,  p0,  interp );
        
        m_temperature[ev] = t0;
        m_pressure[ev] = p0;
        
        // pinfo() << " m_temperature[ev] " <<  t0 << " m_pressure[ev] " <<  p0;
        calculPetT( d0 + 1.e-7,  e0,  t0,  p0,  interp );
        
        if ((p0 -  m_pressure[ev]) > 0.) {
            m_sound_speed[ev] = math::sqrt((p0 - m_pressure[ev])/ 1.e-7);
        } else m_sound_speed[ev] = 1.e4;
        
        calculPetT( d0,  e0 + 1.e-7,  t0,  p0,  interp );
        
        m_dpde[ev] = (p0 -  m_pressure[ev])/1.e-7 ;
    }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/  
void TabuleeEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
        
  if (options()->formulation() == "VE") {
    std::string interp = "lin";
    double d0, t0, p0, e0;
    double epsilon ; 
    
    // pour éviter de sortir de la tabulation
    // d0 = std::max(std::min(m_density[ev], dmax/2), dmin);  
    // e0 = std::max(std::min(m_internal_energy[ev], emax/2), emin);
       
    d0 = m_density[ev];  
    e0 = m_internal_energy[ev]; 
    
    calculPetT( d0,  e0,  t0,  p0,  interp );
    
    m_temperature[ev] = t0;
    m_pressure[ev] = p0;
    
    calculPetT( d0 + 1.e-7,  e0,  t0,  p0,  interp );
    
    if ((p0 -  m_pressure[ev]) > 0.) {
        m_sound_speed[ev] = math::sqrt((p0 - m_pressure[ev])/ 1.e-7);
    } else m_sound_speed[ev] = 1.e4;
    
    calculPetT( d0,  e0 + 1.e-7,  t0,  p0,  interp );
    
    m_dpde[ev] = (p0 -  m_pressure[ev])/1.e-7 ;
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void TabuleeEOSService::Endommagement(IMeshEnvironment* env)
{
  Real damage_thresold = options()->tensionDamageThresold();
  Real density_thresold = options()->densityDamageThresold();
  Real density_max= options()->facteurDensityMaxDamage();
  ENUMERATE_ENVCELL(ienvcell,env)
  {
    EnvCell ev = *ienvcell;  
    Cell cell = ev.globalCell(); 
    // std::cout << " cell " << cell.uniqueId() << m_maille_endo[ev.globalCell()]  << std::endl; 
    if (m_maille_endo[ev.globalCell()] == 0) {
        // Maille saine : verification des seuils 
        // std::cout << " cell " << cell.uniqueId() << m_density[cell]/m_density_0[cell] << std::endl; 
        if (m_pressure[ev] < damage_thresold || m_density[cell]/m_density_0[cell] < density_thresold || m_density[cell]/m_density_0[cell] > density_max) {
            // maille devient endommagée
            m_maille_endo[ev] = 1;
            m_density_fracture[ev] = m_density[ev];
            m_internal_energy_fracture[ev] = m_internal_energy[ev];
            m_pressure[ev] = 0.;
            // std::cout << " cell endommagée " << cell.uniqueId()  << std::endl; 
        }
    } else { 
        // Maille endo : on verifie si elle ne s'est pas recompactée
        if (m_density[ev] > m_density_fracture[ev]) {
            // c'est le cas : on la déclare saine
            m_maille_endo[ev] = 0;
            m_internal_energy[ev] = m_internal_energy_fracture[ev];
        }
    }
  }
}
/*------------------------------------------------------------------*/
/* Initialisartion et rangement de la tabulation dans la map data   */
/*------------------------------------------------------------------*/
void TabuleeEOSService::initdata() {
    
    String  fichier_string = options()->fichierCoeff();
    std::string filename = fichier_string.localstr();
    if (!std::filesystem::exists(filename)) {
        std::cout << "Equation d'état : pas de fichier " << filename << " trouvée " << std::endl;
        data.clear();
        return;
    }
    std::ifstream file(filename);
    std::vector<std::string> lines;
    std::string line;
    
    while (std::getline(file, line)) {
            lines.push_back(line);
    }
    file.close();

    int j0 = 0;
    for (; j0 < lines.size(); ++j0) {
        if (!lines[j0].empty() && lines[j0][0] != '#') {
            break;
        }
    }
    lines.erase(lines.begin(), lines.begin() + j0);
    std::vector<std::string> tab = splitString(lines[0]);
    int npt = std::stoi(tab[2]);
    int npd = std::stoi(tab[4]);
    
    
	std::cout << " ici c'est ok npt =" << npt << " et npd="<< npd << std::endl;

    std::vector<std::vector<double>>  _temp( 1, std::vector<double>(npt, 0.0) );
    std::vector<std::vector<double>>  _dens( 1, std::vector<double>(npd, 0.0) );
    std::vector<std::vector<double>> energietot(npd, std::vector<double>(npt, 0.0));
    std::vector<std::vector<double>> pressiontot(npd, std::vector<double>(npt, 0.0));

	auto& temp = _temp[0];
	auto& dens = _dens[0];
    
    
    try {
        temp = splitAndConvertToDoubleForT(lines[1], ')');
    } catch (const std::invalid_argument& e) {
        std::cout << "Problem reading equation of state (temperature)" << std::endl;
        data.clear();
        exit(1);
        // return;
    }
	std::cout << "Lecture de valeurs jusqu'à " << npd <<  std::endl;
    j0 = 3;
    for (int jd = npd - 1; jd >= 0; --jd) {
        try {
        tab = splitString( lines.at(j0) );
            dens.at(jd) = std::stod(tab.at(2));
            pressiontot.at(jd) = splitAndConvertToDouble(lines.at(j0), ' ');
        } catch (const std::invalid_argument& e) {
            std::cout << "Problem reading equation of state (density and pressure)" << std::endl;
            data.clear();
            exit(1);
        }
        ++j0;
    }
	std::cout << " Passage à l'energie " << std::endl;
    j0 += 1;
    for (int jd = npd - 1; jd >= 0; --jd) {
        try {
        tab = splitString(lines.at(j0));
        energietot.at(jd) = splitAndConvertToDouble(lines.at(j0), ' ');
        } catch (const std::invalid_argument& e) {
            std::cout << "Problem reading equation of state (energy)" << std::endl;
            data.clear();
            exit(1);
        }
        ++j0;
    }
    std::cout << " Conversion et détermination des min et max en énergie et densité" << std::endl;
    emin = FloatInfo < Real >::maxValue();
    emax = 0.;
    dmin = FloatInfo < Real >::maxValue();
    dmax = 0.;
    tmin = FloatInfo < Real >::maxValue();
    tmax = 0.;

    for (int i = 0; i < npd; ++i) {
        for (int j = 0; j < npt; ++j) {
            // on met la pression en Pa (initialement en barye)
            pressiontot[i][j]*= CONVERSION_PRESSION;
            // on met l'energie en J/kg (initialement en erg/g)
            energietot[i][j]*= CONVERSION_ENERGIE;
            emin = math::min(energietot[i][j], emin);
            emax = math::max(energietot[i][j], emax);
            
            tmin = math::min(temp[j], tmin);
            tmax = math::max(temp[j], tmax);
            
        }
        // on met la densite en kg/m3 (initialement en g/cm3)
        dens[i]*= CONVERSION_DENSITE;
        dmin = math::min(dens[i], dmin);
        dmax = math::max(dens[i], dmax);
    }
	std::cout << " Rangement dans data " << std::endl;
        
    data = { {"t", _temp},
        //std::pair<std::string,std::vector<std::vector<double> > > {"tunit", "K"},
		 {"d", _dens},
		//std::pair<std::string,std::vector<std::vector<double> > > {"dunit", "kg/m3"},
		 {"e", energietot},
		//std::pair<std::string,std::vector<std::vector<double> > > {"etotunit", "J/kg"},
		 {"p", pressiontot},
		//std::pair<std::string,std::vector<std::vector<double> > > {"ptotunit", "Pa"}
	};

    std::cout << "File loaded: " << filename << std::endl;
}
/*----------------------------------------------------------------------------------------------*/
/* calculPetT : renvoie la pression et la température en fonction de la densité et de l énergie */
/*----------------------------------------------------------------------------------------------*/
bool TabuleeEOSService::calculPetT(double d0, double e0, double &t0, double &p0, std::string interp = "lin" ) {
  std::vector<std::vector<double>> _d(data["d"].begin(), data["d"].end());
  std::vector<std::vector<double>> _t(data["t"].begin(), data["t"].end());
  std::vector<std::vector<double>> p(data["p"].begin(), data["p"].end());
  std::vector<std::vector<double>> e(data["e"].begin(), data["e"].end());
  bool blocage = true;

  auto d =_d.at(0);
  auto t =_t.at(0);
  
    
//   std::cout << " densité " << d0 << std::endl;
//   std::cout << " energie " << e0 << " emax " << emax << std::endl;
//   std::cout << "data input d: " << d.at(1) << std::endl;
//   std::cout << "data input d : " << d.at(75) << std::endl;
//   std::cout << "data input p : " << p[5][10] << std::endl;
//   std::cout << "data input e: " << e[5][10] << std::endl;

  
  if (blocage) {
    if (e0 > emax || e0 < emin) {
        std::cout << "test in calculPetT energy blocking: " << e0 << " " << emin << " " << emax << std::endl;
        e0 = std::max(std::min(e0, emax), emin);
        exit(1);
    }
    if (d0 > dmax || d0 < dmin) {
        std::cout << "test in calculPetT density blocking: " << d0 << " " << dmin << " " << dmax << std::endl;
        d0 = std::max(std::min(d0, dmax), dmin);
        exit(1);
    }
   }


   if ((e0 - emax) * (e0 - emin) > 0) {
    return false;
   }

   if ((d0 - dmin) * (d0 - dmax) > 0) {
    return false;
   }

   double ad,ae;
   int jt, jd;

   for (jd = 0; jd < d.size() - 1; ++jd) {
    if ((d0 - d[jd]) * (d0 - d[jd + 1]) <= 0) {
        break;
    }
   }
   
   ad = (d0 - d[jd]) / (d[jd + 1] - d[jd]);
   

  std::vector<double> ed(t.size());
 
  for (jt = 0; jt < t.size() - 1; ++jt) {
    ed[jt] =  e[jd][jt] + ad*(e[jd+1][jt]-e[jd][jt]);
  }
   
  if (blocage) {
    auto it = std::minmax_element(ed.begin(), ed.end());
    double ed_min = *it.first;
    double ed_max = *it.second;
    if (e0 > ed_max || e0 < ed_min) {
        std::cout << " minmax_element : test in calculPetT energy blocking: " << e0 << " " << ed_min << " " <<  ed_max << std::endl;
        std::cout << " densité " << d0 << std::endl;
        e0 = std::max(std::min(e0,ed_max ),ed_min );
        exit(1);
    }
  }	
  for (jt = 0; jt < t.size() - 1; ++jt) {
    if ((e0 - ed[jt]) * (e0 - ed[jt + 1]) <= 0) {
        break;
    }
  }

  ae = (e0 - ed[jt]) / (ed[jt + 1] - ed[jt]);

  t0 = t[jt] + ae * (t[jt + 1] - t[jt]);
  p0 = (1.0 - ae) * (1.0 - ad) * p[jd][jt] + (1.0 - ae) * ad * p[jd+1][jt] + ae * (1.0 - ad) * p[jd][jt+1] + ae * ad * p[jd + 1][jt + 1];

  // std::cout << "calculde output: true " << t0 << " " << p0 << std::endl;
  return true; 
}

/*----------------------------------------------------------------------------------------------*/
/* calculPetE : renvoie la pression et l energie en fonction de la densité et de la température */
/*----------------------------------------------------------------------------------------------*/
bool TabuleeEOSService::calculPetE(double d0, double t0, double& p0, double& e0, std::string interp = "lin") {
  std::vector<std::vector<double>> _d(data["d"].begin(), data["d"].end());
  std::vector<std::vector<double>> _t(data["t"].begin(), data["t"].end());
  std::vector<std::vector<double>> p(data["p"].begin(), data["p"].end());
  std::vector<std::vector<double>> e(data["e"].begin(), data["e"].end());
  bool blocage = true;

  auto d =_d.at(0);
  auto t =_t.at(0);
  
//   std::cout << " densité " << d0 << std::endl;
//   std::cout << " température " << t0 << std::endl;
//   std::cout << "data input d: " << d.at(1) << std::endl;
//   std::cout << "data input d : " << d.at(75) << std::endl;
//   std::cout << "data input e : " << e[1][10] << std::endl;
//   std::cout << "data input e: " << e[5][10] << std::endl;
  
  if (blocage) {
    if (t0 > tmax || t0 < tmin) {
        std::cout << "test in calculPetE temperature blocking: " << t0 << " " << tmin << " " << tmax << std::endl;
        t0 = std::max(std::min(t0, tmax), tmin);
        exit(1);
    }
    if (d0 > dmax || d0 < dmin) {
        std::cout << "test in calculPetE density blocking: " << d0 << " " << dmin << " " << dmax << std::endl;
        d0 = std::max(std::min(d0, dmax), dmin);
        exit(1);
    }
   }

   if ((t0 - tmin) * (t0 - tmax) > 0) {
     return false;
   }

   if ((d0 - dmin) * (d0 - dmax) > 0) {
     return false;
   }

   int jt;
   for (jt = 0; jt < t.size() - 1; ++jt) {
    if ((t0 - t[jt]) * (t0 - t[jt + 1]) <= 0) {
        break;
    }
   }

   int jd;
   for (jd = 0; jd < d.size() - 1; ++jd) {
    if ((d0 - d[jd]) * (d0 - d[jd + 1]) <= 0) {
        break;
    }
   }

   double at, ad;
   if (interp == "lin") {
    at = (t0 - t[jt]) / (t[jt + 1] - t[jt]);
    ad = (d0 - d[jd]) / (d[jd + 1] - d[jd]);
   } else if (interp == "log") {
    at = (math::log(t0) - math::log(t[jt])) / (math::log(t[jt + 1]) - math::log(t[jt]));
    ad = (math::log(d0) - math::log(d[jd])) / (math::log(d[jd + 1]) - math::log(d[jd]));
   }
    
   // std::cout << " at et ad "<< at << " "<<  at << std::endl;
   p0 = (1.0 - at) * (1.0 - ad) * p[jd][jt] + (1.0 - at) * ad * p[jd+1][jt] + at * (1.0 - ad) * p[jd][jt+1] + at * ad * p[jd + 1][jt + 1];
   e0 = (1.0 - at) * (1.0 - ad) * e[jd][jt] + (1.0 - at) * ad * e[jd+1][jt] + at * (1.0 - ad) * e[jd][jt+1] + at * ad * e[jd + 1][jt + 1];

   // std::cout << "calcultd output: true " << p0 << " " << e0 << std::endl;
   
   return true;
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real TabuleeEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real TabuleeEOSService::getTensionLimitCst(IMeshEnvironment* env) { return options()->limitTension();}
Real TabuleeEOSService::getSpecificHeatCst(IMeshEnvironment* env) { return options()->specificHeat();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_TABULEEEOS(Tabulee, TabuleeEOSService);
