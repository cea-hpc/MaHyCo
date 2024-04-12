// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "VtkAsciiService.h"

#include <iostream>
#include <fstream>
#include <string>

#include "arcane/IItemFamily.h"
#include "arcane/IVariable.h"
#include "arcane/IMesh.h"
#include "arcane/VariableBuildInfo.h"


using namespace Arcane;


/*---------------------------------------------------------------------------*/
/*! \brief
 * Initialisation des sorties particules :
 * - chargement de la liste des variables à sortir
 */
/*---------------------------------------------------------------------------*/

void VtkAsciiService::initOutput()
{
  m_item_family = mesh()->findItemFamily(eItemKind::IK_Particle, "ActiveParticles");
  for (const String& var_name : options()->getVariable()) {
     IVariable* ivar = m_item_family->findVariable(var_name, true);
     // true => levée d'exception si la variable n'existe pas.
     if (ivar->dataType() == DT_Real) {
       m_variable_list_real.push_back(var_name);
     } else if (ivar->dataType() == DT_Real3) {
       m_variable_list_real3.push_back(var_name);
     }
  }
}

/*---------------------------------------------------------------------------*/
/*! \brief
 * Ecriture des sorties pour les particules en argument 
 *
 * @param particles : groupe de particule à sortir
 */
/*---------------------------------------------------------------------------*/

void VtkAsciiService::writeOutput(ParticleGroup particles)
{
  const Real time = m_global_time();
  ostringstream f;
  f << "output/" << options()->getFilename() << "_" << time << ".vtp";

  info()<< "écriture dans le fichier " << f.str();

  std::ofstream write_file(f.str());

  if (!write_file.is_open()) {
     std::cout << "ERREUR: Impossible d'ouvrir le fichier." << std::endl;
  }

  int size = particles.size();

  // Fonction d'écriture des variables réelles
  std::function<void(const String&)> write_variable_real = [&](const String& var_name){
    write_file << "      <DataArray type=\"Float32\" Name=\""<< var_name << "\" format=\"ascii\">\n";
    write_file << "      ";
    IVariable* ivar = m_item_family->findVariable(var_name);
    VariableParticleReal var(ivar);
    ENUMERATE_PARTICLE (i_part, particles) {
      write_file << var[i_part] << " ";
    }
    write_file << "\n      </DataArray>\n";
  };

  // Fonction d'écriture des variables Real3
  std::function<void(const String&)> write_variable_real3 = [&](const String& var_name){
    write_file << "      <DataArray type=\"Float32\" Name=\""<< var_name << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    write_file << "      ";
    IVariable* ivar = m_item_family->findVariable(var_name);
    VariableParticleReal3 var(ivar);
    ENUMERATE_PARTICLE (i_part, particles) {
      write_file << var[i_part].x << " " << var[i_part].y << " " << var[i_part].z << " \n";
    }
    write_file << "\n      </DataArray>\n";
  };

  // Ecriture du fichier de sortie
  write_file << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\"> \n";
  write_file << "<PolyData>\n";
  write_file << "<Piece NumberOfPoints=\"" << size << "\" NumberOfVerts=\"1\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
  write_file << "    <PointData Scalars=\"Temperature\" Vectors=\"Velocity\">\n";
  // champ vectoriel
    for (const String& name : m_variable_list_real3) {
    write_variable_real3(name);
  }
  // champ scalaire
  for (const String& name : m_variable_list_real) {
    write_variable_real(name);
  }
  write_file << "    </PointData>\n";
  write_file << "    <CellData>\n";
  write_file << "    </CellData>\n";
  write_file << "    <Points>\n";
 /* write_file << "      <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  // Coordonnées des particules
  write_file << "        0 0 0\n";
  write_file << "        1 1 0\n";
  write_file << "        2 1 0\n";
  write_file << "        2 4 0\n";
  write_file << "        2 2 0\n";
  write_file << "      </DataArray>\n";*/
  write_variable_real3("ParticleCoordinates");
  write_file << "    </Points>\n";
  write_file << "    <Verts>\n";
  write_file << "      <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  // liste des particules de 0 à size - 1
  write_file << "        ";
  for (int i=0 ; i < size; i++) {
    write_file << i << " ";
  }
  write_file << "\n      </DataArray>\n";
  write_file << "      <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  write_file << "       " << size << "\n";
  write_file << "      </DataArray>\n";
  write_file << "    </Verts>\n";
  write_file << "  </Piece>\n";
  write_file << "</PolyData>\n";
  write_file << "</VTKFile>\n";

  write_file.close();
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_VTKASCII(VtkAscii, VtkAsciiService);
