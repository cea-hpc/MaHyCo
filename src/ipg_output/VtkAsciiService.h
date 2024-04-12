// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef IPG_OUTPUT_VTKASCIISERVICE_H
#define IPG_OUTPUT_VTKASCIISERVICE_H


#include "ipg_output/IIpgOutput.h"

#include "arccore/base/String.h"
#include "arccore/collections/Array.h"
#include "arcane/ItemTypes.h"

#include "ipg_output/VtkAscii_axl.h"

using namespace Arcane;

/**
 * Service d'écriture des sorties au format vtk
 */
class VtkAsciiService 
: public ArcaneVtkAsciiObject
{
public:
  /** Constructeur de la classe */
  VtkAsciiService(const ServiceBuildInfo & sbi)
    : ArcaneVtkAsciiObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~VtkAsciiService() {};

public:
  void initOutput();
  void writeOutput(ParticleGroup particles);

private:
  /** liste des variables à sortir Real */
  SharedArray<String> m_variable_list_real;
  /** liste des variables à sortir Real3 */
  SharedArray<String> m_variable_list_real3;
  /** famille de particules*/
  IItemFamily* m_item_family;
};

#endif  // IPG_OUTPUT_VTKASCIISERVICE_H
