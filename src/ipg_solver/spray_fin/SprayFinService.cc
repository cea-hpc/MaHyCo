// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "SprayFinService.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/** 
    Initialisation solver particules
*/
/*---------------------------------------------------------------------------*/
void SprayFinService::initSolverParticles() {

  IItemFamily* item_family = mesh()->findItemFamily (eItemKind::IK_Particle, "AllParticles");
  m_particles_family = item_family->toParticleFamily();
  activeParticlesGroup = item_family->findGroup("activeItem");

}

/*---------------------------------------------------------------------------*/
/** 
    correction de la vitesse du fluide (prise en compte de la traînée due aux particules).
*/
/*---------------------------------------------------------------------------*/
void SprayFinService::correctFluidVelocity() {

  const Real dt (m_global_deltat());

  // 1 - trouver les particules appartenant à chaque cellule *duale*

  // 2 - calcul de la force de trainée aux noeuds

  // 3 - MaJ de la vitesse du fluide aux noeuds


  // // // // Arcane::geometric::GeomShapeMng shape_mng = mesh();
  // // Arcane::geometric::GeomShapeMng shape_mng(mesh(), "coucou");
  // // // // Arcane::geometric::GeomShapeMng shape_mng;
  // // // // GeomShapeMng shape_mng;
  // // Arcane::geometric::GeomShapeView shape_view;
  // // shape_mng.initialize(); // il faut initialize avant initShape
  
  // ENUMERATE_CELL(icell,allCells()){

  //   // // Arcane::geometric::GeomShapeMng shape_mng(mesh());
  //   // // Arcane::geometric::GeomShapeView shape_view;
  //   // // shape_mng.initialize(); // il faut initialize avant initShape
  //   // Cell cell = *icell;
  //   // // Initialise la vue \a shape_view sur la maille \a cell
  //   // shape_mng.initShape(shape_view,cell);

  //   // Cell cell = * icell;
  //   // for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) {
  //   //   info() << "Nodes= " << m_node_coord[inode];
  //   // }   
    
   

  // }

  // // VariableNodeReal3& node_coord = {0.,0.,0.};
  // // ENUMERATE_CELL(icell,allCells()){
  // //   Cell cell = *icell;
  // //   info() << "Node0=" << node_coord[cell.node(0)];
  // //   info() << "Node2=" << node_coord[cell.node(2)];
  // //   Real3 middle = (node_coord[cell.node(3)] + node_coord[cell.node(4)]) / 2.0;
  // // }
  
  // // // Arcane::geometric::GeomShapeMng& shape_mng = *mesh();
  // // Arcane::geometric::GeomShapeMng shape_mng = mesh();
  // // // Arcane::geometric::GeomShapeMng shape_mng(mesh());
  // // Arcane::geometric::GeomShapeView shape;

  // // Real3 my_point={0., 0., 0.};
  
  // // ENUMERATE_CELL(icell, allCells()){

  // //   shape_mng.initShape(shape,*icell);

  // //   info() << "mon point " << my_point << " is inside the cell " << shape.node(0) << shape.node(1) << shape.node(2) << shape.node(3) << " ? Answer: "; 

  // //   // bool toto = ProjectionInfo::isInside(shape.node(0),shape.node(1),shape.node(2),my_point);
  // //   // bool toto = isInside(shape.node(0),shape.node(1),shape.node(2),my_point);

  // //   // info() << "mon point " << my_point << " is inside the cell " << shape.node(0) << shape.node(1) << shape.node(2) << shape.node(3) << " ? Answer: " << (isInside(shape.node(0),shape.node(1),shape.node(2),my_point) || isInside(shape.node(3),shape.node(1),shape.node(2),my_point)); 

  // //   // Real3 middle = (shape.node(3) + shape.node(4)) / 2.0;

  // // }


  // ENUMERATE_PARTICLE (part_i, activeParticlesGroup) {

  //   info() << "Particle " << part_i.localId() << " has cell : " << part_i->hasCell() ;

  //   Particle particule = *part_i;

  //   // Cell cell = * icell;

  //   // particule.cell();
    
  //   // setParticleCell(particule, cellule);
    
  //   info() << "Particle " << part_i.localId() << " has cell : " << part_i->hasCell() ;
    
  // }
  
  // ENUMERATE_NODE ( inode, allNodes() ) {
  //   Node node = *inode;
  //   m_velocity[inode] += (dt / m_node_mass[inode]) * m_force[inode];
  // }
  // q. est-ce que je peux erase et utiliser m_force ?
  // ou est-ce que je dois creer une nouvelle variable ?

  
}

/*---------------------------------------------------------------------------*/
/** 
    mise à jour de la vitesse des particules
*/
/*---------------------------------------------------------------------------*/
void SprayFinService::updateParticleVelocity() {

}

/*---------------------------------------------------------------------------*/
/** 
      mise à jour de la position des particules
*/
/*---------------------------------------------------------------------------*/
void SprayFinService::updateParticlePosition() {

  // we need dt^{n+1/2}
  const Real dt ( 0.5 * ( m_global_old_deltat() + m_global_deltat() ) );
  
  ENUMERATE_PARTICLE (part_i, activeParticlesGroup) {
    m_particle_coord[part_i] += m_particle_velocity[part_i] * dt;
    info() << "Particle " << part_i.localId() << " vel = " << m_particle_velocity[part_i];
    info() << "Particle " << part_i.localId() << " coord = " << m_particle_coord[part_i];
  }
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_SPRAYFIN(SprayFin, SprayFinService);
