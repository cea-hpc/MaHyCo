// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Copyright 2000-2024 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0

#ifndef UTILSIPG_H
#define UTILSIPG_H

/*---------------------------------------------------------------------------*/
/*
  Verify if the given coordinates belongs (geometrically) to the given cell.

  WARNING: pour le moment la fonction n'est valable qu'en 2D dans le plan xy
 */
/*---------------------------------------------------------------------------*/

/* bool isCoordInCell(Real3 particule_coord, Real3Array nodes_coord){ */
bool isCoordInCell(Real3 particule_coord, UniqueArray<Real3> nodes_coord){

  bool incell=0;

  // We use a ray-casting method :
  // let start from the particle coordinates and draw an infinite horizontal line to the right of the point.
  // The point belongs to the polygon if the line cross a face of the polygon an odd number of times
  // or if the point belongs to an edge/corner of the polygon.
  // otherwise, it is outside.

  // TODO: version 3D et version 2D dans les plan (y,z) et (x,z)
  //       (le code ci-dessous n'est valable qu'en 2D dans le plan (x,y)

  // TODO: cas particulier où la particule n'est attribuée à aucune cell:
  // la particule est sur la face horizontale d'un élément qui n'a pas d'autre élément en dessous de lui.

  // initialisation: position du dernier node
  Real x1=nodes_coord[nodes_coord.size()-1].x, y1=nodes_coord[nodes_coord.size()-1].y;

  // pour tous les nodes
  for ( Integer inode = 0; inode < nodes_coord.size(); ++inode ) {

    // on recupere la position du node
    Real x2=nodes_coord[inode].x, y2=nodes_coord[inode].y;
    // on connait alors la face (x1,y1) -- (x2,y2)

    // on vérifie si la particule est entre le ymin et le ymax de la face
    if (particule_coord.y > math::min(y1,y2) )
      if (particule_coord.y <= math::max(y1,y2) ){

        if (particule_coord.x <= math::max(x1,x2) ){  // si la particule est à gauche de la face

          double x_intersection = (x2-x1)/(y2-y1)*(particule_coord.y-y1) + x1;  // intersection entre la face et la ligne horizontale à droite de la particule
          // note: le cas particulier y1==y2 n'est pas possible puisque y_particule > y1 et y_particule <= y2

          if ((particule_coord.x <= x_intersection) || (x1==x2) ) // on vérifie s'il y a intersection (ou si la face est vertical)
            incell = !incell;  // true pour un nombre impair d'intersection
        }
      }

    // on va passer à la face suivante
    x1=x2;
    y1=y2;
    
  }

  return incell;
}




// fixme: below, I miss the correct include to have access to m_particle_coor[] and m_node_coord[]

/* /\*---------------------------------------------------------------------------*\/ */
/* /\* */
/*   Associe une particule à un noeud (== à une cellule duale) */
/*  *\/ */
/* /\*---------------------------------------------------------------------------*\/ */

/* Node findNodeOfParticle(Particle particule){ */

/*   /\* /\\* IItemFamily* item_family = mesh()->findItemFamily (eItemKind::IK_Particle, "AllParticles"); *\\/ *\/ */
/*   /\* /\\* IVariable* ivar = item_family->findVariable("ParticleCoordinates", true); *\\/ *\/ */
/*   /\* IVariable* ivar = mesh()->findItemFamily(eItemKind::IK_Particle, "AllParticles")->findVariable("ParticleCoordinates", true); *\/ */
/*   /\* VariableParticleReal3 m_particle_coord(ivar); *\/ */
  
  
/*   Cell cell = particule.cell(); */

/*   // chercher à quel noeud appartient la particule */
/*   // pour le moment, on considère que la particule appartient au noeud le plus proche */
/*   Real distance_min=1e30; */
/*   Node node_proche; */
/*   for ( NodeEnumerator inode ( cell.nodes() ); inode.hasNext(); ++inode ) { */
/*     /\* Real distance = (node_coord[inode] - particle_coord[particule]).normL2(); *\/ */
/*     Real distance = (m_node_coord[inode] - m_particle_coord[particule]).normL2(); */
/*     if (distance < distance_min){ */
/*       distance_min = distance; */
/*       node_proche=*inode; */
/*     } */
/*   } */
/*   return node_proche; */
/* } */












#endif
