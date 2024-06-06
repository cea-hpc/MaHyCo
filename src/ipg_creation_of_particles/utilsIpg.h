// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

/*---------------------------------------------------------------------------*/
/*
  Verify if the given coordinates belongs (geometrically) to the given cell.
 */
/*---------------------------------------------------------------------------*/

/* bool isCoordInCell(Real3 particule_coord, Real3Array nodes_coord){ */
bool isCoordInCell(Real3 particule_coord, UniqueArray<Real3> nodes_coord){

  bool incell=0;

  // We use the following method :
  // let start from the particle coordinates and draw an infinite horizontal line to the right of the point.
  // The point belongs to the polygon if the line cross a face of the polygon an odd number of times
  // or if the point belongs to an edge/corner of the polygon.
  // otherwise, it is outside.


  // TODO: 3D version (le code ci-dessous n'est valable qu'en 2D dans le plan (x,y)
  // A vérifier: il faut juste ajouter if( z_particle>zmin_face et z_particle<=zmax_face ) ?

  // TODO: this method is identical for all service which create particles.
  // So, it would be better to write it only once in a separate file.

  // TODO: cas particulier ou la particule n'est attribuée à aucune cell:
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







/* below: don't work here because we don't have access to *m_node_coord* in this scope and I don't know how to have access to it here. */
/* I think that if I create an interface and service, il will be able to use the variable *m_node_coord*, but I don't know otherise.
Maybe something using findVariable, but I don't have succeded with that approach */

/* /\*---------------------------------------------------------------------------*\/ */
/* /\* */
/*   Verify if the given coordinates belongs (geometrically) to the given cell. */
/*  *\/ */
/* /\*---------------------------------------------------------------------------*\/ */

/* bool isCoordInCell(Real3 particule_coord, Cell cell){ */

/*   bool incell=0; */

/*   // We use the following method : */
/*   // let start from the particle coordinates and draw an infinite horizontal line to the right of the point. */
/*   // The point belongs to the polygon if the line cross a face of the polygon an odd number of times */
/*   // or if the point belongs to an edge/corner of the polygon. */
/*   // otherwise, it is outside. */


/*   // TODO: 3D version (le code ci-dessous n'est valable qu'en 2D dans le plan (x,y) */
/*   // A vérifier: il faut juste ajouter if( z_particle>zmin_face et z_particle<=zmax_face ) ? */

/*   // TODO: this method is identical for all service which create particles. */
/*   // So, it would be better to write it only once in a separate file. */

/*   // TODO: cas particulier ou la particule n'est attribuée à aucune cell: */
/*   // la particule est sur la face horizontale d'un élément qui n'a pas d'autre élément en dessous de lui. */

/*   /\* m_item_family = mesh()->findItemFamily(eItemKind::IK_Particle, "ActiveParticles"); *\/ */
  
/*   /\* IVariable* ivar = mesh()->findItemFamily(eItemKind::IK_Particle, "ActiveParticles")->findVariable("node_coord"); *\/ */

/*   // initialisation: position du dernier node */
/*   Real x1=m_node_coord[cell.node(cell.nbNode()-1)].x, y1=m_node_coord[cell.node(cell.nbNode()-1)].y; */

/*   // pour tous les nodes */
/*   for ( Integer inode = 0; inode < cell.nbNode(); ++inode ) { */

/*     // on recupere la position du node */
/*     Real x2=m_node_coord[cell.node(inode)].x, y2=m_node_coord[cell.node(inode)].y; */
/*     // on connait alors la face (x1,y1) -- (x2,y2) */

/*     // on vérifie si la particule est entre le ymin et le ymax de la face */
/*     if (particule_coord.y > math::min(y1,y2) ) */
/*       if (particule_coord.y <= math::max(y1,y2) ){ */

/*         if (particule_coord.x <= math::max(x1,x2) ){  // si la particule est à gauche de la face */

/*           double x_intersection = (x2-x1)/(y2-y1)*(particule_coord.y-y1) + x1;  // intersection entre la face et la ligne horizontale à droite de la particule */
/*           // note: le cas particulier y1==y2 n'est pas possible puisque y_particule > y1 et y_particule <= y2 */

/*           if ((particule_coord.x <= x_intersection) || (x1==x2) ) // on vérifie s'il y a intersection (ou si la face est vertical) */
/*             incell = !incell;  // true pour un nombre impair d'intersection */
/*         } */
/*       } */

/*     // on va passer à la face suivante */
/*     x1=x2; */
/*     y1=y2; */
    
/*   } */

/*   return incell; */
/* } */


