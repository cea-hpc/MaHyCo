// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
/*---------------------------------------------------------------------------*/
/* FaceDirectionMng.cc                                         (C) 2000-2012 */
/*                                                                           */
/* Infos sur les faces d'une direction X Y ou Z d'un maillage structuré.     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/utils/ArcanePrecomp.h"

#include "arcane/utils/FatalErrorException.h"

#include "arcane/IItemFamily.h"
#include "arcane/ItemGroup.h"

#include "cartesian/FaceDirectionMng.h"
#include "cartesian/ICartesianMesh.h"
#include "cartesian/FactCartDirectionMng.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class FaceDirectionMng::Impl
{
 public:
  Impl() : m_cartesian_mesh(0){}
 public:
  ItemGroup m_inner_all_items;
  ItemGroup m_outer_all_items;
  ItemGroup m_all_items;
  ICartesianMesh* m_cartesian_mesh;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FaceDirectionMng::
FaceDirectionMng()
: m_direction(MD_DirInvalid)
, m_dir(Integer(m_direction))
, m_p (nullptr)
{
  m_null_cell_internal=m_null_cell.internal();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FaceDirectionMng::
FaceDirectionMng(const FaceDirectionMng& rhs)
: m_direction(rhs.m_direction)
, m_dir(rhs.m_dir)
, m_p(rhs.m_p)
, m_p_cart_numb_face_dir(rhs.m_p_cart_numb_face_dir)
, m_p_numb_conv(rhs.m_p_numb_conv)
, m_cell_internals(rhs.m_cell_internals)
, m_nfacesm1_dir(rhs.m_nfacesm1_dir)
, m_delta_prev_cell(rhs.m_delta_prev_cell)
{
  m_null_cell_internal=m_null_cell.internal();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FaceDirectionMng::
~FaceDirectionMng()
{
  // Ne pas détruire le m_p.
  // Le gestionnnaire le fera via destroy()
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FaceDirectionMng::
init(ICartesianMesh* cm,eMeshDirection dir,
    FaceDirectionMng::CartesianGrid* p_cart_grid, FaceDirectionMng::NumbConvFaceCell* p_numb_conv)
{
  if (m_p)
    throw FatalErrorException(A_FUNCINFO,"Initialisation already done");
  m_p = new Impl();
  m_direction = dir;
  m_dir = Integer(m_direction);
  m_p->m_cartesian_mesh = cm;
  m_p_cart_numb_face_dir = p_cart_grid->cartNumFacePtr(m_dir);
  m_p_numb_conv = p_numb_conv;
  m_cell_internals = cm->mesh()->itemsInternal(IK_Cell);
  m_nfacesm1_dir = m_p_cart_numb_face_dir->nbItemDir(m_dir)-1;
  m_delta_prev_cell = p_cart_grid->cartNumCell().deltaDir(m_dir);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FaceDirectionMng::
destroy()
{
  delete m_p;
  m_p = nullptr;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FaceDirectionMng::
computeInnerAndOuterItems(const ItemGroup& items)
{
  IItemFamily* family = items.itemFamily();

  Int32UniqueArray inner_lids;
  Int32UniqueArray outer_lids;
  ENUMERATE_FACE(iitem,items){
    Int32 lid = iitem.itemLocalId();
    Face face = *iitem;
    if (face.nbCell()==1)
      outer_lids.add(lid);
    else
      inner_lids.add(lid);
  }
  int dir = (int)m_direction;
  m_p->m_inner_all_items = family->createGroup(String("AllInnerDirection")+dir,inner_lids,true);
  m_p->m_outer_all_items = family->createGroup(String("AllOuterDirection")+dir,outer_lids,true);
  m_p->m_all_items = items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FaceDirectionMng::
computeInnerAndOuterItemsFromFamily(IItemFamily* face_family)
{
  int dir = (int)m_direction;

  // On s'appuie sur l'implémentation CARTESIENNE pour récupérer les ids
  FactCartDirectionMng fact_cart_dm(m_p->m_cartesian_mesh->mesh());
  CartFaceDirectionMng&& cart_face_dm = fact_cart_dm.faceDirection(dir);

  // D'abord les faces intérieures
  auto&& cart_inner_group = cart_face_dm.innerFaces();
  // Les tableaux des ids peuvent etre pré-alloués
  Int32UniqueArray inner_lids(cart_inner_group.size());
  Int32 inner_cnt=0;
  ENUMERATE_AUTO_FACE(face_i, cart_inner_group) {
    inner_lids[inner_cnt++]=face_i.localId();
  }
  ARCANE_ASSERT(inner_lids.size()==inner_cnt, ("Le tableau inner_lids n'est pas intégralement rempli"));

  // Puis les faces extérieures
  auto&& prev_outer_group = cart_face_dm.previousOuterFaces();
  auto&& next_outer_group = cart_face_dm.nextOuterFaces();
  // Les tableaux des ids peuvent etre pré-alloués
  Int32UniqueArray outer_lids(prev_outer_group.size()+next_outer_group.size());
  Int32 outer_cnt=0;
  ENUMERATE_AUTO_FACE(face_i, prev_outer_group) {
    outer_lids[outer_cnt++]=face_i.localId();
  }
  ENUMERATE_AUTO_FACE(face_i, next_outer_group) {
    outer_lids[outer_cnt++]=face_i.localId();
  }
  ARCANE_ASSERT(outer_lids.size()==outer_cnt, ("Le tableau outer_lids n'est pas intégralement rempli"));

  // Toutes les faces
  auto&& cart_all_group = cart_face_dm.allFaces();
  // Les tableaux des ids peuvent etre pré-alloués
  Int32UniqueArray all_lids(cart_all_group.size());
  Int32 all_cnt=0;
  ENUMERATE_AUTO_FACE(face_i, cart_all_group) {
    all_lids[all_cnt++]=face_i.localId();
  }
  ARCANE_ASSERT(all_lids.size()==all_cnt, ("Le tableau all_lids n'est pas intégralement rempli"));
  ARCANE_ASSERT(all_lids.size()==(inner_lids.size()+outer_lids.size()), ("size(all_lids) != size(inner_lids+outer_lids)"));

  m_p->m_inner_all_items = face_family->createGroup(String("AllInnerDirection")+dir,inner_lids,true);
  m_p->m_outer_all_items = face_family->createGroup(String("AllOuterDirection")+dir,outer_lids,true);
  m_p->m_all_items = face_family->createGroup(String("AllFacesDirection")+dir,all_lids,true);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FaceGroup FaceDirectionMng::
allFaces() const
{
  return m_p->m_all_items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FaceGroup FaceDirectionMng::
innerFaces() const
{
  return m_p->m_inner_all_items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

FaceGroup FaceDirectionMng::
outerFaces() const
{
  return m_p->m_outer_all_items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
