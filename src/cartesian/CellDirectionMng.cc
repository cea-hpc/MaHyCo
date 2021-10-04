// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
/*---------------------------------------------------------------------------*/
/* CellDirectionMng.cc                                         (C) 2000-2017 */
/*                                                                           */
/* Infos sur les mailles d'une direction X Y ou Z d'un maillage structuré.   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/utils/ArcanePrecomp.h"

#include "arcane/utils/FatalErrorException.h"
#include "arcane/utils/ArgumentException.h"
#include "arcane/utils/ITraceMng.h"

#include "arcane/IItemFamily.h"
#include "arcane/ItemGroup.h"
#include "arcane/IMesh.h"

#include "cartesian/CellDirectionMng.h"
#include "cartesian/ICartesianMesh.h"
#include "cartesian/FactCartDirectionMng.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class CellDirectionMng::Impl
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

CellDirectionMng::
CellDirectionMng()
: m_direction(MD_DirInvalid)
, m_dir(Integer(m_direction))
, m_p(nullptr)
, m_next_face_index(-1)
, m_previous_face_index(-1)
, m_sub_domain_offset(-1)
, m_own_nb_cell(-1)
, m_global_nb_cell(-1)
, m_own_cell_offset(-1)
{
  m_null_cell_internal=m_null_cell.internal();
  for( Integer i=0; i<MAX_NB_NODE; ++i )
    m_nodes_indirection[i] = (-1);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

CellDirectionMng::
CellDirectionMng(const CellDirectionMng& rhs)
: m_direction(rhs.m_direction)
, m_dir(rhs.m_dir)
, m_p(rhs.m_p)
, m_next_face_index(rhs.m_next_face_index)
, m_previous_face_index(rhs.m_previous_face_index)
, m_sub_domain_offset(rhs.m_sub_domain_offset)
, m_own_nb_cell(rhs.m_own_nb_cell)
, m_global_nb_cell(rhs.m_global_nb_cell)
, m_own_cell_offset(rhs.m_own_cell_offset)
, m_p_cart_numb_cell(rhs.m_p_cart_numb_cell)
, m_cell_internals(rhs.m_cell_internals)
, m_ncellsm1_dir(rhs.m_ncellsm1_dir)
, m_delta_dir(rhs.m_delta_dir)
{
  m_null_cell_internal=m_null_cell.internal();
  for( Integer i=0; i<MAX_NB_NODE; ++i )
    m_nodes_indirection[i] = rhs.m_nodes_indirection[i];
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

CellDirectionMng::
~CellDirectionMng()
{
  // Ne pas détruire le m_p.
  // Le gestionnnaire le fera via destroy()
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CellDirectionMng::
init(ICartesianMesh* cm,eMeshDirection dir,CellDirectionMng::CartesianNumbering* p_cart_numb_cell)
{
  if (m_p)
    ARCANE_FATAL("Initialisation already done");
  m_p = new Impl();
  m_direction = dir;
  m_dir = Integer(m_direction);
  m_p->m_cartesian_mesh = cm;
  m_p_cart_numb_cell = p_cart_numb_cell;
  m_cell_internals = cm->mesh()->itemsInternal(IK_Cell);
  m_ncellsm1_dir = m_p_cart_numb_cell->nbItemDir(m_dir)-1;
  m_delta_dir = m_p_cart_numb_cell->deltaDir(m_dir);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CellDirectionMng::
destroy()
{
  delete m_p;
  m_p = nullptr;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CellDirectionMng::
computeInnerAndOuterItems(const ItemGroup& items)
{
  int dir = (int)m_direction;
  // On s'appuie sur l'implémentation CARTESIENNE pour récupérer les ids
  FactCartDirectionMng fact_cart_dm(m_p->m_cartesian_mesh->mesh());
  CartCellDirectionMng&& cart_cell_dm = fact_cart_dm.cellDirection(dir);

  // D'abord les mailles intérieures
  auto&& cart_inner_group = cart_cell_dm.innerCells();
  // Les tableaux des ids peuvent etre pré-alloués
  Int32UniqueArray inner_lids(cart_inner_group.size());
  Int32 inner_cnt=0;
  ENUMERATE_AUTO_CELL(cell_i, cart_inner_group) {
    inner_lids[inner_cnt++]=cell_i.localId();
  }
  ARCANE_ASSERT(inner_lids.size()==inner_cnt, ("Le tableau inner_lids n'est pas intégralement rempli"));

  // Puis les mailles extérieures
  auto&& prev_outer_group = cart_cell_dm.previousOuterCells();
  auto&& next_outer_group = cart_cell_dm.nextOuterCells();
  // Les tableaux des ids peuvent etre pré-alloués
  Int32UniqueArray outer_lids(prev_outer_group.size()+next_outer_group.size());
  Int32 outer_cnt=0;
  ENUMERATE_AUTO_CELL(cell_i, prev_outer_group) {
    outer_lids[outer_cnt++]=cell_i.localId();
  }
  ENUMERATE_AUTO_CELL(cell_i, next_outer_group) {
    outer_lids[outer_cnt++]=cell_i.localId();
  }
  ARCANE_ASSERT(outer_lids.size()==outer_cnt, ("Le tableau outer_lids n'est pas intégralement rempli"));
  IItemFamily* family = items.itemFamily();
  m_p->m_inner_all_items = family->createGroup(String("AllInnerDirection")+dir,inner_lids,true);
  m_p->m_outer_all_items = family->createGroup(String("AllOuterDirection")+dir,outer_lids,true);
  m_p->m_all_items = items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

CellGroup CellDirectionMng::
allCells() const
{
  return m_p->m_all_items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

CellGroup CellDirectionMng::
innerCells() const
{
  return m_p->m_inner_all_items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

CellGroup CellDirectionMng::
outerCells() const
{
  return m_p->m_outer_all_items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CellDirectionMng::
setNodesIndirection(Int32ConstArrayView nodes_indirection)
{
  for( Integer i=0; i<MAX_NB_NODE; ++i )
    m_nodes_indirection[i] = nodes_indirection[i];

  ITraceMng* tm = m_p->m_cartesian_mesh->traceMng();

  tm->info() << "Set computed indirection dir=" << (int)m_direction;
  for( Integer i=0; i<MAX_NB_NODE; ++i ){
    tm->info() << "Indirection i=" << i << " v=" << m_nodes_indirection[i];
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
