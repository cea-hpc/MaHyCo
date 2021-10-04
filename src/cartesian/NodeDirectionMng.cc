// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
/*---------------------------------------------------------------------------*/
/* NodeDirectionMng.cc                                         (C) 2000-2012 */
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

#include "cartesian/NodeDirectionMng.h"
#include "cartesian/ICartesianMesh.h"
#include "cartesian/FactCartDirectionMng.h"
#include <algorithm>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Cartesian {

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class NodeDirectionMng::Impl
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

NodeDirectionMng::
NodeDirectionMng()
: m_direction(MD_DirInvalid)
, m_dir(Integer(m_direction))
, m_p(nullptr)
{
  m_null_node_internal=m_null_node.internal();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

NodeDirectionMng::
NodeDirectionMng(const NodeDirectionMng& rhs)
: m_direction(rhs.m_direction)
, m_dir(rhs.m_dir)
, m_p(rhs.m_p)
, m_p_cart_numb_node(rhs.m_p_cart_numb_node)
, m_node_internals(rhs.m_node_internals)
, m_nnodesm1_dir(rhs.m_nnodesm1_dir)
, m_delta_dir(rhs.m_delta_dir)
{
  m_null_node_internal=m_null_node.internal();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

NodeDirectionMng::
~NodeDirectionMng()
{
  // Ne pas détruire le m_p.
  // Le gestionnnaire le fera via destroy()
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NodeDirectionMng::
init(ICartesianMesh* cm,eMeshDirection dir,NodeDirectionMng::CartesianNumbering* p_cart_numb_node)
{
  if (m_p)
    ARCANE_FATAL("Initialisation already done");
  m_p = new Impl();
  m_direction = dir;
  m_dir = Integer(m_direction);
  m_p->m_cartesian_mesh = cm;
  m_p_cart_numb_node = p_cart_numb_node;
  m_node_internals = cm->mesh()->itemsInternal(IK_Node);
  m_nnodesm1_dir = m_p_cart_numb_node->nbItemDir(m_dir)-1;
  m_delta_dir = m_p_cart_numb_node->deltaDir(m_dir);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NodeDirectionMng::
destroy()
{
  delete m_p;
  m_p = nullptr;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void NodeDirectionMng::
computeInnerAndOuterItems(const ItemGroup& items)
{
  int dir = (int)m_direction;
  // On s'appuie sur l'implémentation CARTESIENNE pour récupérer les ids
  FactCartDirectionMng fact_cart_dm(m_p->m_cartesian_mesh->mesh());
  CartNodeDirectionMng&& cart_node_dm = fact_cart_dm.nodeDirection(dir);

  // D'abord les noeuds intérieurs
  auto&& cart_inner_group = cart_node_dm.innerNodes();
  // Les tableaux des ids peuvent etre pré-alloués
  Int32UniqueArray inner_lids(cart_inner_group.size());
  Int32 inner_cnt=0;
  ENUMERATE_AUTO_NODE(node_i, cart_inner_group) {
    inner_lids[inner_cnt++]=node_i.localId();
  }
  ARCANE_ASSERT(inner_lids.size()==inner_cnt, ("Le tableau inner_lids n'est pas intégralement rempli"));

  // Puis les noeuds extérieurs
  auto&& prev_outer_group = cart_node_dm.previousOuterNodes();
  auto&& next_outer_group = cart_node_dm.nextOuterNodes();
  // Les tableaux des ids peuvent etre pré-alloués
  Int32UniqueArray outer_lids(prev_outer_group.size()+next_outer_group.size());
  Int32 outer_cnt=0;
  ENUMERATE_AUTO_NODE(node_i, prev_outer_group) {
    outer_lids[outer_cnt++]=node_i.localId();
  }
  ENUMERATE_AUTO_NODE(node_i, next_outer_group) {
    outer_lids[outer_cnt++]=node_i.localId();
  }
  ARCANE_ASSERT(outer_lids.size()==outer_cnt, ("Le tableau outer_lids n'est pas intégralement rempli"));

  IItemFamily* family = items.itemFamily();
  m_p->m_inner_all_items = family->createGroup(String("AllInnerDirection")+dir,inner_lids,true);
  m_p->m_outer_all_items = family->createGroup(String("AllOuterDirection")+dir,outer_lids,true);
  m_p->m_all_items = items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

NodeGroup NodeDirectionMng::
allNodes() const
{
  return m_p->m_all_items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

NodeGroup NodeDirectionMng::
innerNodes() const
{
  return m_p->m_inner_all_items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

NodeGroup NodeDirectionMng::
outerNodes() const
{
  return m_p->m_outer_all_items;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
