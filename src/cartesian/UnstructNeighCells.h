#ifndef CARTESIAN_UNSTRUCT_NEIGH_CELLS_H
#define CARTESIAN_UNSTRUCT_NEIGH_CELLS_H

#include "arcane/Item.h"
#include "arcane/utils/ArcaneGlobal.h"
#include "arcane/utils/ITraceMng.h"
#include "cartesian/interface/CartesianConnectivity.h"

namespace Cartesian {

/*
 * \brief Caractéristiques pour CartNeighCells en fonction de la dimension
 */
template<Integer DIM>
class UnstructNeighCellsTraits {
 public:
  using eCellPos = Integer;
  static constexpr Integer stencil_sz = -1;
  using ArrayCellType = UniqueArray<Cell>;
 public:
  static void neighCells([[maybe_unused]] const Cell& cell, ArrayCellType& adj_cells, 
      [[maybe_unused]] const CartesianInterface::CartesianConnectivity& cc) {
    ARCANE_ASSERT(false, ("UnstructNeighCellsTraits<DIM>::neighCells(const Cell& cell, ...) doit être spécialisé"));
  }

  static void applyBC([[maybe_unused]] ArrayCellType& adj_cells, 
      [[maybe_unused]] ITraceMng* trace_mng) {
    ARCANE_ASSERT(false, ("UnstructNeighCellsTraits<DIM>::applyBC(...) doit être spécialisé"));
  }
};

// Spécialisation 2D
template<>
class UnstructNeighCellsTraits<2> {
  //! Positions relatives des mailles adjacentes par rapport à une maille centrale
  enum eCellPos2D {
    SW = 0,
    S  = 1,
    SE = 2,
    E  = 3,
    NE = 4,
    N  = 5,
    NW = 6,
    W  = 7,
    C  = 8,  //! C = Center
    Max = 9,
    Invalid = (-1)
  };
 public:
  using eCellPos = eCellPos2D;
  static constexpr Integer stencil_sz = (Integer)Max;
  using ArrayCellType = UniqueArray<Cell>;
 public:
  static void neighCells(const Cell& cell, ArrayCellType& adj_cells, const CartesianInterface::CartesianConnectivity& cc) {
    ARCANE_ASSERT(stencil_sz<=adj_cells.size(), ("adj_cells doit etre au moins de taille stencil_sz"));
    Node node_1 = cc.lowerLeft(cell),
         node_2 = cc.lowerRight(cell),
         node_3 = cc.upperRight(cell),
         node_4 = cc.upperLeft(cell);
    Cell cellSW = cc.lowerLeft(node_1),
         cellW = cc.upperLeft(node_1),
         cellNW = cc.upperLeft(node_4),
         cellNE = cc.upperRight(node_3),
         cellE = cc.lowerRight(node_3),
         cellSE = cc.lowerRight(node_2),
         cellS = cc.lowerRight(node_1),
         cellN = cc.upperRight(node_4);

    adj_cells[0] = cellSW,
      adj_cells[1] = cellS,
      adj_cells[2] = cellSE,
      adj_cells[3] = cellE,
      adj_cells[4] = cellNE,
      adj_cells[5] = cellN,
      adj_cells[6] = cellNW,
      adj_cells[7] = cellW;
    adj_cells[8] = cell;
  }

  static void applyBC(ArrayCellType& adj_cells, [[maybe_unused]] ITraceMng* trace_mng) {
    // Conditions aux limites
    using CP = eCellPos;
    // Bord du bas
    if (adj_cells[CP::S].null())
      adj_cells[CP::SW] = adj_cells[CP::W], adj_cells[CP::S] = adj_cells[CP::C], adj_cells[CP::SE] = adj_cells[CP::E];
    // Bord de droite
    if (adj_cells[CP::E].null())
      adj_cells[CP::SE] = adj_cells[CP::S], adj_cells[CP::E] = adj_cells[CP::C], adj_cells[CP::NE] = adj_cells[CP::N];
    // Bord du haut
    if (adj_cells[CP::N].null())
      adj_cells[CP::NE] = adj_cells[CP::E], adj_cells[CP::N] = adj_cells[CP::C], adj_cells[CP::NW] = adj_cells[CP::W];
    // Bord de gauche
    if (adj_cells[CP::W].null())
      adj_cells[CP::NW] = adj_cells[CP::N], adj_cells[CP::W] = adj_cells[CP::C], adj_cells[CP::SW] = adj_cells[CP::S];
    // Check
    ARCANE_ASSERT( (adj_cells[CP::SW].null() || adj_cells[CP::S].null() || adj_cells[CP::SE].null() || adj_cells[CP::E].null() ||
        adj_cells[CP::NE].null() || adj_cells[CP::N].null() || adj_cells[CP::NW].null() || adj_cells[CP::W].null()),
        ("Revoir l'algorithme : on ne devrait pas se retrouver ici.") );
  }
};

// Spécialisation 3D
template<>
class UnstructNeighCellsTraits<3> {
 public:
  /* Convention : xyz avec x, y, z \in {L, C, R}
   * L(eft)   = -1 dans la direction considérée
   * C(enter) =  0 dans la direction considérée
   * R(ight)  = +1 dans la direction considérée
   *
   *             X, Y, Z
   * Ex : LRC = -1,+1, 0
   */
  enum eCellPos {
    LLL = 0,
    CLL = 1,
    RLL = 2,
    LLC = 3,
    CLC = 4,
    RLC = 5,
    LLR = 6,
    CLR = 7,
    RLR = 8,
    LCL = 9,
    CCL = 10,
    RCL = 11,
    LCC = 12,
    CCC = 13,
    RCC = 14,
    LCR = 15,
    CCR = 16,
    RCR = 17,
    LRL = 18,
    CRL = 19,
    RRL = 20,
    LRC = 21,
    CRC = 22,
    RRC = 23,
    LRR = 24,
    CRR = 25,
    RRR = 26,
    //! C = Center, idem que CC, pourquoi ?
    C   = 27,  
    Max = 28,
    Invalid = (-1)
  };

  static constexpr Integer stencil_sz = (Integer)Max;
  using ArrayCellType = UniqueArray<Cell>;
 public:
  static void neighCells(const Cell& cell, ArrayCellType& adj_cells, const CartesianInterface::CartesianConnectivity& cc) {
    ARCANE_ASSERT(stencil_sz<=adj_cells.size(), ("adj_cells doit etre au moins de taille stencil_sz"));
    // Noeuds de la cellule
    const Node& node_0 = cc.lowerLeft(cell);
    const Node& node_1 = cc.lowerRight(cell);
    const Node& node_2 = cc.upperRight(cell);
    const Node& node_3 = cc.upperLeft(cell);
    const Node& node_4 = cc.topZLowerLeft(cell);
    const Node& node_5 = cc.topZLowerRight(cell);
    const Node& node_6 = cc.topZUpperRight(cell);
    const Node& node_7 = cc.topZUpperLeft(cell);

    // Source : hydrop_fvol_3d.f90
    // Plan avant
    // 1 - sud ouest
    adj_cells[0] = cc.lowerLeft(node_0);
    // 2 - sud
    adj_cells[1] = cc.lowerRight(node_0);
    // 3 - sud est
    adj_cells[2] = cc.lowerRight(node_1);
    // 4 - ouest
    adj_cells[3] = cc.topZLowerLeft(node_0);
    // 5 - centre
    adj_cells[4] = cc.topZLowerRight(node_0);
    // 6 -Est
    adj_cells[5] = cc.topZLowerRight(node_1);
    // 7 - nord ouest
    adj_cells[6] = cc.topZLowerLeft(node_4);
    // 8 - Nord
    adj_cells[7] = cc.topZLowerRight(node_4);
    // 9 - Noerd Est
    adj_cells[8] = cc.topZLowerRight(node_5);

    // Plan de la maille courante
    // 10 - sud ouest
    adj_cells[9] = cc.lowerLeft(node_3);
    // 11 - sud
    adj_cells[10] = cc.lowerRight(node_3);
    // 12 - sud est
    adj_cells[11] = cc.lowerRight(node_2);
    // 13 - ouest
    adj_cells[12] = cc.topZLowerLeft(node_3);
    // 14 - centre
    adj_cells[13] = cc.topZLowerRight(node_3);
    // 15 -Est
    adj_cells[14] = cc.topZLowerRight(node_2);
    // 16 - nord ouest
    adj_cells[15] = cc.topZLowerLeft(node_7);
    // 17 - Nord
    adj_cells[16] = cc.topZLowerRight(node_7);
    // 18 - Noerd Est
    adj_cells[17] = cc.topZLowerRight(node_6);

    // Plan arrière
    // 19 - sud ouest
    adj_cells[18] = cc.upperLeft(node_3);
    // 20 - sud
    adj_cells[19] = cc.upperRight(node_3);
    // 21 - sud est
    adj_cells[20] = cc.upperRight(node_2);
    // 22 - ouest
    adj_cells[21] = cc.topZUpperLeft(node_3);
    // 23 - centre
    adj_cells[22] = cc.topZUpperRight(node_3);
    // 24 -Est
    adj_cells[23] = cc.topZUpperRight(node_2);
    // 25 - nord ouest
    adj_cells[24] = cc.topZUpperLeft(node_7);
    // 26 - Nord
    adj_cells[25] = cc.topZUpperRight(node_7);
    // 27 - Noerd Est
    adj_cells[26] = cc.topZUpperRight(node_6);

    // Pourquoi ? Elle est déjà présente dans adj_cells[13] !?
    adj_cells[27] = cell;
  }

  static void applyBC(ArrayCellType& adj_cells, ITraceMng* ) {
    // Traitement des conditions aux bords
    // Bord du bas
    if (adj_cells[4].null()) {
      for (Integer ii = 0 ; ii < 9 ; ++ii)
        adj_cells[ii] = adj_cells[ii+9];
    }
    // Bord du haut
    if (adj_cells[22].null()) {
      for (Integer ii(18) ; ii < 27 ; ++ii) {
        adj_cells[ii] = adj_cells[ii-9];
      }
    }
    // Bord de droite
    if (adj_cells[14].null()) {
      adj_cells[2] = adj_cells[1];
      adj_cells[5] = adj_cells[4];
      adj_cells[8] = adj_cells[7];
      adj_cells[11] = adj_cells[10];
      adj_cells[14] = adj_cells[13];
      adj_cells[17] = adj_cells[16];
      adj_cells[20] = adj_cells[19];
      adj_cells[23] = adj_cells[22];
      adj_cells[26] = adj_cells[25];
    }
    // Bord de gauche
    if (adj_cells[12].null()) {
      adj_cells[0] = adj_cells[1];
      adj_cells[3] = adj_cells[4];
      adj_cells[6] = adj_cells[7];
      adj_cells[9] = adj_cells[10];
      adj_cells[12] = adj_cells[13];
      adj_cells[15] = adj_cells[16];
      adj_cells[18] = adj_cells[19];
      adj_cells[21] = adj_cells[22];
      adj_cells[24] = adj_cells[25];
    }
    // Bord de devant
    if (adj_cells[10].null()) {
      adj_cells[0] = adj_cells[3];
      adj_cells[1] = adj_cells[4];
      adj_cells[2] = adj_cells[5];
      adj_cells[9] = adj_cells[12];
      adj_cells[10] = adj_cells[13];
      adj_cells[11] = adj_cells[14];
      adj_cells[18] = adj_cells[21];
      adj_cells[19] = adj_cells[22];
      adj_cells[20] = adj_cells[23];
    }
    // Bord de derrière
    if (adj_cells[16].null()) {
      adj_cells[6] = adj_cells[3];
      adj_cells[7] = adj_cells[4];
      adj_cells[8] = adj_cells[5];
      adj_cells[15] = adj_cells[12];
      adj_cells[16] = adj_cells[13];
      adj_cells[17] = adj_cells[14];
      adj_cells[24] = adj_cells[21];
      adj_cells[25] = adj_cells[22];
      adj_cells[26] = adj_cells[23];
    }
  }
};

/*---------------------------------------------------------------------------*/
/*!
 * \brief
 * Permet de retrouver les mailles voisines d'une maille (implémentation non-structurée)
 */
/*---------------------------------------------------------------------------*/
template<Integer DIM>
class UnstructNeighCells {
 public:
  using CellType = Cell;
  using CellEnumeratorType = CellEnumerator;

  //! Type du tableau rempli par neighCells
  using ArrayCellType = UniqueArray<CellType>;

  //! Alias
  using Traits = UnstructNeighCellsTraits<DIM>;

  //! Positions relatives des mailles adjacentes par rapport à une maille centrale
  using eCellPos = typename Traits::eCellPos;

  //! Taille à laquelle doit être alloué le tableau pour neighCells
  static constexpr Integer stencil_sz = Traits::stencil_sz;

 public:
  UnstructNeighCells(const CartesianInterface::CartesianConnectivity& cartesian_connectivity, ITraceMng* trace_mng) 
  : m_cart_conn (cartesian_connectivity), m_trace_mng(trace_mng)
  {
  }

  //! Récupère les ids locaux des mailles adjacente à une maille 
  // Vaut -1 si la maille voisine n'existe pas à la position indiquée
  void neighCells(const CellType& cell, ArrayCellType& adj_cells) const {
    Traits::neighCells(cell, adj_cells, m_cart_conn);
  }

  //! Récupère les ids locaux des mailles adjacente à une maille identifiée par son itérateur
  // Vaut -1 si la maille voisine n'existe pas à la position indiquée
  void neighCells(const CellEnumeratorType& cell_i, ArrayCellType& adj_cells) const {
    neighCells(*cell_i, adj_cells);
  }

  //! Récupère les ids locaux des mailles adjacente à une maille 
  // Applique des conditions aux limites (Boundary Conditions) si au bord du domaine
  void neighCellsBC(const CellType& cell, ArrayCellType& adj_cells) const {
    Traits::neighCells(cell, adj_cells, m_cart_conn);
    Traits::applyBC(adj_cells, m_trace_mng);
  }

  //! Récupère les ids locaux des mailles adjacente à une maille identifiée par son itérateur
  // Applique des conditions aux limites (Boundary Conditions) si au bord du domaine
  void neighCellsBC(const CellEnumeratorType& cell_i, ArrayCellType& adj_cells) const {
    neighCellsBC(*cell_i, adj_cells);
  }

 public:
  const CartesianInterface::CartesianConnectivity& m_cart_conn;
  ITraceMng* m_trace_mng;  //! Tout ça pour un affichage
};

}

#endif

