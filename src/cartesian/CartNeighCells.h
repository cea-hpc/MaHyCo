#ifndef CARTESIAN_CART_NEIGH_CELLS_H
#define CARTESIAN_CART_NEIGH_CELLS_H

#include "cartesian/CartesianNumberingT.h"
#include "cartesian/CartesianGridT.h"
#include "arcane/utils/ArcaneGlobal.h"
#include "arcane/utils/ITraceMng.h"

namespace Cartesian {

// Determine si un decallage en dimension DIM est valide (ie ne deborde pas du domaine de calcul)
template<Integer DIM>
struct ValidDecal {
  static bool isValid(const LocalIdType3& idx3, const LocalIdType3& decal3, const LocalIdType3& ncellsm1) {
    return (ValidDecal<DIM-1>::isValid(idx3, decal3, ncellsm1) &&
        ((decal3[DIM-1]==0) || (decal3[DIM-1]<0 && idx3[DIM-1]>0) || (decal3[DIM-1]>0 && idx3[DIM-1]<ncellsm1[DIM-1]))
        );
  }
};

// Specialisation 1D
template<>
struct ValidDecal<1> {
  static bool isValid(const LocalIdType3& idx3, const LocalIdType3& decal3, const LocalIdType3& ncellsm1) {
    return ((decal3[0]==0) || (decal3[0]<0 && idx3[0]>0) || (decal3[0]>0 && idx3[0]<ncellsm1[0]));
  }
};

/*
 * \brief Caractéristiques pour CartNeighCells en donction de la dimension
 */
template<Integer DIM>
class CartNeighCellsTraits {
 public:
  using eCellPos = Integer;
  static constexpr Integer stencil_sz = -1;
  using DecalType = LocalIdType3[1];
  static constexpr DecalType decal = {{0,0,0}};
};

// Spécialisation 2D
template<>
class CartNeighCellsTraits<2> {
 public:
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
  using DecalType = LocalIdType3[stencil_sz];
  static constexpr DecalType decal = {
    /*SW*/{-1, -1, 0},
    /*S */{ 0, -1, 0},
    /*SE*/{+1, -1, 0},
    /*E */{+1,  0, 0},
    /*NE*/{+1, +1, 0},
    /*N */{ 0, +1, 0},
    /*NW*/{-1, +1, 0},
    /*W */{-1,  0, 0},
    /*C */{ 0,  0, 0}
  };
};

#if 0
// Numérotation en faisant variant X, puis Y, puis Z
struct StencilNumXYZ {
  enum eCellPos {
    // Plan avant, in Front of
    FSW = 0,
    FS  = 1,
    FSE = 2,
    FW  = 3,
    FC  = 4,
    FE  = 5,
    FNW = 6,
    FN  = 7,
    FNE = 8,
    // Plan de la maille Courante
    CSW = 9,
    CS  = 10,
    CSE = 11,
    CW  = 12,
    CC  = 13,
    CE  = 14,
    CNW = 15,
    CN  = 16,
    CNE = 17,
    // Plan arrière, Behind
    BSW = 18,
    BS  = 19,
    BSE = 20,
    BW  = 21,
    BC  = 22,
    BE  = 23,
    BNW = 24,
    BN  = 25,
    BNE = 26,
    //! C = Center, idem que CC, pourquoi ?
    C   = 27,  
    Max = 28,
    Invalid = (-1)
  };
  static constexpr Integer stencil_sz = (Integer)Max;
  using DecalType = LocalIdType3[stencil_sz];

  static constexpr DecalType decal = {
    // Plan avant, in Front of
    /*FSW*/{-1, -1, -1},
    /*FS */{ 0, -1, -1},
    /*FSE*/{+1, -1, -1},
    /*FW */{-1,  0, -1},
    /*FC */{ 0,  0, -1},
    /*FE */{+1,  0, -1},
    /*FNW*/{-1, +1, -1},
    /*FN */{ 0, +1, -1},
    /*FNE*/{+1, +1, -1},
    // Plan de la maille Courante
    /*CSW*/{-1, -1,  0},
    /*CS */{ 0, -1,  0},
    /*CSE*/{+1, -1,  0},
    /*CW */{-1,  0,  0},
    /*CC */{ 0,  0,  0}, // Centre
    /*CE */{+1,  0,  0},
    /*CNW*/{-1, +1,  0},
    /*CN */{ 0, +1,  0},
    /*CNE*/{+1, +1,  0},
    // Plan arrière, Behind
    /*BSW*/{-1, -1, +1},
    /*BS */{ 0, -1, +1},
    /*BSE*/{+1, -1, +1},
    /*BW */{-1,  0, +1},
    /*BC */{ 0,  0, +1},
    /*BE */{+1,  0, +1},
    /*BNW*/{-1, +1, +1},
    /*BN */{ 0, +1, +1},
    /*BNE*/{+1, +1, +1},
    //! C = Center, idem que CC, pourquoi ?
    /*C  */{ 0,  0,  0}  
  };
};
#endif

// Spécialisation 3D
template<>
class CartNeighCellsTraits<3> {
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
  using DecalType = LocalIdType3[stencil_sz];

  static constexpr DecalType decal = {
    // Plan Y = Left
    /*LLL*/{-1, -1, -1},
    /*CLL*/{ 0, -1, -1},
    /*RLL*/{+1, -1, -1},
    /*LLC*/{-1, -1,  0},
    /*CLC*/{ 0, -1,  0},
    /*RLC*/{+1, -1,  0},
    /*LLR*/{-1, -1, +1},
    /*CLR*/{ 0, -1, +1},
    /*RLR*/{+1, -1, +1},
    // Plan Y = Center
    /*LCL*/{-1,  0, -1},
    /*CCL*/{ 0,  0, -1},
    /*RCL*/{+1,  0, -1},
    /*LCC*/{-1,  0,  0},
    /*CCC*/{ 0,  0,  0}, // Centre
    /*RCC*/{+1,  0,  0},
    /*LCR*/{-1,  0, +1},
    /*CCR*/{ 0,  0, +1},
    /*RCR*/{+1,  0, +1},
    // Plan Y = Right
    /*LRL*/{-1, +1, -1}, 
    /*CRL*/{ 0, +1, -1},
    /*RRL*/{+1, +1, -1},
    /*LRC*/{-1, +1,  0},
    /*CRC*/{ 0, +1,  0},
    /*RRC*/{+1, +1,  0},
    /*LRR*/{-1, +1, +1},
    /*CRR*/{ 0, +1, +1},
    /*RRR*/{+1, +1, +1},
    //! C = Center, idem que CC, pourquoi ?
    /*CCC*/{ 0,  0,  0}  
  };
};

// Les conditions aux limites, à spécialiser
template<Integer DIM>
class NeighCellsBC {
 private:
  using ArrayCellType = UniqueArray<LocalIdType>;
 public:
  static void applyBC(ArrayCellType& , ITraceMng* ) {
    ARCANE_ASSERT(false, ("NeighCellsBC<DIM>::applyBC(...) doit être spécialisé"));
  }
};

// Specialisation 2D
template<>
class NeighCellsBC<2> {
 private:
  using ArrayCellType = UniqueArray<LocalIdType>;
  using Traits = CartNeighCellsTraits<2>;
  using CP = Traits::eCellPos;
 public:
  static void applyBC(ArrayCellType& adj_cells, ITraceMng* trace_mng) {

    // Conditions aux limites
    // Bord du bas
    if (ItemId::null(adj_cells[CP::S]))
      adj_cells[CP::SW] = adj_cells[CP::W], adj_cells[CP::S] = adj_cells[CP::C], adj_cells[CP::SE] = adj_cells[CP::E];
    // Bord de droite
    if (ItemId::null(adj_cells[CP::E]))
      adj_cells[CP::SE] = adj_cells[CP::S], adj_cells[CP::E] = adj_cells[CP::C], adj_cells[CP::NE] = adj_cells[CP::N];
    // Bord du haut
    if (ItemId::null(adj_cells[CP::N]))
      adj_cells[CP::NE] = adj_cells[CP::E], adj_cells[CP::N] = adj_cells[CP::C], adj_cells[CP::NW] = adj_cells[CP::W];
    // Bord de gauche
    if (ItemId::null(adj_cells[CP::W]))
      adj_cells[CP::NW] = adj_cells[CP::N], adj_cells[CP::W] = adj_cells[CP::C], adj_cells[CP::SW] = adj_cells[CP::S];
    // Check
    ARCANE_ASSERT( (ItemId::null(adj_cells[CP::SW]) || ItemId::null(adj_cells[CP::S]) || ItemId::null(adj_cells[CP::SE]) || ItemId::null(adj_cells[CP::E]) ||
        ItemId::null(adj_cells[CP::NE]) || ItemId::null(adj_cells[CP::N]) || ItemId::null(adj_cells[CP::NW]) || ItemId::null(adj_cells[CP::W])),
        ("Revoir l'algorithme : on ne devrait pas se retrouver ici.") );
  }
};

// Specialisation 3D
template<>
class NeighCellsBC<3> {
 private:
  using ArrayCellType = UniqueArray<LocalIdType>;
 public:
  static void applyBC(ArrayCellType& adj_cells, ITraceMng* ) {

    // Traitement des conditions aux bords
    // Bord du bas
    if (ItemId::null(adj_cells[4])) {
      for (Integer ii = 0 ; ii < 9 ; ++ii)
        adj_cells[ii] = adj_cells[ii+9];
    }
    // Bord du haut
    if (ItemId::null(adj_cells[22])) {
      for (Integer ii(18) ; ii < 27 ; ++ii) {
        adj_cells[ii] = adj_cells[ii-9];
      }
    }
    // Bord de droite
    if (ItemId::null(adj_cells[14])) {
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
    if (ItemId::null(adj_cells[12])) {
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
    if (ItemId::null(adj_cells[10])) {
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
    if (ItemId::null(adj_cells[16])) {
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
 * Permet de retrouver les mailles voisines d'une maille
 */
/*---------------------------------------------------------------------------*/
template<Integer DIM>
class CartNeighCells {
 public:
  //! Type pour la numérotation cartésienne sur des identifiants locaux
  using CartesianNumbering = CartesianNumberingT<LocalIdType>;

  using CellType = LocalIdType;
  using CellEnumeratorType = CartCellEnumerator;

  //! Type du tableau rempli par neighCells
  using ArrayCellType = UniqueArray<CellType>;

  //! Alias
  using Traits = CartNeighCellsTraits<DIM>;

  //! Positions relatives des mailles adjacentes par rapport à une maille centrale
  using eCellPos = typename Traits::eCellPos;

  //! Taille à laquelle doit être alloué le tableau pour neighCells
  static constexpr Integer stencil_sz = Traits::stencil_sz;

 public:
  CartNeighCells(const CartesianNumbering& cart_numb_cell, ITraceMng* trace_mng) 
  : m_cart_numb_cell (cart_numb_cell), m_trace_mng(trace_mng),
  m_delta3 (m_cart_numb_cell.delta3())
  {
    ARCANE_ASSERT(m_cart_numb_cell.dimension() == DIM, ("Les dimensions sont incoherentes"));
    for(Integer d(0) ; d<DIM ; ++d) {
      m_ncellsm1[d]=m_cart_numb_cell.nbItemDir(d)-1;
    }
    for(Integer d(DIM) ; d<3 ; ++d) {
      m_ncellsm1[d]=0;
    }
    // Pré-calculs des offsets selon les positions des voisins
    for(Integer cp(0) ; cp<stencil_sz ; ++cp) {
      _offset[cp]=Traits::decal[cp][0]*m_delta3[0]+Traits::decal[cp][1]*m_delta3[1]+Traits::decal[cp][2]*m_delta3[2];
    }
  }

  //! Récupère les ids locaux des mailles adjacente à une maille identifiée par son itérateur
  // Vaut -1 si la maille voisine n'existe pas à la position indiquée
  void neighCells(const CellEnumeratorType& cell_i, ArrayCellType& adj_cells) const {
    ARCANE_ASSERT(stencil_sz<=adj_cells.size(), ("adj_cells doit etre au moins de taille stencil_sz"));
    const auto& idx3=cell_i.itemIdx();
    LocalIdType cell_id=cell_i.localId();
    for(Integer cp(0) ; cp<stencil_sz ; ++cp) {
      adj_cells[cp]=(_validCP(cp, idx3) ? cell_id+_offset[cp] : -1);
    }
  }

  //! Récupère les ids locaux des mailles adjacente à une maille identifiée par son itérateur
  // Applique des conditions aux limites (Boundary Conditions) si au bord du domaine
  void neighCellsBC(const CellEnumeratorType& cell_i, ArrayCellType& adj_cells) const {
    neighCells(cell_i, adj_cells);

    NeighCellsBC<DIM>::applyBC(adj_cells, m_trace_mng);
  }

 protected:

  inline bool _validCP(Integer cp, const LocalIdType3& idx3) const {
    return ValidDecal<DIM>::isValid(idx3, Traits::decal[cp], m_ncellsm1);
  }

 protected:
  // cell_adj_id(pos) = cell_id + _offset[pos]
  LocalIdType _offset[stencil_sz];

 public:
  const CartesianNumbering& m_cart_numb_cell;  //! Numérotation cartésienne sur les mailles
  ITraceMng* m_trace_mng;  //! Tout ça pour un affichage
  LocalIdType3 m_ncellsm1; //! Nb de mailles -1 dans chaque direction
  const LocalIdType3& m_delta3; //! +-delta a appliquer sur m_cell_id pour passer a la maille suivante/precedente pour une direction donnée
};

}

#endif

