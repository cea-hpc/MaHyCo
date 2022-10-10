#ifndef ACC_ENV_MULTI_ENV_UTILS_H
#define ACC_ENV_MULTI_ENV_UTILS_H

#include "accenv/AcceleratorUtils.h"
#include "accenv/BufAddrMng.h"

#include "arcane/MeshVariableScalarRef.h"
#include "arcane/MeshVariableArrayRef.h"
#include "arcane/accelerator/VariableViews.h"
#include <arcane/IMesh.h>
#include <arcane/VariableBuildInfo.h>
#include <arcane/utils/IMemoryRessourceMng.h>

#include <arcane/materials/ComponentPartItemVectorView.h>
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/MeshMaterialVariableRef.h"

/*---------------------------------------------------------------------------*/
/* Pour créer une vue sur les valeurs d'un environnement                     */
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::Materials;

template<typename value_type>
ArrayView<value_type> envView(CellMaterialVariableScalarRef<value_type>& var_menv, IMeshEnvironment* env) {
  // Pour rappel, [0] fait référence à .globalVariable() d'où le +1
  return var_menv._internalValue()[env->id()+1];
}

/*---------------------------------------------------------------------------*/
/* Pour repérer un EnvCell (fortement inspiré de MatVarIndex)                */
/*---------------------------------------------------------------------------*/
class EnvVarIndex {
 public:
  EnvVarIndex() :
    m_array_index (-1),
    m_value_index (-1)
  {  }

  ARCCORE_HOST_DEVICE EnvVarIndex(Int32 array_index, Int32 value_index) :
    m_array_index (array_index),
    m_value_index (value_index)
  {
  }
  ARCCORE_HOST_DEVICE Int32 arrayIndex() const { return m_array_index; }
  ARCCORE_HOST_DEVICE Int32 valueIndex() const { return m_value_index; }
 protected:
  Int32 m_array_index;
  Int32 m_value_index;
};

/*---------------------------------------------------------------------------*/
/* Vues sur les variables multi-environnement                                */
/*---------------------------------------------------------------------------*/
template<typename value_type>
class MultiEnvView {
 public:
  MultiEnvView(const Span< Span<value_type> >& spn) :
    m_var_menv_views (spn)
  {
  }

  ARCCORE_HOST_DEVICE MultiEnvView(const MultiEnvView<value_type>& rhs) :
    m_var_menv_views (rhs.m_var_menv_views)
  {
  }

  ARCCORE_HOST_DEVICE value_type operator[](const EnvVarIndex& evi) const {
    return m_var_menv_views[evi.arrayIndex()][evi.valueIndex()];
  }

  ARCCORE_HOST_DEVICE void setValue(const EnvVarIndex& evi, value_type val) const {
    m_var_menv_views[evi.arrayIndex()].setItem(evi.valueIndex(), val);
  }

  ARCCORE_HOST_DEVICE value_type& ref(const EnvVarIndex& evi) const {
    return m_var_menv_views[evi.arrayIndex()].item(evi.valueIndex());
  }

 public:
  Span< Span<value_type> > m_var_menv_views;
};

/*---------------------------------------------------------------------------*/
/* Pour créer des vues sur les variables multi-environnement                 */
/*---------------------------------------------------------------------------*/
template<typename value_type>
class MultiEnvVar {
 public:
  MultiEnvVar(CellMaterialVariableScalarRef<value_type>& var_menv, IMeshMaterialMng* mm) :
   m_var_menv_impl(platform::getAcceleratorHostMemoryAllocator(), mm->environments().size()+1)
  {
    m_var_menv_impl[0].setArray(Span<value_type>(var_menv._internalValue()[0]));
    ENUMERATE_ENV(ienv, mm) {
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();
      m_var_menv_impl[env_id+1].setArray(Span<value_type>(envView(var_menv,env)));
    }
  }

  //! Pour y accéder en lecture/écriture
  auto span() {
    return MultiEnvView<value_type>(m_var_menv_impl);
  }

 protected:
  UniqueArray< Span<value_type> > m_var_menv_impl;
};

/*---------------------------------------------------------------------------*/
/* Accès aux données d'une variable multi-env sans protection                */
/*---------------------------------------------------------------------------*/
template<typename value_type>
class MultiEnvData {
 public:
  MultiEnvData(value_type** dta) :
    m_var_menv_data (dta)
  {
  }

  ARCCORE_HOST_DEVICE MultiEnvData(const MultiEnvData<value_type>& rhs) :
    m_var_menv_data (rhs.m_var_menv_data)
  {
  }

  ARCCORE_HOST_DEVICE value_type operator[](const EnvVarIndex& evi) const {
    return m_var_menv_data[evi.arrayIndex()][evi.valueIndex()];
  }

  ARCCORE_HOST_DEVICE void setValue(const EnvVarIndex& evi, value_type val) const {
    m_var_menv_data[evi.arrayIndex()][evi.valueIndex()] = val;
  }

  ARCCORE_HOST_DEVICE value_type& ref(const EnvVarIndex& evi) const {
    return m_var_menv_data[evi.arrayIndex()][evi.valueIndex()];
  }

  auto data() {
    return m_var_menv_data;
  }

 public:
  value_type** m_var_menv_data; //!< Les données non gardées en 2 dimensions
};

/*---------------------------------------------------------------------------*/
/* Pour créer des vues non protégées sur une variable multi-environnement    */
/*---------------------------------------------------------------------------*/
template<typename value_type>
class MultiEnvDataVar {
 public:
  MultiEnvDataVar(CellMaterialVariableScalarRef<value_type>& var_menv, IMeshMaterialMng* mm,
      eMemoryRessource mem_res = eMemoryRessource::UnifiedMemory) :
   m_buf_addr(platform::getDataMemoryRessourceMng()->getAllocator(mem_res), 
       mm->environments().size()+1)
  {
    m_var_menv_impl = _viewFromBuf(m_buf_addr);

    _initViewFromVar(var_menv, mm);
  }

  // Allocation mais pas encore d'association à une variable multi-env
  MultiEnvDataVar(IMeshMaterialMng* mm, eMemoryRessource mem_res) :
   m_buf_addr(platform::getDataMemoryRessourceMng()->getAllocator(mem_res), 
       mm->environments().size()+1)
  {
    m_var_menv_impl = _viewFromBuf(m_buf_addr);
  }

  // A partir d'un buffer déjà existant
  MultiEnvDataVar(CellMaterialVariableScalarRef<value_type>& var_menv, IMeshMaterialMng* mm,
      ArrayView<Int64> buf_addr) 
  {
    m_var_menv_impl = _viewFromBuf(buf_addr);

    _initViewFromVar(var_menv, mm);
  }

  MultiEnvDataVar(ArrayView<Int64> buf_addr) 
  {
    m_var_menv_impl = _viewFromBuf(buf_addr);
  }

  // Copie asynchrone de h_var_menv (sur hôte)
  void asyncCpy(MultiEnvDataVar<value_type>& h_var_menv, Ref<RunQueue> queue) {
    value_type** h_menv_var_vw = h_var_menv.span().data(); // Vue sur l'hôte
    void* d_dst     = static_cast<void*>(m_var_menv_impl.data()); // adresse sur device
    void* h_src     = static_cast<void*>(h_menv_var_vw);
    Int64 sz = m_var_menv_impl.size()*sizeof(value_type*);
    queue->copyMemory(ax::MemoryCopyArgs(d_dst, h_src, sz).addAsync());
  }

  //! Pour y accéder en lecture/écriture
  auto span() {
    return MultiEnvData<value_type>(m_var_menv_impl.data());
  }

 protected:

  ArrayView<value_type*> _viewFromBuf(ArrayView<Int64> buf_addr) {
    void* base_ptr_v = static_cast<void*>(buf_addr.data());
    value_type** base_ptr = static_cast<value_type**>(base_ptr_v);
    return ArrayView<value_type*>(buf_addr.size(), base_ptr);
  }

  void _initViewFromVar(CellMaterialVariableScalarRef<value_type>& var_menv, 
      IMeshMaterialMng* mm) {

    m_var_menv_impl[0] = var_menv._internalValue()[0].data();
    ENUMERATE_ENV(ienv, mm) {
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();
      m_var_menv_impl[env_id+1] = envView(var_menv,env).data();
    }
  }

 protected:
  UniqueArray<Int64> m_buf_addr;  //<! Int64 va être converti en value_type*
  ArrayView< value_type* > m_var_menv_impl;
};

/*---------------------------------------------------------------------------*/
/* Pour créer des vues non protégées sur une variable multi-environnement en */
/* mémoires Host et Device                                                   */
/*---------------------------------------------------------------------------*/
template<typename value_type>
class MultiEnvVarHD {
 public:
  MultiEnvVarHD(CellMaterialVariableScalarRef<value_type>& var_menv,
      BufAddrMng* bam) :
    m_menv_var_h(var_menv, bam->materialMng(), bam->nextHostView()),
    m_menv_var_d(bam->nextDeviceView())
  {
  }

  //! Vue en mémoire Hôte
  auto spanH() {
    return m_menv_var_h.span();
  }

  //! Vue en mémoire Device
  auto spanD() {
    return m_menv_var_d.span();
  }

 protected:
  MultiEnvDataVar<value_type> m_menv_var_h;  //! View in HOST memory on multi-mat data
  MultiEnvDataVar<value_type> m_menv_var_d;  //! View in DEVICE memory on multi-mat data
};

/*---------------------------------------------------------------------------*/
/* Vue sur stockage du multi-env                                             */
/*---------------------------------------------------------------------------*/
class MultiEnvCellViewIn {
 public:
  MultiEnvCellViewIn(ax::RunCommand& command, Integer max_nb_env,
      const VariableCellInteger& v_nb_env,
      const UniqueArray<Int16>& v_l_env_arrays_idx,
      const VariableCellArrayInteger& v_l_env_values_idx,
      const VariableCellInteger& v_env_id) :
    m_max_nb_env (max_nb_env),
    m_nb_env_in (ax::viewIn(command,v_nb_env)),
    m_l_env_arrays_idx_in (v_l_env_arrays_idx.constSpan()),
    m_l_env_values_idx_in (ax::viewIn(command,v_l_env_values_idx)),
    m_env_id_in (ax::viewIn(command,v_env_id))
  {}

  ARCCORE_HOST_DEVICE MultiEnvCellViewIn(const MultiEnvCellViewIn& rhs) : 
    m_max_nb_env (rhs.m_max_nb_env),
    m_nb_env_in (rhs.m_nb_env_in),
    m_l_env_arrays_idx_in (rhs.m_l_env_arrays_idx_in),
    m_l_env_values_idx_in (rhs.m_l_env_values_idx_in),
    m_env_id_in (rhs.m_env_id_in)
  {}

  ARCCORE_HOST_DEVICE Integer nbEnv(CellLocalId cid) const {
    return m_nb_env_in[cid];
  }

  ARCCORE_HOST_DEVICE EnvVarIndex envCell(CellLocalId cid, Integer ienv) const {
    return EnvVarIndex(
        m_l_env_arrays_idx_in[cid.localId()*m_max_nb_env+ienv],
        m_l_env_values_idx_in[cid][ienv]);
  }

  ARCCORE_HOST_DEVICE Integer envId(CellLocalId cid, Integer ienv) const {
    return (m_env_id_in[cid]>=0 ? 
        m_env_id_in[cid] : 
        m_l_env_arrays_idx_in[cid.localId()*m_max_nb_env+ienv]-1);
  }

 protected:
  Integer m_max_nb_env;
  ax::VariableCellInt32InView m_nb_env_in;
  Span<const Int16> m_l_env_arrays_idx_in;
  ax::ItemVariableArrayInViewT<Cell,Int32> m_l_env_values_idx_in;
  ax::VariableCellInt32InView m_env_id_in;
};

/*---------------------------------------------------------------------------*/
/* Stockage du multi-env                                                     */
/*---------------------------------------------------------------------------*/
class MultiEnvCellStorage {
 public:
  MultiEnvCellStorage(IMeshMaterialMng* mm, AccMemAdviser* acc_mem_adv) :
    m_mesh_material_mng (mm),
    m_max_nb_env (mm->environments().size()),
    m_nb_env(VariableBuildInfo(mm->mesh(), "NbEnv" , IVariable::PNoDump| IVariable::PNoNeedSync)),
    m_l_env_arrays_idx(platform::getAcceleratorHostMemoryAllocator()),
    m_l_env_values_idx(VariableBuildInfo(mm->mesh(), "LEnvValuesIdx" , IVariable::PNoDump| IVariable::PNoNeedSync | IVariable::PSubDomainDepend)),
    m_env_id(VariableBuildInfo(mm->mesh(), "EnvId" , IVariable::PNoDump| IVariable::PNoNeedSync))
  {
    m_l_env_arrays_idx.resize(m_max_nb_env*mm->mesh()->allCells().size());
    acc_mem_adv->setReadMostly(m_l_env_arrays_idx.view());
    m_l_env_values_idx.resize(m_max_nb_env);
  }

  //! Remplissage
  void buildStorage(ax::Runner& runner, Materials::MaterialVariableCellInteger& v_global_cell) {
    PROF_ACC_BEGIN(__FUNCTION__);

    auto queue = makeQueue(runner);
    {
      auto command = makeCommand(queue);

      auto inout_nb_env = ax::viewInOut(command, m_nb_env);

      command << RUNCOMMAND_ENUMERATE(Cell, cid, m_mesh_material_mng->mesh()->allCells()){
        // Init du nb d'env par maille qui va être calculé à la boucle suivante
        inout_nb_env[cid] = 0;
      };
    }

    Integer max_nb_env = m_max_nb_env; // on ne peut pas utiliser un attribut dans le kernel
    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();

      // Mailles mixtes
      {
        auto command = makeCommand(queue);

        auto inout_nb_env = ax::viewInOut(command, m_nb_env);
        auto out_l_env_values_idx = ax::viewOut(command, m_l_env_values_idx);
        auto out_l_env_arrays_idx = m_l_env_arrays_idx.span();

        Span<const Integer> in_global_cell(envView(v_global_cell, env));

	Span<const Int32> in_imp_idx(env->impureEnvItems().valueIndexes());
	Integer nb_imp = in_imp_idx.size();

        command << RUNCOMMAND_LOOP1(iter, nb_imp) {
	  auto imix = in_imp_idx[iter()[0]]; // iter()[0] \in [0,nb_imp[
          CellLocalId cid(in_global_cell[imix]); // on récupère l'identifiant de la maille globale

          Integer index_cell = inout_nb_env[cid];

          // On relève le numéro de l'environnement 
          // et l'indice de la maille dans la liste de mailles mixtes env
          out_l_env_arrays_idx[cid*max_nb_env+index_cell] = env_id+1; // décalage +1 car 0 est pris pour global
          out_l_env_values_idx[cid][index_cell] = imix;

          inout_nb_env[cid] = index_cell+1; // ++ n'est pas supporté
        };
      }

      // Mailles pures
      {
        auto command = makeCommand(queue);

        const auto& pure_env_items = env->pureEnvItems();
        // Pour les mailles pures, valueIndexes() est la liste des ids locaux des mailles
        Span<const Int32> in_cell_id(pure_env_items.valueIndexes());

        auto inout_nb_env = ax::viewInOut(command, m_nb_env);
        auto out_l_env_values_idx = ax::viewOut(command, m_l_env_values_idx);
        auto out_l_env_arrays_idx = m_l_env_arrays_idx.span();

        // Nombre de mailles pures de l'environnement
        Integer nb_pur = pure_env_items.nbItem();

        command << RUNCOMMAND_LOOP1(iter, nb_pur) {
          auto [ipur] = iter(); // ipur \in [0,nb_pur[
          CellLocalId cid(in_cell_id[ipur]); // accés indirect à la valeur de la maille

          Integer index_cell = inout_nb_env[cid];
#ifndef ARCCORE_DEVICE_CODE
          ARCANE_ASSERT(index_cell==0, ("Maille pure mais index_cell!=0"));
#endif
          out_l_env_arrays_idx[cid*max_nb_env+index_cell] = 0; // 0 référence le tableau global
          out_l_env_values_idx[cid][index_cell] = cid.localId();

          // Equivalent à affecter 1
          inout_nb_env[cid] = index_cell+1; // ++ n'est pas supporté
        };
      }
    }

    checkStorage(v_global_cell);
    PROF_ACC_END;
  }

  //! Verification
  void checkStorage([[maybe_unused]] Materials::MaterialVariableCellInteger& v_global_cell) {
#ifdef ARCANE_DEBUG
    PROF_ACC_BEGIN(__FUNCTION__);
    // m_env_id doit être calculé
    MultiEnvVar<Integer> menv_global_cell(v_global_cell, m_mesh_material_mng);
    auto in_menv_global_cell(menv_global_cell.span());

    ENUMERATE_CELL(icell, m_mesh_material_mng->mesh()->allCells()){
      Cell cell = * icell;
      Integer cid = cell.localId();

      Integer nb_env = m_nb_env[icell];
      Integer nb_env_bis = (m_env_id[cell]>=0 ? 1 : -m_env_id[cell]-1);
      ARCANE_ASSERT(nb_env_bis==nb_env, ("nb_env_bis!=nb_env"));

      for(Integer ienv=0 ; ienv<nb_env ; ++ienv) {
        Integer i=m_l_env_arrays_idx[cid*m_max_nb_env+ienv];
        Integer j=m_l_env_values_idx[icell][ienv];
        EnvVarIndex evi(i,j);

        Integer cid_bis=in_menv_global_cell[evi];
        ARCANE_ASSERT(cid_bis==cid, ("cid_bis!=cid"));
      }
    }
    PROF_ACC_END;
#endif
  }

  //! Vue sur le stockage pour utilisation en lecture sur GPU
  MultiEnvCellViewIn viewIn(ax::RunCommand& command) {
    return MultiEnvCellViewIn(command, m_max_nb_env, m_nb_env, m_l_env_arrays_idx, m_l_env_values_idx, m_env_id);
  }

 protected:
  IMeshMaterialMng* m_mesh_material_mng=nullptr;
  Integer m_max_nb_env;
  VariableCellInteger m_nb_env;  //! Nb d'env par maille
  UniqueArray<Int16> m_l_env_arrays_idx; //! liste des indexes des env par maille
  VariableCellArrayInteger m_l_env_values_idx;
  VariableCellInteger m_env_id;
};

#endif

