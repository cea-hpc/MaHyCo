#ifndef ACC_ENV_DEFAULT_OPTIONS_H
#define ACC_ENV_DEFAULT_OPTIONS_H

/*! \brief Manière de choisir le device attaché au processus
 */
enum eDeviceAffinity {
  DA_none = 0,      //! Aucune sélection de device
  DA_world_rank,    //! device = world_rank%device_count
  DA_node_rank,     //! device =  node_rank%device_count
  DA_cu_hwloc       //! cherche à placer sur le device le plus proche (si hwloc détecté)
};

/*! \brief Manière de répartir Host/Device sur les sous-domaines
 */
enum eHeterogPartition {
  HP_none = 0,      //! Répartition homogène (tout Host ou bien tout Device)
  HP_heterog1,      //! 1ere version de répartition hétérogène
};

#endif
