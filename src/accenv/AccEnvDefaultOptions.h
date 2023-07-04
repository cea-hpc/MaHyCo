#ifndef ACC_ENV_DEFAULT_OPTIONS_H
#define ACC_ENV_DEFAULT_OPTIONS_H

/*! \brief Manière de répartir Host/Device sur les sous-domaines
 */
enum eHeterogPartition {
  HP_none = 0,      //! Répartition homogène (tout Host ou bien tout Device)
  HP_heterog1,      //! 1ere version de répartition hétérogène
};

#endif
