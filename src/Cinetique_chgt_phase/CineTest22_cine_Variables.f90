
!-----------------------------------------------------------------------
! ..... Ce fichier définit les constantes et variables pour le modele
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ..... M_CST_TABLE_CINE
! ..... Constantes et variables pour le modele
!-----------------------------------------------------------------------
MODULE M_CST_TABLE_CINE
	INTEGER,PARAMETER        :: nip_max_cine = 6 ! nb de phases au max 

        INTEGER                  :: nip_cine=-1  
	
	REAL(KIND=8), PARAMETER          :: Trans_max_up=1.0
	REAL(KIND=8), PARAMETER          :: Trans_max_dn=0.001

        TYPE TABLE_CINE_CST
            INTEGER              :: type_cinetique
    REAL(KIND=8)                 :: shift_Gii(nip_max_cine)
    REAL(KIND=8)                 :: B_ij(nip_max_cine,nip_max_cine)
    REAL(KIND=8)                 :: nu_ij(nip_max_cine,nip_max_cine)
    REAL(KIND=8)                 :: trans_orig
	END TYPE TABLE_CINE_CST
	
	TYPE TABLE_CINE
    REAL(KIND=8)                 :: Frac_i(nip_max_cine)=-1.
	END TYPE TABLE_CINE

	TYPE TABLE_CINE_LOC
    REAL(KIND=8)                 :: G_i(nip_max_cine)
    REAL(KIND=8)                 :: dFrac_i(nip_max_cine)
    REAL(KIND=8)                 :: R_ij(nip_max_cine,nip_max_cine)
	END TYPE TABLE_CINE_LOC

	TYPE(TABLE_CINE_CST)            :: NPAR_Cine
	

END MODULE M_CST_TABLE_CINE


