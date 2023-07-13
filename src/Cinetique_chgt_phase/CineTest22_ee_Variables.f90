!-----------------------------------------------------------------------
! ..... Ce fichier définit les constantes et variables pour le modele
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ..... M_CST_TABLE_EE
! ..... Constantes et variables pour le modele
!-----------------------------------------------------------------------

MODULE M_CST
  REAL*8, PARAMETER   :: CST_TStd = 300.0
  REAL*8, PARAMETER   :: CST_PStd = 1.e-4
  REAL*8, PARAMETER   :: CST_kb = 1.38062*6.02217e-3 !kJ
  REAL*8, PARAMETER   :: CST_3  =3.   
  
  
  INTEGER,PARAMETER :: CST_Nval_max = 24 ! nb de valeur maximale 
! ..... 1:irho /  2:jtemp/  3:rho/ 4:Temp/  5:Ptot/     
! ..... 6:dP/dv|T/  7:dP/dT|v/  8:etot/     9:de/dv|T/ 10:de/dT|v
! .....11:stot / 12:ds/dv|T/ 13:ds/dT|v/ 14:Gibbs/   15:dG/dv|T/ 16:dG/dT|v/
! .....17:dv/dP|T/ 18:dv/dT|P/ 19:de/dP|T/ 20:de/dT|P/ 21:ds/dP|T/ 22:ds/dT|P/
! .....23:dG/dP|T/ 24:dG/dT|P/


! ..... Parametres de convergence pour l'inversion via F_JWL_X_rt
  INTEGER, PARAMETER   :: ZB_PT_ITMAX=100               ! iteration max 
  REAL*8   , PARAMETER   :: ZB_PT_TOL  =1.e-10            ! tolerance P
  REAL*8   , PARAMETER   :: ZB_PT_EPS  =3.e-10            ! tolerance

! ..... Grille d'interpolation T(ve) en e(vT)
INTEGER, PARAMETER  :: CST_NxTgril   = 2   ! Facteur multiplicatif
REAL*8   , PARAMETER  :: CST_shiftEint = 1. ! shift en E0
REAL*8,PARAMETER      :: CST_deltaemin= 1.e-5  ! delta energie min sur la grille
! ..... Grille d'interpolation v(PT) en P(vT)
INTEGER, PARAMETER  :: CST_NxPgril   =2   ! Facteur multiplicatif
REAL*8   , PARAMETER  :: CST_shiftPtot = 1. ! shift en P0
REAL*8,PARAMETER      :: CST_deltaPmin= 1.e-5  ! delta energie min sur la grille

! ..... Grille d'inversion en PT ee_calc_PT 
REAL*8   , PARAMETER   :: CST_PT_TOL  =1.e-8            ! tolerance P
 

! ..... Valeur d'initialisation en PT 
REAL*8   , PARAMETER   :: CST_INIT_P  =1.e-8            ! si 0+-CST_INIT_P init
REAL*8   , PARAMETER   :: CST_INIT_T  =1.e-8            ! si 0+-CST_INIT_T init
REAL*8   , PARAMETER   :: CST_INIT_P0 =1.e-4            ! valeur P0
REAL*8   , PARAMETER   :: CST_INIT_T0 =300.             ! valeur T0
 

END MODULE M_CST



MODULE M_CST_TABLE_EE
	INTEGER,PARAMETER        :: nip_max_ee = 6 ! nb de phases au max 

	INTEGER,PARAMETER        :: nb_polymax =50! nb de points de polygone max 

	REAL*8,PARAMETER         :: my_HUGE    =1.e38 ! Valeur max 
	!REAL*8                   :: my_HUGE_loc    = 1.e50

! ..... Parametres pour le Newton dans S_SEARCH_NEWTON_ISOPT
	TYPE CST_NEWTON
	    INTEGER              :: ITMAX=100 ! Nb d'iteration max 
	    REAL*8                 :: dP=1.e-12 !1.e-5 !  convergence sur P 
	    REAL*8                 :: dT=1.e-10 !1.e-5 !  convergence sur T 
	    REAL*8                 :: de=1.e-6  ! 1.e-2 ! valeur relative modif au lieu de 1.e-6  ! convergence sur e 
	    REAL*8                 :: dv=1.e-4  ! 1.e-2 ! valeur relative modif au lieu de 1.e-4  ! convergence sur v 
	    REAL*8                 :: TMIN=0.
	    REAL*8                 :: VMIN=0.
	END TYPE CST_NEWTON

	TYPE BORNES_VT
	    REAL*8                 :: TMIN=0.
	    REAL*8                 :: VMIN=0.
	    REAL*8                 :: TMAX=0.
	    REAL*8                 :: VMAX=0.
	END TYPE BORNES_VT
! ..... Utilise dans S_CALC_CINE_PT cette grandeur pour les E/S 
! ..... Allocation n'est pas dynamique pour gagner du temps lors des appels 
	TYPE Thermo_Newton
	    INTEGER              :: iconv
	    REAL*8                 :: Ptot
	    REAL*8                 :: Temp
	    REAL*8                 :: rho
	    REAL*8                 :: ene
	    REAL*8                 :: v_j(nip_max_ee)
	    REAL*8                 :: e_j(nip_max_ee)
	    REAL*8                 :: dedP_T_j(nip_max_ee)
	    REAL*8                 :: dvdP_T_j(nip_max_ee)
	    REAL*8                 :: dedT_P_j(nip_max_ee)
	    REAL*8                 :: dvdT_P_j(nip_max_ee)
	    REAL*8                 :: Gibbs(nip_max_ee)
	END TYPE Thermo_Newton

! ..... Definition de l'état Std 
	TYPE Etat_Std 
	  REAL*8            :: Temp 
	  REAL*8            :: Ptot 
	  REAL*8            :: rho
	  REAL*8            :: M_mol
	END TYPE Etat_Std 
	
! ..... Definition de la table en rho-T
	TYPE TABLE_RT
		INTEGER              :: nr
		INTEGER              :: nT
		REAL*8,    ALLOCATABLE :: rho(:) 
		REAL*8,    ALLOCATABLE :: Temp(:) 
		REAL*8,    ALLOCATABLE :: Ptot(:) 
		REAL*8,    ALLOCATABLE :: etot(:) 
		REAL*8,    ALLOCATABLE :: stot(:) 
		REAL*8,    ALLOCATABLE :: Gibbs(:) 

		REAL*8,    ALLOCATABLE :: cbf_etot(:) 

		TYPE(Etat_Std)       :: Std           
	END TYPE TABLE_RT

! ..... Definition de la table de conversion en T(rho-e) 
	TYPE TABLE_RE
		REAL*8                 :: Tmin
		REAL*8                 :: Tmax
		REAL*8                 :: Cv_3kb

		INTEGER              :: nr
		INTEGER              :: ne
		REAL*8,    ALLOCATABLE :: rho(:) 
		REAL*8,    ALLOCATABLE :: ene(:) 
		REAL*8,    ALLOCATABLE :: Temp(:) 

		REAL*8,    ALLOCATABLE :: cbf_etot(:) 
	END TYPE TABLE_RE

! ..... Definition de la table en de conversion en rho(P-T) 
	TYPE TABLE_PT
		REAL*8                 :: Pmin
		REAL*8                 :: rmin
		REAL*8                 :: rmax

		INTEGER              :: nP
		INTEGER              :: nT
		REAL*8,    ALLOCATABLE :: Ptot(:) 
		REAL*8,    ALLOCATABLE :: Temp(:) 
		REAL*8,    ALLOCATABLE :: rho(:) 

	END TYPE TABLE_PT

	TYPE TAB_APPEL_EE
	    REAL*8          :: rho=-1.
	    REAL*8          :: ene
	    REAL*8          :: Ptot
	    REAL*8          :: Temp=-1.
	    REAL*8          :: v_j(nip_max_ee)=-1
	    REAL*8          :: e_j(nip_max_ee)
	END TYPE TAB_APPEL_EE

	TYPE POLYGONE_EE
 	    INTEGER              :: nb_pts
	    REAL*8                 :: PT(nb_polymax,2) ! grille en PT
	END TYPE POLYGONE_EE



	REAL*8, ALLOCATABLE              :: Frac_init(:)        ! % de phase stable à l'init
	REAL*8,    PARAMETER             :: eps_init=1.e-1 ! phase stable à l'init

	INTEGER                        :: nip_ee=-1    
	TYPE(TABLE_RT),    ALLOCATABLE :: TB_rT(:) 
	TYPE(TABLE_RE),    ALLOCATABLE :: TB_re(:) 
	TYPE(TABLE_PT),    ALLOCATABLE :: TB_PT(:) 
	
	TYPE(CST_NEWTON)               :: NPAR_CONV_NEWTON
	
	TYPE(BORNES_VT),     ALLOCATABLE :: LIMIT_TABLE(:)
	TYPE(POLYGONE_EE),   ALLOCATABLE :: POLYGONE(:)

	INTEGER              :: nb_restart
        REAL*8,   ALLOCATABLE  :: restart_PT(:,:) ! grille en PT

	INTEGER, PARAMETER   :: iFlag_conv=-800
	
	
! ..... Parametres par défaut si solutions non physiques dans couplage_ve
	REAL*8, PARAMETER      :: NPAR_TEMP_MIN=0.1     ! Temp Min 
	REAL*8, PARAMETER      :: NPAR_DPDE_MIN=1.e-5   ! DPDE Min 
	REAL*8, PARAMETER      :: NPAR_CSO2_MIN=1.e-12  ! CSO2 Min 
	
END MODULE M_CST_TABLE_EE


