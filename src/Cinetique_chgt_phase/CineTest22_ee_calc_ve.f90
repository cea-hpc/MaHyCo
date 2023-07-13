
MODULE M_calc_table_re

CONTAINS 

REAL*8 FUNCTION F_T_Tab_re(rho,ene,Tab_re,Tab_rT) 	               	 
! ..... Calcul de T pour un ve
USE M_CST_TABLE_EE
USE M_CST,         ONLY : CST_shiftEint
USE M_CST,         ONLY : CST_Nval_max

USE M_interpol_EOS, ONLY : locate_eos
USE M_interpol_EOS, ONLY : S_INTERPOLV2
USE M_interpol_EOS, ONLY : S_INTERPOLE_RATF1_eos

USE M_calc_table_rt, ONLY : S_calc_thermo_one_rt

USE M_init_tools,    ONLY : F_R_XT_ZBRENT

IMPLICIT NONE																	      
INTEGER, PARAMETER:: nbmail=1              ! une seule maille à la fois 
INTEGER, PARAMETER:: imail=1               ! une seule maille à la fois 
TYPE(TABLE_RE) :: Tab_re
TYPE(TABLE_RT) :: Tab_rT
REAL*8           :: rho(nbmail),ene(nbmail)   ! densite,energie 

REAL*8     :: valeurs(CST_Nval_max)

! ..... Grandeurs locales 
REAL*8     :: rho_loc(nbmail)          ! rho locale
REAL*8     :: ene_loc(nbmail)          ! energie locale
REAL*8     :: ene_loc2(nbmail)          ! energie locale
REAL*8     :: ene_loclog(nbmail)          ! energie locale

INTEGER  :: irhoi(nbmail)          ! localisation r 
INTEGER  :: jenei(nbmail)          ! localisation e
INTEGER, PARAMETER     :: jdelta=2

REAL*8     :: temp_m(nbmail)          ! energie locale
REAL*8     :: dTdr_tmp(nbmail)          ! energie locale
REAL*8     :: dTde_tmp(nbmail)          ! energie locale

INTEGER  :: ii,jj
INTEGER,PARAMETER  :: iprop=8  !etot

REAL*8, PARAMETER        :: DTMIN  =1.       !ecart estimé absolue sur la température K entre ein et eout ==> si non respecté newtown
REAL*8, PARAMETER        :: DTPCMIN=0.005    !ecart estimé ratio   sur la température CV*DTPCMINentre ein et eout ==> si non respecté newtown

REAL*8                   :: xtmp(1),dxtmp(1)


do ii=1,nbmail													 
  rho_loc(ii)=rho(ii)
  ene_loc2(ii)=ene(ii)
enddo
														 
do ii=1,nbmail													 
  CALL S_INTERPOLE_RATF1_eos (Tab_re%rho(1:Tab_re%nr),Tab_re%cbf_etot(1:Tab_re%nr)&			 
  ,Tab_re%nr,rho_loc(ii),xtmp(1),dxtmp(1))									 
  
  
  if (ene_loc2(ii).lt.xtmp(1)) then 										 
    ene_loc2(ii)=xtmp(1)											 
    ene_loclog(ii)=log(CST_shiftEint)  ! ATTENTION en LOG+Shift e0							 
  else														 
    ene_loclog(ii)=log(ene_loc2(ii)-xtmp(1)+CST_shiftEint)  ! ATTENTION en LOG+Shift e0				 
  endif 													 
enddo

do ii=1,nbmail
  ! -------------------------------------------------------- 
  !localisation table r 											
  call locate_eos(Tab_re%rho(1:Tab_re%nr),Tab_re%nr,rho_loc(ii),irhoi(ii))				
  if (irhoi(ii).lt.2)	irhoi(ii)=2										
  if (irhoi(ii).gt.Tab_re%nr-2) irhoi(ii)=Tab_re%nr-2							
  !localisation table e 											
  call locate_eos(Tab_re%ene(1:Tab_re%ne),Tab_re%ne,ene_loclog(ii),jenei(ii))			
  if (jenei(ii).lt.2)	jenei(ii)=2
  if (jenei(ii).gt.Tab_re%ne-2)   jenei(ii)=Tab_re%ne-2
  ! -------------------------------------------------------- 							
enddo														
														
														
! ..... conversion temperature 											
CALL S_INTERPOLV2(Tab_re%nr,Tab_re%ne,Tab_re%rho(1:Tab_re%nr), &					
Tab_re%ene(1:Tab_re%ne),Tab_re%Temp(1:Tab_re%nr*Tab_re%ne),nbmail ,&				
rho_loc(1:nbmail),ene_loclog(1:nbmail), temp_m(1:nbmail), dTdr_tmp(1:nbmail), & 				
dTde_tmp(1:nbmail), irhoi(1:nbmail), jenei(1:nbmail));								

														
DO ii=1,nbmail													
! ..... GESTION DE L'EXTERAPOLATION HORS GRILLE 									
  if (temp_m(ii).LT.Tab_re%Tmin) then 
    temp_m(ii)=Tab_re%Tmin 								
  endif  
  if (temp_m(ii).GT.Tab_re%Tmax) then
    temp_m(ii)=Tab_re%Tmax							
  endif  

! ..... Il peut avoir des divergences entre vT et ve
! ..... On appel la table en vT
  CALL S_calc_thermo_one_rt(rho=rho_loc(ii:ii),Temp=temp_m(ii:ii),&
                          Tab_rT=Tab_rT,prop='etot',&
			  irho_old=irhoi(ii:ii),&
			  val_out=valeurs)
  ene_loc(ii)=valeurs(iprop) ! etot
  if (abs(ene_loc(ii)-ene_loc2(ii)).gt.TB_re(ii)%Cv_3kb*MAX(temp_m(ii)*DTPCMIN,DTMIN)) then 
! ..... localisation table T
    call locate_eos(Tab_rT%Temp(1:Tab_rT%nT),Tab_rT%nT,temp_m(ii),jj)
    if (jj.lt.1+jdelta)           jj=1+jdelta
    if (jj.gt.Tab_rT%nT-jdelta)   jj=Tab_rT%nT-jdelta
    temp_m(ii)=F_r_XT_ZBRENT(Tab_rT=Tab_rT,X_ref=ene_loc2(ii),V_ref=rho_loc(ii),&
                r_min=Tab_rT%Temp(jj-jdelta),r_max=Tab_rT%Temp(jj+jdelta),&
                iprop=iprop,itype=2,SAFE=1)
  endif
ENDDO														

F_T_Tab_re=temp_m(1)






END FUNCTION F_T_Tab_re

!================================================================== 					             		                      	 
!================================================================== 					             		                      	 
 
SUBROUTINE S_calc_thermo_one_re(rho,ene,Tab_rE,Tab_rT,prop,val_out)
  USE M_CST_TABLE_EE
  USE M_calc_table_rt, ONLY : S_calc_thermo_one_rt

  IMPLICIT NONE 
  INTEGER, PARAMETER:: nbmail=1                  ! une seule maille à la fois 
  INTEGER, PARAMETER:: imail=1                   ! une seule maille à la fois 
  REAL*8              :: rho(nbmail),ene(nbmail)   ! densite,energie interne 
  REAL*8              :: val_out(:)                ! valeur de sortie 
  ! ..... si present calcul de la grandeur et sortie 
  ! ..... la taille de val_out depend de prop 
  ! ..... /!\ la MaJ ne se fait que sur la partie du tableau concernee
  ! ..... prop='irho'    : calcul de la position en rho       -> val_out(1)
  ! ..... prop='jtemp'   : calcul de la position en temp      -> val_out(2)
  ! ..... prop='rho'     : calcul de rho                      -> val_out(3)
  ! ..... prop='Temp'    : calcul de Temp                     -> val_out(4)
  ! ..... prop='Ptot'    : calcul de Ptot, dP/dv|T et dP/dT|v -> val_out(5:7)
  ! ..... prop='etot'    : calcul de etot, de/dv|T et de/dT|v -> val_out(8:10)
  ! ..... prop='stot'    : calcul de stot, ds/dv|T et ds/dT|v -> val_out(11:13)
  ! ..... prop='Gibbs'   : calcul de G, dG/dv|T et dG/dT|v    -> val_out(14:16)
  ! ..... autre          : précedent (dans l'ordre) 
  ! ..... autre          : + dv/dP|T,dv/dT|P,de/dP|T,de/dT|P -> val_out(17:20)
  ! ..... autre          : + ds/dP|T,ds/dT|P,dG/dP|T,dG/dT|P -> val_out(21:24)


  TYPE(TABLE_RT) :: Tab_rT
  TYPE(TABLE_RE) :: Tab_rE
  CHARACTER(LEN=*),  OPTIONAL :: prop 
  REAL*8                        :: Temp(nbmail)  ! Temperature


  Temp(imail)=F_T_Tab_re(rho(imail),ene(imail),Tab_rE,Tab_rT)
  
  
  CALL S_calc_thermo_one_rt(rho=rho(imail),Temp=Temp(imail),Tab_rT=Tab_rT,prop=prop,val_out=val_out)

END SUBROUTINE S_calc_thermo_one_re


SUBROUTINE S_calc_thermo_one_ve(vol,ene,Tab_rE,Tab_rT,prop,val_out)
  USE M_CST_TABLE_EE
  USE M_calc_table_rt, ONLY : S_calc_thermo_one_rt

  IMPLICIT NONE 
  INTEGER, PARAMETER:: nbmail=1                  ! une seule maille à la fois 
  INTEGER, PARAMETER:: imail=1                   ! une seule maille à la fois 
  REAL*8              :: vol(nbmail)               ! volume massique 
  REAL*8              :: rho(nbmail),ene(nbmail)   ! densite,energie interne 
  REAL*8              :: val_out(:)                ! valeur de sortie 
  ! ..... si present calcul de la grandeur et sortie 
  ! ..... la taille de val_out depend de prop 
  ! ..... /!\ la MaJ ne se fait que sur la partie du tableau concernee
  ! ..... prop='irho'    : calcul de la position en rho       -> val_out(1)
  ! ..... prop='jtemp'   : calcul de la position en temp      -> val_out(2)
  ! ..... prop='rho'     : calcul de rho                      -> val_out(3)
  ! ..... prop='Temp'    : calcul de Temp                     -> val_out(4)
  ! ..... prop='Ptot'    : calcul de Ptot, dP/dv|T et dP/dT|v -> val_out(5:7)
  ! ..... prop='etot'    : calcul de etot, de/dv|T et de/dT|v -> val_out(8:10)
  ! ..... prop='stot'    : calcul de stot, ds/dv|T et ds/dT|v -> val_out(11:13)
  ! ..... prop='Gibbs'   : calcul de G, dG/dv|T et dG/dT|v    -> val_out(14:16)
  ! ..... autre          : précedent (dans l'ordre) 
  ! ..... autre          : + dv/dP|T,dv/dT|P,de/dP|T,de/dT|P -> val_out(17:20)
  ! ..... autre          : + ds/dP|T,ds/dT|P,dG/dP|T,dG/dT|P -> val_out(21:24)


  TYPE(TABLE_RT) :: Tab_rT
  TYPE(TABLE_RE) :: Tab_rE
  CHARACTER(LEN=*),  OPTIONAL :: prop 
  REAL*8                        :: Temp(nbmail)  ! Temperature


  rho(imail)=1./vol(imail)
  Temp(imail)=F_T_Tab_re(rho(imail),ene(imail),Tab_rE,Tab_rT)
  

  CALL S_calc_thermo_one_rt(rho=rho(imail),Temp=Temp(imail),Tab_rT=Tab_rT,prop=prop,val_out=val_out)

END SUBROUTINE S_calc_thermo_one_ve


END MODULE M_calc_table_re




