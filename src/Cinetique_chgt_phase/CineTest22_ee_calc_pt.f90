
MODULE M_calc_table_PT

CONTAINS 

REAL*8 FUNCTION F_T_TabPT(Pre,Temp,Tab_PT,Tab_rT) 	               	 
! ..... Calcul de T pour un ve
USE M_CST_TABLE_EE
USE M_CST,         ONLY : CST_shiftPtot
USE M_CST,         ONLY : CST_Nval_max
USE M_CST,         ONLY : CST_PT_TOL

USE M_interpol_EOS, ONLY : locate_eos
USE M_interpol_EOS, ONLY : S_INTERPOLV2
USE M_interpol_EOS, ONLY : S_INTERPOLE_RATF1_eos

USE M_calc_table_rt, ONLY : S_calc_thermo_one_rt

USE M_init_tools,    ONLY : F_R_XT_ZBRENT

IMPLICIT NONE																	      
INTEGER, PARAMETER:: nbmail=1              ! une seule maille à la fois 
INTEGER, PARAMETER:: imail=1               ! une seule maille à la fois 
TYPE(TABLE_PT) :: Tab_PT
TYPE(TABLE_RT) :: Tab_rT
REAL*8           :: Pre(nbmail),Temp(nbmail)   ! densite,pression

REAL*8     :: valeurs(CST_Nval_max)

! ..... Grandeurs locales 
REAL*8     :: Temp_loc(nbmail)          ! temperature locale
REAL*8     :: pre_loc(nbmail)          ! pression locale
REAL*8     :: pre_loc2(nbmail)          ! pression locale
REAL*8     :: pre_loclog(nbmail)          ! pression locale

INTEGER  :: iprei(nbmail)         ! localisation P
INTEGER  :: jtempi(nbmail)          ! localisation T
INTEGER, PARAMETER     :: jdelta=2

REAL*8     :: rho_m(nbmail)             ! densite locale
REAL*8     :: drdT_tmp(nbmail)          ! drdT|P locale
REAL*8     :: drdP_tmp(nbmail)          ! drdP|T locale

INTEGER  :: ii,jj
INTEGER,PARAMETER  :: iprop=5  !Ptot

REAL*8, PARAMETER        :: DTMIN  =1.       !ecart estimé absolue sur la température K entre ein et eout ==> si non respecté newtown
REAL*8, PARAMETER        :: DTPCMIN=0.005    !ecart estimé ratio   sur la température CV*DTPCMINentre ein et eout ==> si non respecté newtown

REAL*8                   :: xtmp(1),dxtmp(1)


do ii=1,nbmail													 
  Temp_loc(ii)=Temp(ii)
  pre_loc2(ii)=pre(ii)
enddo
														 
do ii=1,nbmail													 
  if (pre_loc2(ii).lt.Tab_PT%Pmin) then 										 
    pre_loc2(ii)=Tab_PT%Pmin											 
    pre_loclog(ii)=log(CST_shiftPtot)  ! ATTENTION en LOG+CST_shiftPtot						 
  else														 
    pre_loclog(ii)=log(pre_loc2(ii)-Tab_PT%Pmin+CST_shiftPtot)  ! ATTENTION en LOG+CST_shiftPtot		 
  endif 													 
enddo


do ii=1,nbmail
  ! -------------------------------------------------------- 
  !localisation table P 											
  call locate_eos(Tab_PT%Ptot(1:Tab_PT%nP),Tab_PT%nP,pre_loclog(ii),iprei(ii))				
  if (iprei(ii).lt.2)		iprei(ii)=2										
  if (iprei(ii).gt.Tab_PT%nP-2) iprei(ii)=Tab_PT%nP-2							
  !localisation table T 											
  call locate_eos(Tab_PT%Temp(1:Tab_PT%nT),Tab_PT%nT,Temp_loc(ii),jtempi(ii))			
  if (jtempi(ii).lt.2)		   jtempi(ii)=2
  if (jtempi(ii).gt.Tab_PT%nT-2)   jtempi(ii)=Tab_PT%nT-2
  ! -------------------------------------------------------- 							
enddo														
														
														
! ..... conversion densite  											
CALL S_INTERPOLV2(Tab_PT%nP,Tab_PT%nT,Tab_PT%Ptot(1:Tab_PT%nP), &					
Tab_PT%Temp(1:Tab_PT%nT),Tab_PT%rho(1:Tab_PT%nP*Tab_PT%nT),nbmail ,&				
pre_loclog(1:nbmail),Temp_loc(1:nbmail), rho_m(1:nbmail), drdP_tmp(1:nbmail), & 				
drdT_tmp(1:nbmail), iprei(1:nbmail), jtempi(1:nbmail));								


														
DO ii=1,nbmail													
! ..... GESTION DE L'EXTERAPOLATION HORS GRILLE 									
  if (rho_m(ii).LT.Tab_PT%rmin) then 
    rho_m(ii)=Tab_PT%rmin 								
  endif  
  if (rho_m(ii).GT.Tab_PT%rmax) then
    rho_m(ii)=Tab_PT%rmax							
  endif  

! ..... Il peut avoir des divergences entre vT et ve
! ..... On appel la table en vT
  CALL S_calc_thermo_one_rt(rho=rho_m(ii:ii),Temp=Temp_loc(ii:ii),&
                          Tab_rT=Tab_rT,prop='ptot',&
			  val_out=valeurs)
   pre_loc(ii)=valeurs(iprop) ! ptot
  if (abs(pre_loc(ii)-pre_loc2(ii)).gt.CST_PT_TOL) then 
! ..... localisation table r 
    call locate_eos(Tab_rT%rho(1:Tab_rT%nr),Tab_rT%nr,rho_m(ii),jj)
    if (jj.lt.1+jdelta)           jj=1+jdelta
    if (jj.gt.Tab_rT%nr-jdelta)   jj=Tab_rT%nr-jdelta
    rho_m(ii)=F_R_XT_ZBRENT(Tab_rT=Tab_rT,X_ref=pre_loc2(ii),V_ref=Temp_loc(ii),&
                r_min=Tab_rT%rho(jj-jdelta),r_max=Tab_rT%rho(jj+jdelta),&
                iprop=iprop,itype=1,SAFE=1)
  endif
ENDDO														

F_T_Tabpt=rho_m(1)






END FUNCTION F_T_TabPT

!================================================================== 					             		                      	 
!================================================================== 					             		                      	 
 
!================================================================== 					             		                      	 
!================================================================== 					             		                      	 
 
SUBROUTINE S_calc_thermo_one_PT(Pre,Temp,Tab_PT,Tab_rT,prop,val_out)
  USE M_CST_TABLE_EE
  USE M_calc_table_rt, ONLY : S_calc_thermo_one_rt

  IMPLICIT NONE 
  INTEGER, PARAMETER:: nbmail=1               ! une seule maille à la fois 
  INTEGER, PARAMETER:: imail=1               ! une seule maille à la fois 
  REAL*8              :: Pre(nbmail),Temp(nbmail)   ! Pression Temperature 
  REAL*8              :: val_out(:)                 ! valeur de sortie 
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
  TYPE(TABLE_PT) :: Tab_PT
  CHARACTER(LEN=*),  OPTIONAL :: prop 
  REAL*8                        :: rho(nbmail)  ! Masse volumique 

  rho(imail)=F_T_TabPT(Pre=Pre(imail),Temp=Temp(imail),Tab_PT=Tab_PT,Tab_rT=Tab_rT)
  
  CALL S_calc_thermo_one_rt(rho=rho(imail),Temp=Temp(imail),Tab_rT=Tab_rT,prop=prop,val_out=val_out)

END SUBROUTINE S_calc_thermo_one_PT

END MODULE M_calc_table_PT




