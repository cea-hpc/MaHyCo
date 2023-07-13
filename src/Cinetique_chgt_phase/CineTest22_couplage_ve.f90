MODULE SCALC
CONTAINS 

SUBROUTINE S_CALC_CINE_VE(nb,dtime_in,rho,ene,&
							 Pres,Temp,&
                             frac_1,frac_2,frac_3,frac_4,frac_5,frac_6,&
			                 dpde,cs2,conv)

USE M_CST_TABLE_EE,       ONLY : TAB_APPEL_EE  ! type pour les propriétés thermo d'entree
USE M_CST_TABLE_EE,       ONLY : Thermo_Newton ! type pour les propriétés thermo de sortie 
USE M_CST_TABLE_EE,       ONLY : nip_max_ee
USE M_CST_TABLE_EE,       ONLY : Frac_init
USE M_CST_TABLE_EE,       ONLY : eps_init

USE M_CST_TABLE_EE,       ONLY : NPAR_TEMP_MIN
USE M_CST_TABLE_EE,       ONLY : NPAR_DPDE_MIN
USE M_CST_TABLE_EE,       ONLY : NPAR_CSO2_MIN

USE M_CST_TABLE_CINE,     ONLY : NPAR_Cine

USE M_COUPLAGE_TOOLS,     ONLY : S_CALC_CINE_VE_ONE
USE M_COUPLAGE_TOOLS,     ONLY : S_COPIE_THERMO_IN
USE M_COUPLAGE_TOOLS,     ONLY : S_COPIE_THERMO_OUT

USE M_CST_TABLE_EE,       ONLY : nip_ee
USE M_CST_TABLE_EE,       ONLY : iFlag_conv


IMPLICIT NONE 
!.... INPUT:
INTEGER   :: nb
REAL*8     :: dtime_in(nb)
REAL*8      :: rho(nb)
REAL*8      :: ene(nb)


!.... INOUT:
REAL*8     :: Pres(nb)
REAL*8     :: Temp(nb)
REAL*8     :: frac_1(nb)
REAL*8     :: frac_2(nb)
REAL*8     :: frac_3(nb)
REAL*8     :: frac_4(nb)
REAL*8     :: frac_5(nb)
REAL*8     :: frac_6(nb)

REAL*8 	 :: dpde(nb)
REAL*8 	 :: cs2(nb)
REAL*8 	 :: conv(nb)



!.... OUTPUT:
REAL*8 	 :: Pres_o(nb)
REAL*8 	 :: Temp_o(nb)

!.... LOCAL:
INTEGER                          :: ii,jj,kk
INTEGER                          :: niter

REAL*8 	 :: de_dP_T(nb)
REAL*8 	 :: dv_dP_T(nb)
REAL*8 	 :: de_dT_P(nb)
REAL*8 	 :: dv_dT_P(nb)

REAL*8 	 :: dP_dT_v(nb)
REAL*8 	 :: de_dT_v(nb)
REAL*8 	 :: dP_de_v(nb)

REAL*8 	 :: dP_dT_e(nb)
REAL*8 	 :: dv_dT_e(nb)
REAL*8 	 :: dP_dv_e(nb)

REAL*8 	 :: dP_dv_s(nb)

REAL*8	 :: dpde_o(nb)
REAL*8	 :: cs2_o(nb)
REAL*8	 :: conv_o(nb)


TYPE(TAB_APPEL_EE)               :: Thermo_in(nb)
REAL*8                             :: Frac_in(nb,nip_max_ee) ! doit prendre la valeur max , ne prend pas de vecteur

TYPE(Thermo_Newton)              :: Thermo_ou_loc(nb)
REAL*8                             :: Frac_ou_loc(nb,nip_max_ee) ! doit prendre la valeur max , ne prend pas de vecteur

TYPE(Thermo_Newton)              :: Thermo_ou(nb)
REAL*8                             :: Frac_ou(nb,nip_max_ee) ! doit prendre la valeur max , ne prend pas de vecteur


REAL*8                             :: conv_in(nb)      



TYPE(TAB_APPEL_EE)               :: Thermo_in_sav(nb)
REAL*8                           :: Frac_in_sav(nb,nip_max_ee) 

INTEGER                          :: iconverge
REAL*8 	                         :: taux_trans            ! taux de transition maximal

INTEGER,PARAMETER                :: irelaunch=10             ! relance max sinon stop
REAL*8,PARAMETER                 :: taux_trans_ratio=0.9    ! diminution du sous cyclage par ce ratio à chaque non converggence 

! INITIALISATION DU CALCUL
! Affectation des fraction des mailles
Frac_in(1:nb,1)=frac_1(1:nb)
Frac_in(1:nb,2)=frac_2(1:nb)
Frac_in(1:nb,3)=frac_3(1:nb)
if (nip_ee.gt.3) then
Frac_in(1:nb,4)=frac_4(1:nb)
Frac_in(1:nb,5)=frac_5(1:nb)
Frac_in(1:nb,6)=frac_6(1:nb)
endif
Frac_ou_loc(1:nb,1)=0.
Frac_ou_loc(1:nb,2)=0.
Frac_ou_loc(1:nb,3)=0.
if (nip_ee.gt.3) then
Frac_ou_loc(1:nb,4)=0.
Frac_ou_loc(1:nb,5)=0.
Frac_ou_loc(1:nb,6)=0.
endif
Frac_ou(1:nb,1)=0.
Frac_ou(1:nb,2)=0.
Frac_ou(1:nb,3)=0.
if (nip_ee.gt.3) then
Frac_ou(1:nb,4)=0.
Frac_ou(1:nb,5)=0.
Frac_ou(1:nb,6)=0.
endif

! .... Reussi pas à mettre une seule maille avec EOS
taux_trans=NPAR_Cine%trans_orig
! Verification de l'état initial 
! Si la somme des fractions est <eps : initialisation 
do ii=1,nb
  if (SUM(Frac_in(ii,1:nip_ee)).lt.eps_init) then
    Frac_in(ii,1:nip_ee)=Frac_init(1:nip_ee)
  endif
enddo

conv_in(:)=1. ! par defaut tout sera calculé


!write(*,*) 'A l entrée de l fonction'
!write(*,*) '  ',rho(10)
!write(*,*) '  ',ene(10)
!write(*,*) '  ',Temp(10)
!write(*,*) '  ',Pres(10)

! Affectation des grandeurs thermo
Thermo_in(1:nb)%rho=rho(1:nb)
Thermo_in(1:nb)%ene=ene(1:nb)

! /!\ Attention aux inits a verifier si mieux si on force à l'état init ou pas via le Flag T<0 et v_j(:)<0 ?? 
! /!\ on peut passer les Pn-1 et Tn-1 au besoin 
Thermo_in(1:nb)%Ptot=Pres(1:nb)
Thermo_in(1:nb)%Temp=Temp(1:nb)
do jj=1,nip_ee
  Thermo_in(1:nb)%v_j(jj)=-1.
  Thermo_in(1:nb)%e_j(jj)=-1.
enddo



Thermo_ou(1:nb)%iconv=777

! on garde les infos d'origine 
CALL S_COPIE_THERMO_IN(nb,nip_ee,Thermo_in,Frac_in,Thermo_in_sav,Frac_in_sav)


do ii=1,irelaunch
  !write(*,*) " ********************* LANCEMENT ", ii, " **************************"
  
  if (ii.ne.1) then
    CALL S_COPIE_THERMO_IN(nb,nip_ee,Thermo_in_sav,Frac_in_sav,Thermo_in,Frac_in)
  endif 
  ! Bloc de calcul cinetique 
  CALL S_CALC_CINE_VE_ONE(nb,nip_ee,dtime_in,taux_trans,&
			  Thermo_in(1:nb),&
			  Frac_in(1:nb,1:nip_ee),&
			  conv_in(1:nb),&
			  Thermo_ou_loc(1:nb),&
			  Frac_ou_loc(1:nb,1:nip_ee),&
			  conv_o(1:nb))
  
  ! on ne copie que les grandeurs convergées FLAG conv_o(1:nb)=0
  CALL S_COPIE_THERMO_OUT(nb,nip_ee,conv_o(1:nb),&
                          Thermo_ou_loc(1:nb),Frac_ou_loc(1:nb,1:nip_ee),&
                          Thermo_ou(1:nb),Frac_ou(1:nb,1:nip_ee)&
	  		  )
  conv_in(1:nb)=conv_o(1:nb)
  do kk=1,nb
    ! on a convergé avant et on n'a rien recalculé dans S_CALC_CINE_VE_ONE 
    if (int(conv_in(kk)).eq.iFlag_conv) then 
      conv_in(kk)=0.
      Thermo_ou(kk)%iconv=0
    endif  
  enddo
  iconverge=0
  do kk=1,nb
    if (Thermo_ou(kk)%iconv.ne.0) iconverge=-1 
  enddo
  if (iconverge.eq.0) exit
  taux_trans=taux_trans*taux_trans_ratio
enddo
if (iconverge.ne.0) then
    do kk=1,nb
	if (Thermo_ou(kk)%iconv.ne.0) then
	    write(*,'(a,2i4,20e14.6)') 'kk,FLAG,rho_in,ene_in, ',kk,Thermo_ou(kk)%iconv,Thermo_in(kk)%rho,Thermo_in(kk)%ene
	    write(*,*) ' Ptot ' , Thermo_in(kk)%Ptot
	    write(*,*) ' temp ' , Thermo_in(kk)%Temp
	    write(*,*) ' frac out ', Frac_ou(kk,1:nip_ee)
	    write(*,*) ' frac in ', Frac_in(kk,1:nip_ee)
            !conv_o(kk)=1.*Thermo_ou(kk)%iconv
            ! on garde les resultats à n
            Thermo_ou(kk)%Ptot = Thermo_in(kk)%Ptot
            Thermo_ou(kk)%Temp = Thermo_in(kk)%Temp
            Thermo_ou(kk)%dedP_T_j(1:nip_ee) = Thermo_ou_loc(kk)%dedP_T_j(1:nip_ee)
            Thermo_ou(kk)%dvdP_T_j(1:nip_ee) = Thermo_ou_loc(kk)%dvdP_T_j(1:nip_ee)
            Thermo_ou(kk)%dedT_P_j(1:nip_ee) = Thermo_ou_loc(kk)%dedT_P_j(1:nip_ee)
            Thermo_ou(kk)%dvdT_P_j(1:nip_ee) = Thermo_ou_loc(kk)%dvdT_P_j(1:nip_ee)
            Frac_ou(kk,1:nip_ee) = Frac_in(kk,1:nip_ee)
            conv_o(kk)=0.
            write(*,*) ' Problème de convergence maille : ', kk
	endif
    enddo
!STOP
else
conv_o(1:nb)=0.
endif


! ..... Calcul des dérivées. Certaines sont inutiles 
do kk=1,nb
	de_dP_T(kk)=0.
	dv_dP_T(kk)=0.
	de_dT_P(kk)=0.
	dv_dT_P(kk)=0.
	do ii=1,nip_ee
	  de_dP_T(kk)=de_dP_T(kk)+Thermo_ou(kk)%dedP_T_j(ii)*Frac_ou(kk,ii)
	  dv_dP_T(kk)=dv_dP_T(kk)+Thermo_ou(kk)%dvdP_T_j(ii)*Frac_ou(kk,ii) ! ~ 1/KT
	  de_dT_P(kk)=de_dT_P(kk)+Thermo_ou(kk)%dedT_P_j(ii)*Frac_ou(kk,ii) ! ~ Cp
	  dv_dT_P(kk)=dv_dT_P(kk)+Thermo_ou(kk)%dvdT_P_j(ii)*Frac_ou(kk,ii) ! ~ Beta 
	enddo
! ..... Calculs des derivées a v constant
	dP_dT_v(kk)=-dv_dT_P(kk)/dv_dP_T(kk)               ! ~ Beta*Kt ou dS/dV|T
	de_dT_v(kk)=dP_dT_v(kk)*de_dP_T(kk)+de_dT_P(kk)    ! ~ Cv
	dP_de_v(kk)=dP_dT_v(kk)/de_dT_v(kk)                
! ..... Calculs des derivées a e constant
	dP_dT_e(kk)=-de_dT_P(kk)/de_dP_T(kk)  
	dv_dT_e(kk)=dP_dT_e(kk)*dv_dP_T(kk)+dv_dT_P(kk)    
	dP_dv_e(kk)=dP_dT_e(kk)/dv_dT_e(kk)
! ..... Calculs de -Ks/rho|s 
	dP_dv_s(kk)=dP_dv_e(kk)-Thermo_ou(kk)%Ptot*dP_de_v(kk)

	! Allocation des variables de sortie  
	Pres_o(kk) = Thermo_ou(kk)%Ptot
	Temp_o(kk) = Thermo_ou(kk)%Temp
	dpde_o(kk) = dP_de_v(kk)
	!dtde_o(kk)=1./de_dT_v(kk)
	cs2_o(kk)  =-dP_dv_s(kk)/(Thermo_ou(kk)%rho*Thermo_ou(kk)%rho)

	!   Qq verif et ajout de grandeurs 
	if (Temp_o(kk).lt.NPAR_TEMP_MIN)  conv_o(kk)=-40
	if (dpde_o(kk).lt.NPAR_DPDE_MIN) conv_o(kk)=-41
	if (cs2_o(kk).lt.NPAR_CSO2_MIN) conv_o(kk)=-42

    if (conv_o(kk).ne.0) then
        Temp_o(kk) =NPAR_TEMP_MIN ! grandeur faible et normalement non physique pouvant servir de flag pour la convergence 
        dpde_o(kk)=NPAR_DPDE_MIN
        cs2_o(kk)=NPAR_CSO2_MIN ! grandeur faible et normalement non physique pouvant servir de flag pour la convergence  
    endif

enddo

frac_1(1:nb) = Frac_ou(1:nb,1)
frac_2(1:nb) = Frac_ou(1:nb,2)
frac_3(1:nb) = Frac_ou(1:nb,3)
if (nip_ee.gt.3) then
frac_4(1:nb) = Frac_ou(1:nb,4)
frac_5(1:nb) = Frac_ou(1:nb,5)
frac_6(1:nb) = Frac_ou(1:nb,6)
endif
Pres(1:nb)=Pres_o(1:nb)
Temp(1:nb)=Temp_o(1:nb)

dpde(1:nb)=dpde_o(1:nb)
!dtde(1:nb)=dtde_o(1:nb)
cs2(1:nb) =cs2_o(1:nb)
conv(1:nb)=conv_o(1:nb)


END SUBROUTINE S_CALC_CINE_VE

END MODULE SCALC
