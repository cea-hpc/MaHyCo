MODULE M_COUPLAGE_TOOLS
CONTAINS

SUBROUTINE S_COPIE_THERMO_IN(nb,nip,Thermo_in,Frac_in,Thermo_in_sav,Frac_in_sav)
USE M_CST_TABLE_EE,       ONLY : TAB_APPEL_EE  ! type pour les propriétés thermo d'entree

IMPLICIT NONE 
INTEGER                          :: nb,nip
TYPE(TAB_APPEL_EE)               :: Thermo_in(nb)
TYPE(TAB_APPEL_EE)               :: Thermo_in_sav(nb)
REAL*8                            :: Frac_in(nb,nip)
REAL*8                             :: Frac_in_sav(nb,nip)
INTEGER                          :: kk

do kk=1,nb
  Thermo_in_sav(kk)%rho       =Thermo_in(kk)%rho
  Thermo_in_sav(kk)%ene       =Thermo_in(kk)%ene
  Thermo_in_sav(kk)%Ptot      =Thermo_in(kk)%Ptot
  Thermo_in_sav(kk)%Temp      =Thermo_in(kk)%Temp
  Thermo_in_sav(kk)%v_j(1:nip)=Thermo_in(kk)%v_j(1:nip)
  Thermo_in_sav(kk)%e_j(1:nip)=Thermo_in(kk)%e_j(1:nip)  
  Frac_in_sav(kk,1:nip)       =Frac_in(kk,1:nip)
enddo
END SUBROUTINE S_COPIE_THERMO_IN
SUBROUTINE S_COPIE_THERMO_OUT(nb,nip,rcopie,Thermo_ou_loc,Frac_ou_loc,Thermo_ou,Frac_ou)
USE M_CST_TABLE_EE,       ONLY : Thermo_Newton  ! type pour les propriétés thermo d'entree

IMPLICIT NONE 
INTEGER                          :: nb,nip
REAL*8                             :: rcopie(nb)
TYPE(Thermo_Newton)              :: Thermo_ou(nb)
TYPE(Thermo_Newton)              :: Thermo_ou_loc(nb)
REAL*8                    :: Frac_ou(nb,nip)
REAL*8                     :: Frac_ou_loc(nb,nip)
INTEGER                          :: kk
do kk=1,nb
  Thermo_ou(kk)%iconv         =int(rcopie(kk))
  
  if (int(rcopie(kk)).eq.0) then
  Thermo_ou(kk)%Ptot           =Thermo_ou_loc(kk)%Ptot
  Thermo_ou(kk)%Temp           =Thermo_ou_loc(kk)%Temp
  Thermo_ou(kk)%rho            =Thermo_ou_loc(kk)%rho
  Thermo_ou(kk)%ene            =Thermo_ou_loc(kk)%ene
  Thermo_ou(kk)%v_j(1:nip)     =Thermo_ou_loc(kk)%v_j(1:nip)
  Thermo_ou(kk)%e_j(1:nip)     =Thermo_ou_loc(kk)%e_j(1:nip)  
  Thermo_ou(kk)%dedP_T_j(1:nip)=Thermo_ou_loc(kk)%dedP_T_j(1:nip)
  Thermo_ou(kk)%dvdP_T_j(1:nip)=Thermo_ou_loc(kk)%dvdP_T_j(1:nip)  
  Thermo_ou(kk)%dedT_P_j(1:nip)=Thermo_ou_loc(kk)%dedT_P_j(1:nip)
  Thermo_ou(kk)%dvdT_P_j(1:nip)=Thermo_ou_loc(kk)%dvdT_P_j(1:nip)  
  Thermo_ou(kk)%Gibbs(1:nip)   =Thermo_ou_loc(kk)%Gibbs(1:nip)
  Frac_ou(kk,1:nip)            =Frac_ou_loc(kk,1:nip)
  endif
enddo
END SUBROUTINE S_COPIE_THERMO_OUT

SUBROUTINE S_CALC_CINE_VE_ONE(nb,nip,dtime,trans_max,&
			  Ther_in,&
              frac_in,&
              conv_in,&
			  Ther_ou,&
			  frac_ou,&
			  conv_ou)

USE M_CST_TABLE_EE,       ONLY : TAB_APPEL_EE  ! type pour les propriétés thermo d'entree
USE M_CST_TABLE_EE,       ONLY : Thermo_Newton ! type pour les propriétés thermo de sortie 
USE M_CST_TABLE_EE,       ONLY : nip_max_ee
USE M_CST_TABLE_EE,       ONLY : Frac_init
USE M_CST_TABLE_EE,       ONLY : my_HUGE
USE M_CST_TABLE_EE,       ONLY : iFlag_conv
USE M_CST_TABLE_CINE,     ONLY : TABLE_CINE_LOC    ! type pour les propriétés cine d'appel 
USE M_CST_TABLE_CINE,     ONLY : NPAR_Cine

USE M_CST_TABLE_EE,       ONLY : NPAR_TEMP_MIN
USE M_CST_TABLE_EE,       ONLY : NPAR_DPDE_MIN
USE M_CST_TABLE_EE,       ONLY : NPAR_CSO2_MIN

USE M_calc_equi_PT,       ONLY : S_CALC_CINE_PT

USE M_polygone,      ONLY : F_in_polygone
USE M_CST_TABLE_EE,  ONLY : POLYGONE

USE M_calc_equi_PT,  ONLY : S_CHECK_FRACTION



IMPLICIT NONE 
!.... INPUT:
INTEGER   :: nb
INTEGER   :: nip
REAL*8      :: dtime(nb)
REAL*8      :: trans_max
REAL*8      :: frac_in(nb,nip)
REAL*8      :: conv_in(nb)


!.... OUTPUT:
TYPE(TAB_APPEL_EE)  :: Ther_in(nb)
TYPE(Thermo_Newton) :: Ther_ou(nb)
REAL*8    			 :: frac_ou(nb,nip)
REAL*8 	 			 :: conv_ou(nb)

!.... LOCAL:
INTEGER                          :: ii,jj,kk
INTEGER                          :: niter
REAL*8                             ::fracdfrac

!.... PARTIE CINETIQUE 
REAL*8                            :: dtime_som(nb)
REAL*8                            :: dtime_loc(nb)
REAL*8                             :: R_ij_dt
INTEGER                          :: iconvergence

REAL*8                             :: DG_ij
REAL*8                             :: Fac
REAL*8, PARAMETER                  :: Fac_max=10.
TYPE(TABLE_CINE_LOC)             :: Cine_Loc

INTEGER                          :: icalc(nb)
REAL*8                          :: Frac_i(nb,nip)
REAL*8                             :: ddtime

INTEGER, PARAMETER               :: niter_max=5
INTEGER                          :: niter_max_eff ! nb iteration max = niter_max/trans_max

! INITIALISATION DU CALCUL

! Affectation des fraction des mailles
Frac_i(1:nb,1:nip)=frac_in(1:nb,1:nip)  


dtime_som(1:nb)=0
dtime_loc(1:nb)=dtime(1:nb)


do kk=1,nb
  if (conv_in(kk).eq.0) then ! convergence avant 
    dtime_som(kk)=2.*dtime(kk)   ! on a convergé 
  endif
enddo

niter_max_eff=niter_max*int(1./trans_max)

niter=0
iconvergence=-1
do while (iconvergence.ne.0) 
    ! on verifie que l'on a convergé sur le pas de temps 
    iconvergence=0
    do kk=1,nb
      icalc(kk)=1
      !verifie si le pas de temps > dtime du code hydro
      if (dtime_som(kk).ge.dtime(kk)) then
        icalc(kk)=0
      endif   
      if (icalc(kk).ne.0) then
        iconvergence=-1
      endif
      conv_ou(kk)=icalc(kk)
    enddo
    ! tout est ok, on sort c'est la sortie normale    
    if (iconvergence.eq.0) exit 
    ! Appel à l'équilibre PT pour les Frac_ni
    ! Seuls les points ical!=0 seront calculées 
    ! La verification des fractions entre 0 et 1 
    ! se fait dans S_CHECK_FRACTION inclus dans S_CALC_CINE_PT
    
    !write(*,*) " on repart avec pression", Ther_in(1)%Ptot, " Temp ", Ther_in(1)%Temp
    !write(*,*) " on repart avec energie", Ther_in(1)%ene, " densité ", Ther_in(1)%rho
    CALL S_CALC_CINE_PT(nbmail=nb,nip=nip,&
			Prop_pas_ni=Ther_in(1:nb),&
			Prop_pas_nf=Ther_ou(1:nb),&
			Frac_ni=Frac_i(1:nb,1:nip),&
			icalc=icalc(1:nb)&
			)
			
    !write(*,*) " on resort avec pression", Ther_ou(1)%Ptot, " Temp ", Ther_ou(1)%Temp
    !write(*,*) " on resort avec energie", Ther_ou(1)%ene, " densité ", Ther_ou(1)%rho
    ! Affectation de l'énergie libre 
    do kk=1,nb
        !if (kk.eq.1) write(*,*) 'iteration cv', niter, ' Fraction à l entrée pour la maille ', kk, " =" , frac_in(kk,1:nip)
        !if (kk.eq.1) write(*,*) 'avec iconv (in )', Ther_in(kk)
        !if (kk.eq.1) write(*,*) 'avec iconv (out)', Ther_ou(kk)%iconv
        !if (kk.eq.1) write(*,*) 'cinetique (0 hayes, 1 Greeff) :', NPAR_cine%type_cinetique
 
      if (Ther_ou(kk)%iconv.lt.0) then
        ! on ne recalcul plus ca ne sert à rien
        icalc(kk)=0
        dtime_som(kk)=2*dtime(kk)
      else
        ! n'appartient pas à la zone d''existence du polygone 
        do jj=1,nip
          if ((F_in_polygone(POLYGONE(jj)%nb_pts,POLYGONE(jj)%PT(1:POLYGONE(jj)%nb_pts,1:2),&
           Ther_ou(kk)%Ptot,Ther_ou(kk)%Temp))) then 
          else  
            Ther_ou(kk)%Gibbs(jj)=my_HUGE    
          endif
        enddo
        
        do ii=1,nip
          do jj=1,nip           
            Cine_Loc%R_ij(ii,jj)=0.
            DG_ij=Ther_ou(kk)%Gibbs(ii)-Ther_ou(kk)%Gibbs(jj) + NPAR_cine%shift_Gii(ii)
            if (NPAR_cine%type_cinetique.eq.1) then 
                ! Formulation de Greeff Rij=nu_ij*H(DG)*(DG/B_ij)*exp((DG/B_ij)²)
                if (DG_ij.gt.0.) then 
                    Fac=MAX(MIN(DG_ij/NPAR_Cine%B_ij(ii,jj),Fac_Max),-Fac_Max)
                    Cine_Loc%R_ij(ii,jj)=NPAR_Cine%nu_ij(ii,jj)*Fac*exp(Fac*Fac)
                    !if (kk.eq.1) write(*,*) 'Greeff Fac de phase ', ii, ' et ', jj, '=', Fac, ' et DG_ij ', DG_ij   
                    !if (kk.eq.1) write(*,*) ' Cine_Loc%R_ij(ii,jj) = ' , Cine_Loc%R_ij(ii,jj)
                    !if (kk.eq.1) write(*,*) 'Ther_ou(kk)%Gibbs(ii)', Ther_ou(kk)%Gibbs(ii)
                    !if (kk.eq.1) write(*,*) 'Ther_ou(kk)%Gibbs(jj)', Ther_ou(kk)%Gibbs(jj) 
                    !if (kk.eq.1) write(*,*) 'NPAR_cine%shift_Gii(ii)', NPAR_cine%shift_Gii(ii)
                endif
            else 
                ! Formulation de Hayes Rij=nu_ij*H(DG)*(DG/B_ij)
                if (DG_ij.gt.0.) then 
                    Fac= DG_ij/NPAR_Cine%B_ij(ii,jj)
                    Cine_Loc%R_ij(ii,jj)=NPAR_Cine%nu_ij(ii,jj)*Fac
                    !if (kk.eq.1) write(*,*) 'Hayes Fac de phase ', ii, ' et ', jj, '=', Fac, ' et DG_ij ', DG_ij   
                    !if (kk.eq.1) write(*,*) ' Cine_Loc%R_ij(ii,jj) = ' , Cine_Loc%R_ij(ii,jj)
                endif
            endif   
          enddo
        enddo
        ! calcul des dfrac sans le temps
        do ii=1,nip
        Cine_Loc%dFrac_i(ii)=0.
          do jj=1,nip
            Cine_Loc%dFrac_i(ii)=Cine_Loc%dFrac_i(ii)+&
                        (Frac_i(kk,jj)*Cine_Loc%R_ij(jj,ii)-&
                        Frac_i(kk,ii)*Cine_Loc%R_ij(ii,jj)&
                        )
          enddo 
        enddo
        
        ! calcul du temps local 
        dtime_loc(kk)=dtime(kk)-dtime_som(kk)
        
        !if (kk.eq.1) write(*,*) ' AV dtime_loc ' , dtime_loc(kk), ' ddtime ' , ddtime, ' et trans_max' , trans_max
        do ii=1,nip
            !if (kk.eq.1) write(*,*) ' cineloc*dtloc = ' , abs(Cine_Loc%dFrac_i(ii))*dtime_loc(kk)
            !if (kk.eq.1) write(*,*) ' avec cine ' , Cine_Loc%dFrac_i(ii)
          if (abs(Cine_Loc%dFrac_i(ii))*dtime_loc(kk).gt.trans_max) then
            !ddtime=max(trans_max/abs(Cine_Loc%dFrac_i(ii)),1.e-24) ! borne indispensable à paramétrer
            ddtime=trans_max/abs(Cine_Loc%dFrac_i(ii))
            dtime_loc(kk)=min(ddtime,dtime_loc(kk))
          endif
        enddo
        !if (kk.eq.1) write(*,*) ' AP dtime_loc ' , dtime_loc(kk), ' ddtime ' , ddtime
        do ii=1,nip
          Frac_i(kk,ii)=Frac_i(kk,ii)+dtime_loc(kk)*Cine_Loc%dFrac_i(ii)
        enddo
        dtime_som(kk)=dtime_som(kk)+dtime_loc(kk)
        !if (kk.eq.1) write(*,*) ' dtime_som ' , dtime_som (kk), ' dtime', dtime(kk) , ' dtime_loc ' , dtime_loc(kk)
        icalc(kk)=Ther_ou(kk)%iconv
        Ther_in(kk)%Ptot=Ther_ou(kk)%Ptot
        Ther_in(kk)%Temp=Ther_ou(kk)%Temp
        Ther_in(kk)%v_j(1:nip)=Ther_ou(kk)%v_j(1:nip)
        Ther_in(kk)%e_j(1:nip)=Ther_ou(kk)%e_j(1:nip)
        !if (kk.eq.1) then 
        !   write(*,'(a, 6e22.14)') 'DFraction pour la maille 1 ', Cine_Loc%dFrac_i
        !   write(*,'(a, 6e22.14)') 'FRac+DFraction ', Frac_i(kk,1)+dtime_loc(kk)*Cine_Loc%dFrac_i(1)
        !   write(*,*) ' et dtloc ', dtime_loc(kk) , ' et conv ' , Ther_ou(kk)%iconv, ' et nip', nip
        !endif    
      endif
    enddo
niter=niter+1
if (niter.gt.niter_max_eff) exit 
enddo

conv_ou(1:nb)=Ther_ou(1:nb)%iconv
do kk=1,nb
  if (conv_in(kk).eq.0) then ! convergence avant 
    conv_ou(kk)=iFlag_conv         ! on avait convergé, on ne change rien 
  endif
enddo


do kk=1,nb
  if (conv_ou(kk).eq.0) then 
    ! renormalisation pour la sortie 
    CALL S_CHECK_FRACTION(nip,Frac_i(kk,1:nip))
    ! Affectation des fraction des mailles
    frac_ou(kk,1:nip)=Frac_i(kk,1:nip)
    ! conv_ou Ther_in et Ther_ou est déja renseigné en cours de calcul si Ther_ou(1:nb)%iconv =0
    !if (kk.eq.1) write(*,'(a, i4, a, 6e22.14)') 'iteration cv', niter , ' Fraction pour la maille 1 ', frac_ou(1,1:nip)
  endif
enddo

!write(*,*) 'fin de S_CALC_CINE_VE_ONE' , ' iteration ', niter

END SUBROUTINE S_CALC_CINE_VE_ONE


END MODULE M_COUPLAGE_TOOLS
