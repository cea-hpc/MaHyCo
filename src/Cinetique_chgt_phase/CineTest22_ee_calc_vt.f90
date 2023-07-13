MODULE M_calc_table_rt
CONTAINS 

INTEGER FUNCTION F_loc_X(irho,XX,nbX,XT)
! ..... Fonction qui renvoit la position dans le tableau 
USE M_interpol_EOS, ONLY : locate_eos
IMPLICIT NONE 
INTEGER irho,nbX
REAL*8    XX,XT(:)
! ..... 
F_loc_X=irho
if (F_loc_X.le.0) then 
    call locate_eos(XT(1:nbX),nbX,XX,F_loc_X)
    if (F_loc_X.lt.2)	          F_loc_X=2
    if (F_loc_X.gt.nbX-2)         F_loc_X=nbX-2
endif
return 

END FUNCTION F_loc_X


SUBROUTINE S_calc_thermo_one_rt(rho,Temp,Tab_rT,prop,irho_old,jtemp_old,val_out)
! ..... SUBROUTINE de calcul pour un point  
  USE M_interpol_EOS, ONLY : locate_eos
  USE M_interpol_EOS, ONLY : S_INTERPOLV2
  USE M_CST_TABLE_EE
  IMPLICIT NONE 
  INTEGER, PARAMETER:: nbmail=1               ! une seule maille à la fois 
  INTEGER, PARAMETER:: imail=1               ! une seule maille à la fois 

  TYPE(TABLE_RT) :: Tab_rT

  REAL*8          :: rho(nbmail),Temp(nbmail)   ! densite,temperature
  REAL*8          :: val_out(:)                 ! valeur de sortie 
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

  CHARACTER(LEN=*),  OPTIONAL :: prop 
! ..... si présent et positif, pas de locate 
  INTEGER,           OPTIONAL :: irho_old(nbmail),jtemp_old(nbmail)

  INTEGER           :: irhoi(nbmail)
  INTEGER           :: jtempi(nbmail)
  REAL*8              :: X_m(nbmail), dXdv_m(nbmail), dXdt_m(nbmail)
  INTEGER           :: ij_prop

  INTEGER           :: i_dPdv_T,i_dPdT_r,i_dedv_T,i_dedT_r,i_dsdv_T,i_dsdT_r,i_dGdv_T,i_dGdT_r,i_dvdT_P


  
  if (prop.eq.'irho') then
    !localisation table r
    ij_prop=1
    irhoi(imail)=-1
    if (present(irho_old)) irhoi(imail)=irho_old(imail)
    irhoi(imail)=F_loc_X(irhoi(imail),rho(imail),Tab_rT%nr,Tab_rT%rho(1:Tab_rT%nr))
    val_out(ij_prop)=1.*irhoi(imail)
    return
  elseif (prop.eq.'jtemp') then
    !localisation table T
    ij_prop=2
    jtempi(imail)=-1
    if (present(jtemp_old)) jtempi(imail)=jtemp_old(imail)
    jtempi(imail)=F_loc_X(jtempi(imail),Temp(imail),Tab_rT%nT,Tab_rT%temp(1:Tab_rT%nT))
    val_out(ij_prop)=1.*jtempi(imail)
    return
  else 
    !localisation table r
    ij_prop=1
    irhoi(imail)=-1
    if (present(irho_old)) irhoi(imail)=irho_old(imail)
    irhoi(imail)=F_loc_X(irhoi(imail),rho(imail),Tab_rT%nr,Tab_rT%rho(1:Tab_rT%nr))
    val_out(ij_prop)=1.*irhoi(imail)
    !localisation table T
    ij_prop=ij_prop+1
    jtempi(imail)=-1
    if (present(jtemp_old)) jtempi(imail)=jtemp_old(imail)
    jtempi(imail)=F_loc_X(jtempi(imail),Temp(imail),Tab_rT%nT,Tab_rT%temp(1:Tab_rT%nT))
    val_out(ij_prop)=1.*jtempi(imail)
  endif

! ..... On renvoit toujours rho et T
    ij_prop=3
    val_out(ij_prop)=rho(imail) 
    ij_prop=4
    val_out(ij_prop)=Temp(imail)

  if (prop.eq.'Ptot') then
! ..... Ptot
    ij_prop=5
    CALL S_INTERPOLV2(Tab_rT%nr,Tab_rT%nT,Tab_rT%rho(1:Tab_rT%nr),Tab_rT%temp(1:Tab_rT%nT),&
		      Tab_rT%Ptot(1:Tab_rT%nr*Tab_rT%nT), &
		      nbmail,rho(imail),Temp(imail),&
		      X_m(imail), dXdv_m(imail), dXdT_m(imail),&
		      irhoi(imail), jtempi(imail));
    dXdv_m(imail)=-dXdv_m(imail)*rho(imail)*rho(imail)
    val_out(ij_prop)=X_m(imail)    ; ij_prop=ij_prop+1 ! 5 Ptot
    val_out(ij_prop)=dXdv_m(imail) ; ij_prop=ij_prop+1 ! 6 dP_dv|T
    val_out(ij_prop)=dXdT_m(imail) ; ij_prop=ij_prop+1 ! 7 dP_dT|v
  elseif (prop.eq.'etot') then
! ..... etot
    ij_prop=8
    CALL S_INTERPOLV2(Tab_rT%nr,Tab_rT%nT,Tab_rT%rho(1:Tab_rT%nr),Tab_rT%temp(1:Tab_rT%nT),&
		      Tab_rT%etot(1:Tab_rT%nr*Tab_rT%nT), &
		      nbmail,rho(imail),Temp(imail),&
		      X_m(imail), dXdv_m(imail), dXdT_m(imail),&
		      irhoi(imail), jtempi(imail));
    dXdv_m(imail)=-dXdv_m(imail)*rho(imail)*rho(imail)
    val_out(ij_prop)=X_m(imail)    ; ij_prop=ij_prop+1 ! 8 etot
    val_out(ij_prop)=dXdv_m(imail) ; ij_prop=ij_prop+1 ! 9 de_dv|T
    val_out(ij_prop)=dXdT_m(imail) ; ij_prop=ij_prop+1 !10 de_dT|v
  elseif (prop.eq.'stot') then
! ..... stot
    ij_prop=11
    CALL S_INTERPOLV2(Tab_rT%nr,Tab_rT%nT,Tab_rT%rho(1:Tab_rT%nr),Tab_rT%temp(1:Tab_rT%nT),&
		      Tab_rT%stot(1:Tab_rT%nr*Tab_rT%nT), &
		      nbmail,rho(imail),Temp(imail),&
		      X_m(imail), dXdv_m(imail), dXdT_m(imail),&
		      irhoi(imail), jtempi(imail));
    dXdv_m(imail)=-dXdv_m(imail)*rho(imail)*rho(imail)
    val_out(ij_prop)=X_m(imail)    ; ij_prop=ij_prop+1 ! 11 stot
    val_out(ij_prop)=dXdv_m(imail) ; ij_prop=ij_prop+1 ! 12 ds_dv|T
    val_out(ij_prop)=dXdT_m(imail) ; ij_prop=ij_prop+1 ! 13 ds_dT|v
  elseif (prop.eq.'Gibbs') then
! .....Gibbs
    ij_prop=14
    CALL S_INTERPOLV2(Tab_rT%nr,Tab_rT%nT,Tab_rT%rho(1:Tab_rT%nr),Tab_rT%temp(1:Tab_rT%nT),&
		      Tab_rT%Gibbs(1:Tab_rT%nr*Tab_rT%nT), &
		      nbmail,rho(imail),Temp(imail),&
		      X_m(imail), dXdv_m(imail), dXdT_m(imail),&
		      irhoi(imail), jtempi(imail));
    dXdv_m(imail)=-dXdv_m(imail)*rho(imail)*rho(imail)
    val_out(ij_prop)=X_m(imail)    ; ij_prop=ij_prop+1 ! 14 Gibbs  
    val_out(ij_prop)=dXdv_m(imail) ; ij_prop=ij_prop+1 ! 15 dG_dv|T
    val_out(ij_prop)=dXdT_m(imail) ; ij_prop=ij_prop+1 ! 16 dG_dT|v
  else
    ij_prop=5
! ..... Ptot
    CALL S_INTERPOLV2(Tab_rT%nr,Tab_rT%nT,Tab_rT%rho(1:Tab_rT%nr),Tab_rT%temp(1:Tab_rT%nT),&
		      Tab_rT%Ptot(1:Tab_rT%nr*Tab_rT%nT), &
		      nbmail,rho(imail),Temp(imail),&
		      X_m(imail), dXdv_m(imail), dXdT_m(imail),&
		      irhoi(imail), jtempi(imail));
    dXdv_m(imail)=-dXdv_m(imail)*rho(imail)*rho(imail)
    val_out(ij_prop)=X_m(imail)    ; ij_prop=ij_prop+1 ! 5 Ptot
    val_out(ij_prop)=dXdv_m(imail) ; ij_prop=ij_prop+1 ! 6 dP_dv|T
    val_out(ij_prop)=dXdT_m(imail) ; ij_prop=ij_prop+1 ! 7 dP_dT|v
! ..... etot
    CALL S_INTERPOLV2(Tab_rT%nr,Tab_rT%nT,Tab_rT%rho(1:Tab_rT%nr),Tab_rT%temp(1:Tab_rT%nT),&
		      Tab_rT%etot(1:Tab_rT%nr*Tab_rT%nT), &
		      nbmail,rho(imail),Temp(imail),&
		      X_m(imail), dXdv_m(imail), dXdT_m(imail),&
		      irhoi(imail), jtempi(imail));
    dXdv_m(imail)=-dXdv_m(imail)*rho(imail)*rho(imail)
    val_out(ij_prop)=X_m(imail)    ; ij_prop=ij_prop+1 ! 8 etot
    val_out(ij_prop)=dXdv_m(imail) ; ij_prop=ij_prop+1 ! 9 de_dv|T
    val_out(ij_prop)=dXdT_m(imail) ; ij_prop=ij_prop+1 !10 de_dT|v
! ..... stot
    CALL S_INTERPOLV2(Tab_rT%nr,Tab_rT%nT,Tab_rT%rho(1:Tab_rT%nr),Tab_rT%temp(1:Tab_rT%nT),&
		      Tab_rT%stot(1:Tab_rT%nr*Tab_rT%nT), &
		      nbmail,rho(imail),Temp(imail),&
		      X_m(imail), dXdv_m(imail), dXdT_m(imail),&
		      irhoi(imail), jtempi(imail));
    dXdv_m(imail)=-dXdv_m(imail)*rho(imail)*rho(imail)
    val_out(ij_prop)=X_m(imail)    ; ij_prop=ij_prop+1 ! 11 stot
    val_out(ij_prop)=dXdv_m(imail) ; ij_prop=ij_prop+1 ! 12 ds_dv|T
    val_out(ij_prop)=dXdT_m(imail) ; ij_prop=ij_prop+1 ! 13 ds_dT|v
    CALL S_INTERPOLV2(Tab_rT%nr,Tab_rT%nT,Tab_rT%rho(1:Tab_rT%nr),Tab_rT%temp(1:Tab_rT%nT),&
		      Tab_rT%Gibbs(1:Tab_rT%nr*Tab_rT%nT), &
		      nbmail,rho(imail),Temp(imail),&
		      X_m(imail), dXdv_m(imail), dXdT_m(imail),&
		      irhoi(imail), jtempi(imail));
    dXdv_m(imail)=-dXdv_m(imail)*rho(imail)*rho(imail)
    val_out(ij_prop)=X_m(imail)    ; ij_prop=ij_prop+1 ! 14 Gibbs  
    val_out(ij_prop)=dXdv_m(imail) ; ij_prop=ij_prop+1 ! 15 dG_dv|T
    val_out(ij_prop)=dXdT_m(imail) ; ij_prop=ij_prop+1 ! 16 dG_dT|v
    
  ! ..... prop='Ptot'    : calcul de Ptot, dP/dv|T et dP/dT|v -> val_out(5:7)
  ! ..... prop='etot'    : calcul de etot, de/dv|T et de/dT|v -> val_out(8:10)
  ! ..... prop='stot'    : calcul de stot, ds/dv|T et ds/dT|v -> val_out(11:13)
  ! ..... prop='Gibbs'   : calcul de G, dG/dv|T et dG/dT|v    -> val_out(14:16)
  ! ..... autre          : précedent (dans l'ordre) 
  ! ..... autre          : + dv/dP|T,dv/dT|P,de/dP|T,de/dT|P -> val_out(17:20)
  ! ..... autre          : + ds/dP|T,ds/dT|P,dG/dP|T,dG/dT|P -> val_out(21:24)
   
    i_dPdv_T=6
    i_dPdT_r=7
    i_dedv_T=9
    i_dedT_r=10
    i_dsdv_T=12
    i_dsdT_r=13
    i_dGdv_T=15
    i_dGdT_r=16
    i_dvdT_P=18
    
! ..... 17/ dvdP_T=1./dPdv_T   
    val_out(ij_prop)=1./val_out(i_dPdv_T)                                  ; ij_prop=ij_prop+1
! ..... 18/ dvdT_P=-dPdT_r/dPdv_T
    val_out(ij_prop)=-val_out(i_dPdT_r)/val_out(i_dPdv_T)                  ; ij_prop=ij_prop+1

! ..... 19/ dedP_T=dedv_T/dPdv_T
    val_out(ij_prop)=val_out(i_dedv_T)/val_out(i_dPdv_T)                   ; ij_prop=ij_prop+1
! ..... 20/ dedT_P=dedT_r+dvdT_P*dedv_T
    val_out(ij_prop)=val_out(i_dedT_r)+val_out(i_dvdT_P)*val_out(i_dedv_T) ; ij_prop=ij_prop+1

! ..... 21/ dsdP_T=dsdv_T/dPdv_T
    val_out(ij_prop)=val_out(i_dsdv_T)/val_out(i_dPdv_T)                    ; ij_prop=ij_prop+1
! ..... 22/ dsdT_P=dsdT_r+dvdT_P*dsdv_T
    val_out(ij_prop)=val_out(i_dsdT_r)+val_out(i_dvdT_P)*val_out(i_dsdv_T)  ; ij_prop=ij_prop+1

! ..... 23/ dGdP_T=dGdv_T/dPdv_T
    val_out(ij_prop)=val_out(i_dGdv_T)/val_out(i_dPdv_T)                    ; ij_prop=ij_prop+1
! ..... 24/ dGdT_P=dGdT_r+dvdT_P*dGdv_T
    val_out(ij_prop)=val_out(i_dGdT_r)+val_out(i_dvdT_P)*val_out(i_dGdv_T)  ; ij_prop=ij_prop+1


  endif
    

END SUBROUTINE S_calc_thermo_one_rt



END MODULE M_calc_table_rt



