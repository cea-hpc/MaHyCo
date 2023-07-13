MODULE M_ee_Sinit
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Routine d'initialisation appelée par le fichier .st
! ..... Module à supprimer et NCONST à mettre à 1 sans passage si modèle
! ..... composite
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE CineTest22_ee_Sinit(fcoeff,rho_std,ene_std,Nconst)
  USE M_CST
  USE M_CST_TABLE_CINE, ONLY : nip_cine

  USE M_CST_TABLE_EE,  ONLY : nip_max_ee,TB_rT,TB_re,TB_PT,LIMIT_TABLE
  USE M_CST_TABLE_EE,  ONLY : my_HUGE
  USE M_CST_TABLE_EE,  ONLY : POLYGONE
  USE M_CST_TABLE_EE,  ONLY : nip_ee,Frac_init
  USE M_CST_TABLE_EE,  ONLY : eps_init
  
  USE M_CST_TABLE_EE,  ONLY : nb_restart,restart_PT
  
  USE M_init_tools,    ONLY : F_R_PT
  USE M_init_tools,    ONLY : S_DEF_GIBBS
  USE M_init_tools,    ONLY : S_DEF_CBF
  USE M_init_tools,    ONLY : S_DEFTAB_RE
  USE M_init_tools,    ONLY : S_DEFTAB_PT
  
  USE M_calc_table_rt, ONLY : S_calc_thermo_one_rt
  USE M_calc_table_re, ONLY : S_calc_thermo_one_re
  USE M_calc_table_pt, ONLY : S_calc_thermo_one_PT

  USE M_polygone,      ONLY : F_in_polygone

  USE M_calc_equi_PT,  ONLY : S_CHECK_FRACTION

  IMPLICIT NONE


! ..... Buffers : reel
  REAL*8        ::  fcoeff(*)
  INTEGER     ::  Nconst

! ..... Compteurs 
  INTEGER     :: ii,jj,kk,ll
  INTEGER     :: djj
! ..... Masse molaire local --> copie dans TB_rt(ii)%Std%M_mol une fois le tableau crée 
  REAL*8  :: M_mol

  LOGICAL :: Linit_frac
  REAL*8    :: norm

  INTEGER,PARAMETER         :: i_rho=3 ! emplacement du Gibbs dans val_out
  INTEGER,PARAMETER         :: i_ene=8 ! emplacement du Gibbs dans val_out
  INTEGER,PARAMETER         :: i_Gibbs=14 ! emplacement du Gibbs dans val_out
  REAL*8  		    :: valeurs(CST_Nval_max)
 
  REAL*8                      :: Pstd(1)
  REAL*8                      :: Tstd(1)
  REAL*8                      :: Gibbs(nip_max_ee)
  REAL*8                      :: vol_j(nip_max_ee)
  REAL*8                      :: ene_j(nip_max_ee)

  REAL*8                      :: rho_std,ene_std

  write(*,*) 'Rentré dans CineTest22_ee_Sinit'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Initialisation par les jeux de coefficients inclus dans le buffer
! ..... Une adapation est à prévoir entre les indices pour fcoeff   
! ..... Nb de phase 
  nip_ee=int(fcoeff(Nconst));Nconst=Nconst+1
  write(*,*) 'nip_ee', nip_ee



! ..... Verif nb phases max 
  IF (nip_ee.gt.nip_max_ee) THEN
  write(*,'(a)') 'Erreur nombre de phases stop'
  stop
  ENDIF

! ..... Verif coherence nb de phases ee-cine
  if ((nip_ee.gt.0).and.(nip_cine.gt.0)) then
    if (nip_ee.ne.nip_cine) then
      write(*,'(a)') 'Erreur definition non coherente du nb de phases STOP'
      STOP
    endif
  endif
  
 ! ..... Masse molaire 
  M_mol=fcoeff(Nconst);Nconst=Nconst+1
  write(*,*) 'M_mol', M_mol

! ..... Allocation mémoire  
! ..... Enregistrement des fractions initiales 
  IF (ALLOCATED(Frac_init))      DEALLOCATE(Frac_init);
  ALLOCATE(Frac_init(nip_ee))
  !Frac_init(:)=0.
  norm=0.
  DO ii=1,nip_ee
    Frac_init(ii)=fcoeff(Nconst)
    norm=norm+Frac_init(ii)
    Nconst=Nconst+1
  ENDDO
  if (MINVAL(Frac_init).lt.0.) then
  ! il y a une valeur < 0.
    Linit_frac=.TRUE.
    Frac_init(:)=0.
  elseif   (SUM(Frac_init(1:nip_ee)).lt.eps_init) then 
    Linit_frac=.TRUE.
    Frac_init(:)=0.
  else
  ! il y a aucune valeur < 0.
    Frac_init(:)=Frac_init(:)/norm
    Linit_frac=.FALSE.
  endif

! ..... Allocation mémoire  
! ..... Enregistrement des polygones 
  IF (ALLOCATED(POLYGONE))      DEALLOCATE(POLYGONE);
  ALLOCATE(POLYGONE(nip_ee))
  DO ii=1,nip_ee
    POLYGONE(ii)%nb_pts=int(fcoeff(Nconst))
    Nconst=Nconst+1
  ENDDO
  
  DO kk=1,2
    DO ii=1,nip_ee
      DO jj=1,POLYGONE(ii)%nb_pts
        POLYGONE(ii)%PT(jj,kk)=fcoeff(Nconst)
	Nconst=Nconst+1
      ENDDO
    ENDDO
  ENDDO

! ..... Point pour repartir si prb convergence 
  nb_restart=int(fcoeff(Nconst))
  Nconst=Nconst+1
  IF (ALLOCATED(restart_PT))      DEALLOCATE(restart_PT);
  ALLOCATE(restart_PT(nb_restart,2))
  
  DO kk=1,2
      DO jj=1,nb_restart
        restart_PT(jj,kk)=fcoeff(Nconst)
	Nconst=Nconst+1
    ENDDO
  ENDDO
 
 ! ..... Allocation mémoire  
  IF (ALLOCATED(TB_rt))      DEALLOCATE(TB_rt);
  ALLOCATE(TB_rt(nip_ee))
! ..... Taille tableau de chaque phase 
  DO ii=1,nip_ee
    TB_rt(ii)%nr=int(fcoeff(Nconst))
    Nconst=Nconst+1
  ENDDO
  DO ii=1,nip_ee
    TB_rt(ii)%nT=int(fcoeff(Nconst))
    Nconst=Nconst+1
  ENDDO
  

! ..... Allocation mémoire tableaux 
  DO ii=1,nip_ee
    IF (ALLOCATED(TB_rt(ii)%rho))   DEALLOCATE(TB_rt(ii)%rho);
    ALLOCATE(TB_rt(ii)%rho(TB_rt(ii)%nr))
    TB_rt(ii)%rho(:)=0.
    IF (ALLOCATED(TB_rt(ii)%temp))   DEALLOCATE(TB_rt(ii)%temp);
    ALLOCATE(TB_rt(ii)%temp(TB_rt(ii)%nT))
    TB_rt(ii)%temp(:)=0.
    IF (ALLOCATED(TB_rt(ii)%Ptot))   DEALLOCATE(TB_rt(ii)%Ptot);
    ALLOCATE(TB_rt(ii)%Ptot(TB_rt(ii)%nr*TB_rt(ii)%nT))
    TB_rt(ii)%Ptot(:)=0.
    IF (ALLOCATED(TB_rt(ii)%etot))   DEALLOCATE(TB_rt(ii)%etot);
    ALLOCATE(TB_rt(ii)%etot(TB_rt(ii)%nr*TB_rt(ii)%nT))
    TB_rt(ii)%etot(:)=0.
    IF (ALLOCATED(TB_rt(ii)%stot))   DEALLOCATE(TB_rt(ii)%stot);
    ALLOCATE(TB_rt(ii)%stot(TB_rt(ii)%nr*TB_rt(ii)%nT))
    TB_rt(ii)%stot(:)=0.
  ENDDO

! ..... Allocations tableaux en densite 
  DO ii=1,nip_ee
    djj=TB_rt(ii)%nr
    TB_rt(ii)%rho(1:djj)=fcoeff(Nconst:Nconst+djj-1)
    Nconst=Nconst+djj
  ENDDO
! ..... Allocations tableaux en temperature
  DO ii=1,nip_ee
    djj=TB_rt(ii)%nT
    TB_rt(ii)%temp(1:djj)=fcoeff(Nconst:Nconst+djj-1)
    Nconst=Nconst+djj
  ENDDO
! ..... Allocations tableaux en Pression
  DO ii=1,nip_ee
    djj=TB_rt(ii)%nr*TB_rt(ii)%nT
    TB_rt(ii)%Ptot(1:djj)=fcoeff(Nconst:Nconst+djj-1)
    Nconst=Nconst+djj
  ENDDO
! ..... Allocations tableaux en energie interne 
  DO ii=1,nip_ee
    djj=TB_rt(ii)%nr*TB_rt(ii)%nT
    TB_rt(ii)%etot(1:djj)=fcoeff(Nconst:Nconst+djj-1)
    Nconst=Nconst+djj
  ENDDO
! ..... Allocations tableaux en entropie
  DO ii=1,nip_ee
    djj=TB_rt(ii)%nr*TB_rt(ii)%nT
    TB_rt(ii)%stot(1:djj)=fcoeff(Nconst:Nconst+djj-1)
    Nconst=Nconst+djj
  ENDDO
! ..... Fin de la lecture des coefficients 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  write(*,*) 'ici c est bon 5' 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Allocations tableaux et calcul du Gibbs
  DO ii=1,nip_ee
    CALL S_DEF_GIBBS(Tab_rT=TB_rt(ii))
  ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Allocations de la courbe de reference 
  DO ii=1,nip_ee
    CALL S_DEF_CBF(Tab_rT=TB_rt(ii))
  ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Initialisation du P et T état standard
  DO ii=1,nip_ee
    TB_rt(ii)%Std%Ptot =CST_Pstd
    TB_rt(ii)%Std%Temp =CST_Tstd
    TB_rt(ii)%Std%M_mol=M_mol
  ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Recherche du r0 de l'état standard P0 et T0
  DO ii=1,nip_ee
    TB_rt(ii)%Std%rho=F_R_PT(P_ref=TB_rt(ii)%Std%Ptot,T_ref=TB_rt(ii)%Std%Temp,Tab_rT=TB_rt(ii))
  ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Construction de la table de transcription T(ve)
! ..... Allocation mémoire  
  IF (ALLOCATED(TB_re))      DEALLOCATE(TB_re);
  ALLOCATE(TB_re(nip_ee))
! ..... Allocation mémoire et conversion en T(ve)
  DO ii=1,nip_ee
    TB_re(ii)%Cv_3kb=CST_3*CST_kb/TB_rt(ii)%Std%M_mol
    CALL S_DEFTAB_RE(TB_rt(ii),TB_re(ii))
  ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Construction de la table de transcription P(vT)
! ..... Allocation mémoire  
  IF (ALLOCATED(TB_PT))      DEALLOCATE(TB_PT);
  ALLOCATE(TB_PT(nip_ee))
! ..... Allcoation mémoire et conversion en P(vT)
  DO ii=1,nip_ee
    CALL S_DEFTAB_PT(TB_rt(ii),TB_PT(ii))
  ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Definition des bornes v et T 
! ..... Allocation mémoire  
  IF (ALLOCATED(LIMIT_TABLE))      DEALLOCATE(LIMIT_TABLE);
  ALLOCATE(LIMIT_TABLE(nip_ee))
! ..... Allcoation mémoire et conversion en P(vT)
  DO ii=1,nip_ee
  
    LIMIT_TABLE(ii)%TMIN=TB_rt(ii)%Temp(1)
    LIMIT_TABLE(ii)%TMAX=TB_rt(ii)%Temp(TB_rt(ii)%nT)
    LIMIT_TABLE(ii)%VMIN=1./TB_rt(ii)%rho(TB_rt(ii)%nr)
    LIMIT_TABLE(ii)%VMAX=1./TB_rt(ii)%rho(1)
  ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pstd(1)=CST_PStd
    Tstd(1)=CST_TStd
    do jj=1,nip_ee
	CALL S_calc_thermo_one_PT(Pre=Pstd(1),Temp=Tstd(1),Tab_PT=TB_PT(jj),Tab_rT=TB_rT(jj),prop='all',val_out=valeurs)
    enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Calcul de l'état Std
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Pstd(1)=CST_PStd
Tstd(1)=CST_TStd

do jj=1,nip_ee
  Gibbs(jj)=my_HUGE
  if (F_in_polygone(POLYGONE(jj)%nb_pts,POLYGONE(jj)%PT(1:POLYGONE(jj)%nb_pts,1:2),Pstd(1),Tstd(1))) then 
    CALL S_calc_thermo_one_PT(Pre=Pstd(1),Temp=Tstd(1),Tab_PT=TB_PT(jj),Tab_rT=TB_rT(jj),prop='all',val_out=valeurs)
    Gibbs(jj)=valeurs(i_Gibbs)
    ene_j(jj)=valeurs(i_ene)
    vol_j(jj)=1./valeurs(i_rho)
  endif
enddo


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Calcul des fractions de l'état Std au besoin 
if (Linit_frac) then 
    Frac_init(MINLOC(Gibbs(1:nip_ee),DIM=1))=1.
endif
CALL S_CHECK_FRACTION(nip_ee,dble(Frac_init(1:nip_ee)))

  write(*,*) 'ici c est bon 6', Frac_init(1)
  
ene_std=0.
rho_std=0.
do jj=1,nip_ee
ene_std=ene_std+Frac_init(jj)*ene_j(jj)
rho_std=rho_std+Frac_init(jj)*vol_j(jj)
enddo
rho_std=1./rho_std
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... La table est connue dans les formulations en ve,vT et PT
! ..... Normalement tout est prêt pour la cinétique 
! ..... Les phases stables pour l'état standard sont dans Frac_init 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write(*,*) 'Valeur standard', ene_std, rho_std

END SUBROUTINE CineTest22_ee_Sinit
 
END MODULE M_ee_Sinit
