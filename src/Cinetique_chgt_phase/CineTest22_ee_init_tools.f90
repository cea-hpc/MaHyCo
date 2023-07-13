MODULE M_init_tools
CONTAINS 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Fonction de calcul de rho pour XREF et TREF via un zbrent numrec
! ..... XREF est defini par iprop 
! ..... 1:irho /  2:jtemp/  3:rho/ 4:Temp/  5:Ptot/     
! ..... 6:dP/dv|T/  7:dP/dT|v/  8:etot/     9:de/dv|T/ 10:de/dT|v
! .....11:stot / 12:ds/dv|T/ 13:ds/dT|v/ 14:Gibbs/   15:dG/dv|T/ 16:dG/dT|v/
! .....17:dv/dP|T/ 18:dv/dT|P/ 19:de/dP|T/ 20:de/dT|P/ 21:ds/dP|T/ 22:ds/dT|P/
! .....23:dG/dP|T/ 24:dG/dT|P/
! ..... G. Robert 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL*8 FUNCTION F_R_XT_ZBRENT(Tab_rT,X_ref,V_ref,r_min,r_max,iprop,SAFE,rguess,Delta_r,itype)
	USE M_CST,           ONLY : ZB_PT_ITMAX,ZB_PT_TOL,ZB_PT_EPS
	USE M_CST,           ONLY : CST_Nval_max  
	USE M_calc_table_rt, ONLY : S_calc_thermo_one_rt
	USE M_CST_TABLE_EE
	IMPLICIT NONE 
    TYPE(TABLE_RT)     :: Tab_rT
	INTEGER,OPTIONAL   :: SAFE        ! Si present et pas de solution renvoit une borne min ou max sinon -1 
	INTEGER            :: iprop       ! grandeur retournée
	REAL*8               :: X_ref,V_ref ! valeurs de references 
	REAL*8,OPTIONAL      :: r_min,r_max ! bornes min et max. Ne sont pas utilisé si rguess et Delta_r
	REAL*8,OPTIONAL      :: rguess      ! Si present, on utilise ce Vguess
	REAL*8,OPTIONAL      :: Delta_r     ! Si present, on utilise ce Delta_r pour fixer les bornes 
        INTEGER            :: itype       ! si 1 V_ref=T_ref sinon V_ref=r_ref

	INTEGER iter
	REAL*8 c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
	REAL*8 	    ::  Valeurs_a(CST_Nval_max)
	REAL*8 	    ::  Valeurs_b(CST_Nval_max)

	REAL*8 rmin(1),rmax(1)
	REAL*8 VFIXE(1)
	
	
! ..... Initialisation 
	if (PRESENT(rguess).AND.PRESENT(Delta_r)) then 
	  rmin(1)=rguess*(1.-Delta_r)
	  rmax(1)=rguess*(1.+Delta_r)
	elseif (PRESENT(r_min).AND.PRESENT(r_max)) then   
	  rmin(1)=r_min
	  rmax(1)=r_max
	else
	  write(*,'(a)') 'Il doit avoir des bornes dans F_R_XT_ZBRENT'
	  write(*,'(a)') 'r_min,r_max ou rguess,Delta_r'
	  stop
	endif
        
	
	VFIXE(1)=V_ref
	if (itype.eq.1) then 
	  call S_calc_thermo_one_rt(rho=rmin,Temp=VFIXE,Tab_rT=Tab_rT,prop='All',val_out=Valeurs_a)
	  fa=Valeurs_a(iprop)-X_ref
	  call S_calc_thermo_one_rt(rho=rmax,Temp=VFIXE,Tab_rT=Tab_rT,prop='All',val_out=Valeurs_b)
	  fb=Valeurs_b(iprop)-X_ref
	else
	  call S_calc_thermo_one_rt(rho=VFIXE,Temp=rmin,Tab_rT=Tab_rT,prop='All',val_out=Valeurs_a)
	  fa=Valeurs_a(iprop)-X_ref
	  call S_calc_thermo_one_rt(rho=VFIXE,Temp=rmax,Tab_rT=Tab_rT,prop='All',val_out=Valeurs_b)
	  fb=Valeurs_b(iprop)-X_ref
	endif
	
! ..... Pt de départ
	if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then 
		if (abs(fa).lt.abs(fb)) then 
			if (present(SAFE)) then 
				F_R_XT_ZBRENT=rmin(1)
			else
				F_R_XT_ZBRENT=-1.
			endif
		else
			if (present(SAFE)) then 
				F_R_XT_ZBRENT=rmax(1)
			else
				F_R_XT_ZBRENT=-1.
			endif
		endif
		return
	endif 
	c=rmax(1)
	fc=fb
	do iter=1,ZB_PT_ITMAX
		if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
			c=rmin(1)
			fc=fa
			d=rmax(1)-rmin(1)
			e=d
		endif
		if(abs(fc).lt.abs(fb)) then
			rmin(1)=rmax(1)
			rmax(1)=c
			c=rmin(1)
			fa=fb
			fb=fc
			fc=fa
		endif
		tol1=2.*ZB_PT_EPS*abs(rmax(1))+0.5*ZB_PT_TOL
		xm=.5*(c-rmax(1))
		if(abs(xm).le.tol1 .or. fb.eq.0.)then
			F_R_XT_ZBRENT=rmax(1)
			return
		endif
		if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
			s=fb/fa
			if(rmin(1).eq.c) then
				p=2.*xm*s
				q=1.-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2.*xm*q*(q-r)-(rmax(1)-rmin(1))*(r-1.))
				q=(q-1.)*(r-1.)*(s-1.)
			endif
			if(p.gt.0.) q=-q
			p=abs(p)
			if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
				e=d
				d=p/q
			else
				d=xm
				e=d
			endif
		else
			d=xm
			e=d
		endif
		rmin(1)=rmax(1)
		fa=fb
		if(abs(d) .gt. tol1) then
			rmax(1)=rmax(1)+d
		else
			rmax(1)=rmax(1)+sign(tol1,xm)
		endif
		if (itype.eq.1) then 
		  call S_calc_thermo_one_rt(rho=rmax,Temp=VFIXE,Tab_rT=Tab_rT,prop='All',val_out=Valeurs_b)
		  fb=Valeurs_b(iprop)-X_ref
		else
		  call S_calc_thermo_one_rt(rho=VFIXE,Temp=rmax,Tab_rT=Tab_rT,prop='All',val_out=Valeurs_b)
		  fb=Valeurs_b(iprop)-X_ref
		endif
	enddo
	F_R_XT_ZBRENT=rmax(1)
	return
END FUNCTION F_R_XT_ZBRENT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL*8 FUNCTION F_calc_ro_raw(Pstd,Tstd,Tab_rT)
! ..... Recherche pour un P et un T un guess de départ a la louche
USE M_interpol_EOS, ONLY : locate_eos
USE M_CST_TABLE_EE
IMPLICIT NONE
REAL*8	         :: Pstd,Tstd
TYPE(TABLE_RT) :: Tab_rT

INTEGER     ::  jtemi
INTEGER     ::  ii,ll

REAL*8	  ::  xtmp_r(Tab_rT%nr)

call locate_eos(Tab_rT%Temp(1:Tab_rT%nT),Tab_rT%nT,Tstd,jtemi)
if (jtemi.lt.2)	                   jtemi=2
if (jtemi.gt.Tab_rT%nT-2)                 jtemi=Tab_rT%nT-2
ll=(jtemi-1)*Tab_rT%nr
xtmp_r(1:Tab_rT%nr)=Tab_rT%Ptot(ll+1:ll+Tab_rT%nr)      
ii=MINLOC(abs(xtmp_r(:)-Pstd),DIM=1)
F_calc_ro_raw=Tab_rT%rho(ii)


END FUNCTION F_calc_ro_raw
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL*8 FUNCTION F_R_PT(P_ref,T_ref,Tab_rT)
! ..... Calcul de r pour un P et un T et une table en RT
! ..... Tab_rT est la table en rT
USE M_CST_TABLE_EE
IMPLICIT NONE 
TYPE(TABLE_RT) :: Tab_rT
REAL*8         :: P_ref,T_ref
INTEGER        :: iprop
F_R_PT=F_calc_ro_raw(P_ref,T_ref,Tab_rT)
iprop =5  !Ptot avec T=cst
F_R_PT=F_R_XT_ZBRENT(Tab_rT=Tab_rT,X_ref=P_ref,V_ref=T_ref,iprop=iprop,rguess=F_R_PT,Delta_r=0.2d0,itype=1)
END FUNCTION F_R_PT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE S_DEF_CBF(Tab_rT)
! ..... Definition de la courbe froide en rT
USE M_CST_TABLE_EE
IMPLICIT NONE 
TYPE(TABLE_RT) :: Tab_rT

! ..... Allocation mémoire 
IF (ALLOCATED(Tab_rT%cbf_etot))   DEALLOCATE(Tab_rT%cbf_etot);
ALLOCATE(Tab_rT%cbf_etot(Tab_rT%nr))
! ..... Definition de l'énergie à la temperature la plus basse 
Tab_rT%cbf_etot(1:Tab_rT%nr)=Tab_rT%etot(1:Tab_rT%nr)
END SUBROUTINE S_DEF_CBF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE S_DEF_GIBBS(Tab_rT)
! ..... Calcul de l'energie libre de Gibbs G=e-Ts+Pv
USE M_CST_TABLE_EE
IMPLICIT NONE 
TYPE(TABLE_RT) :: Tab_rT
INTEGER        :: ii,jj,kk

! ..... Allocation mémoire 
IF (ALLOCATED(Tab_rT%Gibbs))   DEALLOCATE(Tab_rT%Gibbs);
ALLOCATE(Tab_rT%Gibbs(Tab_rT%nr*Tab_rT%nT))

! ..... Calcul de l'energie de Gibbs 
ii=1
DO kk=1,Tab_rT%nT
  DO jj=1,Tab_rT%nr
    Tab_rT%Gibbs(ii)=Tab_rT%etot(ii)-Tab_rT%Temp(kk)*Tab_rT%stot(ii)&
		    +Tab_rT%Ptot(ii)/Tab_rT%rho(jj)
  ii=ii+1
  ENDDO
ENDDO
  
END SUBROUTINE S_DEF_GIBBS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE S_DEFTAB_rE(Tab_rT,Tab_re)
! ..... Construction d'une grille de correspondance de e(vT)-T(ve)
USE M_CST_TABLE_EE
USE M_CST,         ONLY : CST_shiftEint
USE M_CST,         ONLY : CST_NxTgril
USE M_CST,         ONLY : CST_deltaemin
USE M_CST,         ONLY : CST_Nval_max

USE M_interpol_EOS, ONLY : locate_eos
USE M_interpol_EOS, ONLY : sort1,sort_m
USE M_interpol_EOS, ONLY : S_INTERPOLE_RATF1_eos

USE M_calc_table_rt, ONLY : S_calc_thermo_one_rt

IMPLICIT NONE
TYPE(TABLE_RT) :: Tab_rT
TYPE(TABLE_RE) :: Tab_re


! ..... grilles temporaires  
REAL*8, ALLOCATABLE ::  grille_xtmp(:)
REAL*8, ALLOCATABLE ::  grille_etmp(:,:)
REAL*8, ALLOCATABLE ::  grille_etmpA(:,:)

! ..... Locales 
INTEGER           :: ii,jj,kk
INTEGER           :: irho0
INTEGER           :: iprop
INTEGER           :: ielocA
REAL*8              :: valeurs(CST_Nval_max)
REAL*8              :: dtmp

! ..... Creation de la grille temporaire 
IF (ALLOCATED(grille_xtmp))    DEALLOCATE(grille_xtmp);    ;
ALLOCATE(grille_xtmp(Tab_rT%nr*Tab_rT%nT))

! ..... On limite la concavite en construisant une grille e(v,T)-ecbf(v)
kk=1
do jj=1,Tab_rT%nT
do ii=1,Tab_rT%nr
  grille_xtmp(kk)=log(Tab_rT%etot(kk)-Tab_rT%cbf_etot(ii)+CST_shiftEint)  ! ATTENTION en LOG+Shift e0
  if (grille_xtmp(kk).lt.0.) grille_xtmp(kk)=0.
  kk=kk+1
enddo
enddo




! ..... On cherche les énergies min et max afin de 
! ..... pouvoir construire une grille regulière 
! ..... Pour cela, on suppose l'isochore proche de l'équilibre 
! ..... et celle de la plus forte compression 
! ..... Creation de la grille temporaire de CST_NxTgril*Tab_rT%nT
IF (ALLOCATED(grille_etmp))    DEALLOCATE(grille_etmp);   
ALLOCATE(grille_etmp(CST_NxTgril*Tab_rT%nT,1))
IF (ALLOCATED(grille_etmpA))   DEALLOCATE(grille_etmpA);   
ALLOCATE(grille_etmpA(CST_NxTgril*Tab_rT%nT,1))

! ..... localisation table r de rho0
call locate_eos(Tab_rT%rho(1:Tab_rT%nr),Tab_rT%nr,Tab_rT%Std%rho,irho0)
if (irho0.lt.2)			irho0=2
if (irho0.gt.Tab_rT%nr-2) 	irho0=Tab_rT%nr-2
! ..... grille_etmp est l'énergie sur cette isochore sans la cbf 
ii=irho0
do jj=1,Tab_rT%nT
  kk=ii+(jj-1)*Tab_rT%nr
  grille_etmp(jj,1)=grille_xtmp(kk)
enddo  
! ..... grille_etmp est l'énergie sur cette isochore sans la cbf 
ii=Tab_rT%nr
do jj=1,Tab_rT%nT
  kk=ii+(jj-1)*Tab_rT%nr
  grille_etmp(jj+Tab_rT%nT,1)=grille_xtmp(kk)
enddo
! ..... Trie en ordre croissant   
CALL sort1(CST_NxTgril*Tab_rT%nT,grille_etmp(1:CST_NxTgril*Tab_rT%nT,1))  

! on supprime les solutions doubles ie denergie > CST_deltaemin 
Tab_re%ne=1
grille_etmpA(Tab_re%ne,1)=grille_etmp(Tab_re%ne,1)
do ii=2,CST_NxTgril*Tab_rT%nT
  if ((grille_etmp(ii,1)-grille_etmp(ii-1,1)).gt.CST_deltaemin) then 
    Tab_re%ne=Tab_re%ne+1
    grille_etmpA(Tab_re%ne,1)=grille_etmp(ii,1)
  endif
enddo


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Allocation mémoire 
Tab_re%nr=Tab_rT%nr

IF (ALLOCATED(Tab_re%rho))   DEALLOCATE(Tab_re%rho);
ALLOCATE(Tab_re%rho(Tab_re%nr))
IF (ALLOCATED(Tab_re%cbf_etot))   DEALLOCATE(Tab_re%cbf_etot);
ALLOCATE(Tab_re%cbf_etot(Tab_re%nr))
IF (ALLOCATED(Tab_re%ene))   DEALLOCATE(Tab_re%ene);
ALLOCATE(Tab_re%ene(Tab_re%ne))
IF (ALLOCATED(Tab_re%Temp))   DEALLOCATE(Tab_re%Temp);
ALLOCATE(Tab_re%Temp(Tab_re%nr*Tab_re%ne))


! ..... Definition des bornes Tmin et Tmax 
Tab_re%Tmin=Tab_rT%Temp(1)
Tab_re%Tmax=Tab_rT%Temp(Tab_rT%nT)
! ..... Definition de la cbfroide 
Tab_re%cbf_etot(1:Tab_rT%nr)=Tab_rT%cbf_etot(1:Tab_rT%nr)
! ..... Definition de la grille en densite
Tab_re%rho(1:Tab_re%nr)     =Tab_rT%rho(1:Tab_rT%nr)
! ..... Definition de la grille en energie 
Tab_re%ene(1:Tab_re%ne)     =grille_etmpA(1:Tab_re%ne,1)


! ..... Reallocation grille_etmp
IF (ALLOCATED(grille_etmp))    DEALLOCATE(grille_etmp);   
ALLOCATE(grille_etmp(Tab_rT%nT,2))   ! energie et T
IF (ALLOCATED(grille_etmpA))    DEALLOCATE(grille_etmpA);   
ALLOCATE(grille_etmpA(Tab_rT%nT,2))   ! energie et T

! ..... Calcul pour toute les isochores des énergies 
! ..... Pour cela on cherche la correspondance e(T) puis on inverse et interpole via ratfson

! ..... Creation de la grille temporaire qui va contenir les températures pour les Tab_re%ene
IF (ALLOCATED(grille_xtmp))    DEALLOCATE(grille_xtmp);    ;
ALLOCATE(grille_xtmp(Tab_re%ne))
iprop=8
do ii=1,Tab_rT%nr
  grille_xtmp(:)=-1.
  ! ..... Calcul pour toutes les T de e* qui est décalé de la cbfroide 
  do jj=1,Tab_rT%nT
    CALL S_calc_thermo_one_rt(rho=Tab_rT%rho(ii),Temp=Tab_rT%Temp(jj),Tab_rT=Tab_rT,prop='etot',val_out=valeurs)
    grille_etmp(jj,1)=log(valeurs(iprop)-Tab_rT%cbf_etot(ii)+CST_shiftEint)  ! ATTENTION en LOG+Shift e0
    grille_etmp(jj,2)=Tab_rT%Temp(jj)
    if (grille_etmp(jj,1).lt.0.) grille_etmp(jj,1)=0.
  enddo
  ! ..... Suppression de doublons
  ielocA=1
  grille_etmpA(ielocA,1:2)=grille_etmp(ielocA,1:2)
  do jj=2,Tab_rT%nT
    if ((grille_etmp(jj,1)-grille_etmp(jj-1,1)).gt.CST_deltaemin) then 
      grille_etmpA(ielocA+1,1:2)=grille_etmp(jj,1:2)
      ielocA=ielocA+1
    endif
  enddo
  ! ..... Interpolation 
  DO jj=1,Tab_re%ne
     ! ..... energie entre les deux, on interpole
    if ((Tab_re%ene(jj).ge.grille_etmpA(1,1)).and.(Tab_re%ene(jj).le.grille_etmpA(ielocA,1))) then 
    CALL S_INTERPOLE_RATF1_eos(grille_etmpA(1:ielocA,1),grille_etmpA(1:ielocA,2),ielocA,&
    Tab_re%ene(jj),grille_xtmp(jj),dtmp)
     ! ..... energie < borne min --> Tmin
    elseif (Tab_re%ene(jj).lt.grille_etmpA(1,1)) then 
      grille_xtmp(jj)=grille_etmpA(1,2)   !Tmin
     ! ..... energie > borne min --> Tmax
    else
      grille_xtmp(jj)=grille_etmpA(ielocA,2)	!Tmax
    endif
     ! ..... verification 
    if (grille_xtmp(jj).gt.grille_etmpA(ielocA,2)) grille_xtmp(jj)=grille_etmpA(ielocA,2)
    if (grille_xtmp(jj).lt.grille_etmpA(1,2))      grille_xtmp(jj)=grille_etmpA(1,2)
  ENDDO  

  DO jj=1,Tab_re%ne
    kk=(jj-1)*Tab_rT%nr+ii ! indice tableau grille_ttot_re
    Tab_re%Temp(kk)=grille_xtmp(jj)
  ENDDO
enddo
! ..... On a fini la conversion 


END SUBROUTINE S_DEFTAB_RE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE S_DEFTAB_PT(Tab_rT,Tab_PT)
! ..... Construction d'une grille de correspondance de P(vT)-v(PT)
USE M_CST_TABLE_EE
USE M_CST,         ONLY : CST_shiftPtot
USE M_CST,         ONLY : CST_NxPgril
USE M_CST,         ONLY : CST_deltaPmin
USE M_CST,         ONLY : CST_Nval_max

USE M_interpol_EOS, ONLY : locate_eos
USE M_interpol_EOS, ONLY : sort1,sort_m
USE M_interpol_EOS, ONLY : S_INTERPOLE_RATF1_eos

USE M_calc_table_rt, ONLY : S_calc_thermo_one_rt

IMPLICIT NONE
TYPE(TABLE_RT) :: Tab_rT
TYPE(TABLE_PT) :: Tab_PT


! ..... grilles temporaires  
REAL*8, ALLOCATABLE ::  grille_xtmp(:)
REAL*8, ALLOCATABLE ::  grille_Ptmp(:,:)
REAL*8, ALLOCATABLE ::  grille_PtmpA(:,:)

! ..... Locales 
INTEGER           :: ii,jj,kk
INTEGER           :: irho0
INTEGER           :: iprop
INTEGER           :: ielocA
REAL*8              :: valeurs(CST_Nval_max)
REAL*8              :: dtmp

! ..... Creation de la grille temporaire 
IF (ALLOCATED(grille_xtmp))    DEALLOCATE(grille_xtmp);    ;
ALLOCATE(grille_xtmp(Tab_rT%nr*Tab_rT%nT))

! ..... On limite la concavite en construisant une grille P(v,T)-Pmin
! ..... On va supposer que dP/dv|T < 0  soit dP/dr|T > 0

Tab_PT%Pmin=MINVAL(Tab_rT%Ptot)

kk=1
do jj=1,Tab_rT%nT
do ii=1,Tab_rT%nr
  grille_xtmp(kk)=log(Tab_rT%Ptot(kk)-Tab_PT%Pmin+CST_shiftPtot)  ! ATTENTION en LOG+Shift e0
  if (grille_xtmp(kk).lt.0.) grille_xtmp(kk)=0.
  kk=kk+1
enddo
enddo



! ..... On cherche les pressions min et max afin de 
! ..... pouvoir construire une grille regulière 
! ..... Pour cela, on suppose l'isotherme min et max 
! ..... Creation de la grille temporaire de CST_NxPgril*Tab_rT%nr
IF (ALLOCATED(grille_Ptmp))    DEALLOCATE(grille_Ptmp);   
ALLOCATE(grille_Ptmp(CST_NxPgril*Tab_rT%nr,1))
IF (ALLOCATED(grille_PtmpA))   DEALLOCATE(grille_PtmpA);   
ALLOCATE(grille_PtmpA(CST_NxPgril*Tab_rT%nr,1))

! ..... grille_Ptmp est l'énergie sur cette isotherme avec décalage Tab_PT%Pmin et log
jj=1
do ii=1,Tab_rT%nr
  kk=ii+(jj-1)*Tab_rT%nr
  grille_Ptmp(ii,1)=grille_xtmp(kk)
enddo  
! ..... grille_Ptmp est l'énergie sur cette isochore sans la cbf 
jj=Tab_rT%nT
do ii=1,Tab_rT%nr
  kk=ii+(jj-1)*Tab_rT%nr
  grille_Ptmp(ii+Tab_rT%nr,1)=grille_xtmp(kk)
enddo
! ..... Trie en ordre croissant   
CALL sort1(CST_NxPgril*Tab_rT%nr,grille_Ptmp(1:CST_NxPgril*Tab_rT%nr,1))  

! on supprime les solutions doubles ie denergie > CST_deltaemin 
Tab_PT%nP=1
grille_PtmpA(Tab_PT%nP,1)=grille_Ptmp(Tab_PT%nP,1)
do ii=2,CST_NxPgril*Tab_rT%nr
  if ((grille_Ptmp(ii,1)-grille_Ptmp(ii-1,1)).gt.CST_deltaPmin) then 
    Tab_PT%nP=Tab_PT%nP+1
    grille_PtmpA(Tab_PT%nP,1)=grille_Ptmp(ii,1)
  endif
enddo


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Allocation mémoire 
Tab_PT%nT=Tab_rT%nT

IF (ALLOCATED(Tab_PT%Ptot))   DEALLOCATE(Tab_PT%Ptot);
ALLOCATE(Tab_PT%Ptot(Tab_PT%nP))
IF (ALLOCATED(Tab_PT%Temp))   DEALLOCATE(Tab_PT%Temp);
ALLOCATE(Tab_PT%Temp(Tab_PT%nT))
IF (ALLOCATED(Tab_PT%rho))   DEALLOCATE(Tab_PT%rho);
ALLOCATE(Tab_PT%rho(Tab_PT%nP*Tab_PT%nT))

! ..... Definition des bornes Tmin et Tmax 
Tab_PT%rmin=Tab_rT%rho(1)
Tab_PT%rmax=Tab_rT%rho(Tab_rT%nr)
! ..... Definition de la grille en temperature
Tab_PT%temp(1:Tab_PT%nT)     =Tab_rT%temp(1:Tab_rT%nT)
! ..... Definition de la grille en pression
Tab_PT%Ptot(1:Tab_PT%nP)     =grille_PtmpA(1:Tab_PT%nP,1)

! ..... Reallocation grille_Ptmp
IF (ALLOCATED(grille_Ptmp))    DEALLOCATE(grille_Ptmp);   
ALLOCATE(grille_Ptmp(Tab_rT%nr,2))   ! pression et volume
IF (ALLOCATED(grille_PtmpA))    DEALLOCATE(grille_PtmpA);   
ALLOCATE(grille_PtmpA(Tab_rT%nr,2))   ! pression et volume

! ..... Calcul pour toute les isotermes des pressions
! ..... Pour cela on cherche la correspondance P(v) puis on inverse et interpole via ratfson

! ..... Creation de la grille temporaire qui va contenir les volumes pour les Tab_PT%Ptot
IF (ALLOCATED(grille_xtmp))    DEALLOCATE(grille_xtmp);    ;
ALLOCATE(grille_xtmp(Tab_PT%nP))
iprop=5
do jj=1,Tab_rT%nT
  grille_xtmp(:)=-1.
  ! ..... Calcul pour toutes les T de P* qui est décalé de la constante
  do ii=1,Tab_rT%nr
    CALL S_calc_thermo_one_rt(rho=Tab_rT%rho(ii),Temp=Tab_rT%Temp(jj),Tab_rT=Tab_rT,prop='Ptot',val_out=valeurs)
    grille_Ptmp(ii,1)=log(valeurs(iprop)-Tab_PT%Pmin+CST_shiftPtot)  ! ATTENTION en LOG+Shift e0
    grille_Ptmp(ii,2)=Tab_rT%rho(ii)
    if (grille_Ptmp(ii,1).lt.0.) grille_Ptmp(ii,1)=0.
  enddo
  ! ..... Suppression de doublons
  ielocA=1
  grille_PtmpA(ielocA,1:2)=grille_Ptmp(ielocA,1:2)
  do ii=2,Tab_rT%nr
    if ((grille_Ptmp(ii,1)-grille_Ptmp(ii-1,1)).gt.CST_deltaPmin) then 
      grille_PtmpA(ielocA+1,1:2)=grille_Ptmp(ii,1:2)
      ielocA=ielocA+1
    endif
  enddo
  ! ..... Interpolation 
  DO ii=1,Tab_PT%nP
     ! ..... energie entre les deux, on interpole
    if ((Tab_PT%Ptot(ii).ge.grille_PtmpA(1,1)).and.(Tab_PT%Ptot(ii).le.grille_PtmpA(ielocA,1))) then 
    CALL S_INTERPOLE_RATF1_eos(grille_PtmpA(1:ielocA,1),grille_PtmpA(1:ielocA,2),ielocA,&
    Tab_PT%Ptot(ii),grille_xtmp(ii),dtmp)
     ! ..... energie < borne min --> Pmin
    elseif (Tab_PT%Ptot(ii).lt.grille_PtmpA(1,1)) then 
      grille_xtmp(ii)=grille_PtmpA(1,2)   !Pmin
     ! ..... energie > borne min --> Pmax
    else
      grille_xtmp(ii)=grille_PtmpA(ielocA,2)	!Pmax
    endif
     ! ..... verification 
    if (grille_xtmp(ii).gt.grille_PtmpA(ielocA,2)) grille_xtmp(ii)=grille_PtmpA(ielocA,2)
    if (grille_xtmp(ii).lt.grille_PtmpA(1,2))      grille_xtmp(ii)=grille_PtmpA(1,2)
  ENDDO  
  DO ii=1,Tab_PT%nP
    kk=(jj-1)*Tab_PT%nP+ii ! indice tableau v(P,T)
    Tab_PT%rho(kk)=grille_xtmp(ii)
  ENDDO
enddo
! ..... On a fini la conversion 
END SUBROUTINE S_DEFTAB_PT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE M_init_tools



