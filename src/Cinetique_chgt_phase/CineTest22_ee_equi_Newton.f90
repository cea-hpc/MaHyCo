MODULE M_calc_equi_Newton

CONTAINS 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Rescaling des grandeurs Xi selon le ratio X_ref/Xloc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE S_rescall(X_ref,Xloc,Xi)
	IMPLICIT NONE
	REAL*8                  :: 	X_ref,Xloc
	REAL*8                  :: 	Xi(:)
	REAL*8                  :: 	xfac
	xfac=X_ref/Xloc
	Xi(:)=Xi(:)*xfac
	RETURN
END SUBROUTINE S_rescall


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Affichage et arret si va est ISNAN ou HUGE 
! ..... G. Robert 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LOGICAL FUNCTION L_ERR_ISNAN_HUGE(va,ii,xx,yy,msg,aff)
	IMPLICIT NONE 
	REAL*8              :: va
	REAL*8              :: xx
	REAL*8              :: yy
	INTEGER           :: ii
	CHARACTER(len=*)  :: msg
	REAL*8, PARAMETER :: vamax=0.9*HUGE(va)
	INTEGER, OPTIONAL :: aff
	L_ERR_ISNAN_HUGE=.FALSE. 
	if (ISNAN(va).OR.(abs(va).GE.HUGE(vamax))) then
		if (present(aff)) write(*,'(a,i4,3e16.8)') msg,ii,xx,yy,va;
		L_ERR_ISNAN_HUGE=.TRUE.
	endif
END FUNCTION L_ERR_ISNAN_HUGE

	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Fonction de calcul de l'état P et T via le double newton
! ..... On peut implanter à moindre cout l'equilibre isentropique 
! ..... et evoluer d'un equilibre isentropique vers isotherme pour l'explo
! ..... Papier de C. Cranfill : 
! ..... EOS of a material Mixture in Pressure Equilibrium
! ..... LA13661 Jan. 2000
! ..... 
! ..... Subroutine d'appel au newton
! ..... Unité en SESAME 
! ..... ENTREE 
! ..... vol_in(nb),ene_in(nb)    : volume, energie
! ..... Fra_in(nb,nphase)        : fraction des phases 
! ..... iconv(nb)                : Flag de convergence : OK si =0, tjs <0 dans le cas contraire 
! ..... P_guess(nb)       	 : Pression initiale
! ..... T_guess(nb)       	 : Température initiale
! ..... SORTIE 
! ..... P_calc(nb)               : Pression 
! ..... T_calc(nb)               : Température 
! ..... OPTIONS 
! ..... 
! ..... G. Robert 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE S_SEARCH_NEWTON_ISOPT(nb,nip,Prop_ni,Prop_no,Fra_in,icalc)
				
	USE M_CST,              ONLY : CST_Nval_max
	USE M_CST_TABLE_EE,     ONLY : TB_rT,TB_PT,TB_re

	USE M_CST_TABLE_EE,     ONLY : NPAR_CONV_NEWTON
	USE M_CST_TABLE_EE,     ONLY : LIMIT_TABLE

	USE M_calc_table_re, ONLY : S_calc_thermo_one_ve
	USE M_calc_table_pt, ONLY : S_calc_thermo_one_PT
	
	USE M_CST_TABLE_EE,     ONLY : Thermo_Newton ! type pour les propriétés thermo
	USE M_CST_TABLE_EE,     ONLY : TAB_APPEL_EE

	USE M_CST_TABLE_EE,  ONLY : nb_restart,restart_PT

	IMPLICIT NONE 
! ..... Entrée/Sortie
	INTEGER   :: nb
	INTEGER   :: nip
	REAL*8      :: Fra_in(nb,nip)          ! taille nphases, nb

	TYPE(TAB_APPEL_EE)  :: Prop_ni(nb)
	TYPE(Thermo_Newton) :: Prop_no(nb)
	INTEGER             :: icalc(nb)   ! si 0 on calcul rien 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Local
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	REAL*8                :: rho_in(nb),ene_in(nb)  ! taille nb 
	REAL*8                :: P_guess(nb)           ! taille nb 
	REAL*8                :: T_guess(nb)           ! taille nb 

	REAL*8                :: vol_ref,ene_ref
	INTEGER             :: ii,jj,kk
! ..... Tableau de valeurs 
	REAL*8  		    :: valeurs(CST_Nval_max)
	INTEGER,PARAMETER   :: i_rho=3 ! emplacement de la densité dans le tableau valeurs
	INTEGER,PARAMETER   :: i_ene=8 ! emplacement de l'énergie dans le tableau valeurs
	REAL*8                :: vloc,eloc,xfac

! ..... Grandeurs locales Newton 
	REAL*8                :: vi_j(nip),ei_j(nip),Fra_j(nip)
	REAL*8                :: Pre_j(nip),Temp_j(nip)
	REAL*8                :: dvdP_T_j(nip),dvdT_P_j(nip)
	REAL*8                :: dedP_T_j(nip),dedT_P_j(nip)
	REAL*8                :: dedT_v_j(nip),dPdT_v_j(nip),dedv_T_j(nip),dPdv_T_j(nip)
	REAL*8                :: Gibbs_j(nip)
	REAL*8                :: Av,Bv,Cv,Ae,Be,Ce
! ..... Grandeurs locales Convergence 
	REAL*8                :: Pold(1),Told(1)
	REAL*8                :: ploc(1),tloc(1)
	REAL*8                :: dprec_P(1),dprec_T(1)
	REAL*8                :: dprec_v(1),dprec_e(1)

	REAL*8                :: dv_loc(nip)
	REAL*8                :: de_loc(nip)
	REAL*8                :: mu_v(nip),mu_e(nip)

	LOGICAL             :: icalc_PT

	INTEGER             :: NX_XINIT

	rho_in(1:nb) =Prop_ni(1:nb)%rho
	ene_in(1:nb) =Prop_ni(1:nb)%ene
	P_guess(1:nb)=Prop_ni(1:nb)%Ptot
	T_guess(1:nb)=Prop_ni(1:nb)%Temp


	do ii=1,nb
	  !if (ii.eq.1) write(*,*) " dans le newton" 
          !if (ii.eq.1) write(*,*)  'densite' , rho_in(ii) 
          !if (ii.eq.1) write(*,*)  'energie ', ene_in(ii)
          !if (ii.eq.1) write(*,*)  'pression ' , P_guess(ii)
          !if (ii.eq.1) write(*,*)  'temperature ', T_guess(ii)
	 
	  if (icalc(ii).ne.0) then 
	    Prop_no(ii)%iconv=1
	    ! ..... Affectation d'une valeur par défaut au cas de non convergence 	
	    pold(1)=P_guess(ii)
	    Told(1)=T_guess(ii)
	    ! ..... Affectation d'une valeur débile au cas de non convergence 	
	    dprec_P(1)=HUGE(dprec_P)
	    dprec_T(1)=HUGE(dprec_T)

	    ! on affecte les fractions 
	    do jj=1,nip
	      Fra_j(jj)=Fra_in(ii,jj)
	    enddo
	    icalc_PT=.TRUE.
	    IF (MINVAL(Prop_ni(ii)%v_j(1:nip)).GT.NPAR_CONV_NEWTON%VMIN) icalc_PT=.FALSE.
	    vol_ref=1./rho_in(ii)
	    ene_ref=ene_in(ii)
	    
	    IF (icalc_PT) then 
		! appel au PT pour trouver les vi et ei 
		! ces vi ei sont rescalées et utilisé comme point de départ 
		do jj=1,nip
		    CALL S_calc_thermo_one_PT(Pre=pold(1),Temp=Told(1),Tab_PT=TB_PT(jj),Tab_rT=TB_rT(jj),prop='all',val_out=valeurs)
		    vi_j(jj)=1./valeurs(i_rho)
		    ei_j(jj)=valeurs(i_ene)
		enddo
	    else
	      vi_j(1:nip)=Prop_ni(ii)%v_j(1:nip)
	      ei_j(1:nip)=Prop_ni(ii)%e_j(1:nip)
	    ENDIF
	    vloc=0.;eloc=0.;
    	    do jj=1,nip
		vloc=vloc+Fra_j(jj)*vi_j(jj)
		eloc=eloc+Fra_j(jj)*ei_j(jj)
	    enddo
	    ! conversion en volume massique 
	    ! on rescale les volumes et les energies. 
	    ! La méthode proposée par Cranfill suppose un volume et une énergie constante
	    CALL S_rescall(vol_ref,vloc,vi_j)
	    CALL S_rescall(ene_ref,eloc,ei_j)
	    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    ! ..... Debut du newton
	    ! ..... 
	    ! ..... Ce Newton assure en principe l'équilibre PT pour un vmoyen et emoyen
	    ! ..... les sommes des volumes et energies partielles vol_nwt et ein_nwt 
	    ! ..... ne doivent pas être modifiées
	    ! ..... Utilise l'approche de C. W. Cranfill : EOS of a material mixture 
	    ! ..... in pressure equilibrium LA 13661 de 2000 
	    ! ..... G. Robert 
	    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    kk=1
	    NX_XINIT=1
	    do while ((kk.lt.NPAR_CONV_NEWTON%ITMAX).AND.((dprec_P(1).ge.NPAR_CONV_NEWTON%dP).OR.(dprec_T(1).ge.NPAR_CONV_NEWTON%dT)))
		if (Prop_no(ii)%iconv.LT.0) then
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		! ..... 	Le newton conduit à une valeur de T<NPAR_PHI_TMIN soit <0K ou un NAN
		! .....     	On recherche les v_p et e_p (p pour partiel) pour pold et told 
		! .....     	puis on renormalise pour garder le vmoyen et emoyen 
		! .....     	hypothèse de Newton pour méthode de cranfill
		! .....     	et on repart de ces points 
		! .....		Une correction pourrait être 
		! .....		de partir de calculer v(ijwl,Pold,Told) et pour chacun de ces cas avoir le vguess de départ
		! .....		puis mettre dprec_P a huge pour avoir un nouveau point de départ  
		! .....		C'est une méthode 'force brute' mais la plus stable que j'ai trouvé 
		! ..... 	G. Robert 
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    ! Il y a un probleme dans la convergence
		    ! On repart d'un autre point de départ 
		    Prop_no(ii)%iconv=1
		    ! on affecte les grandeurs définies dans le module de 
		    ! variable et les v et e initiaux 
		    Pold(1)=restart_PT(NX_XINIT,1)
		    Told(1)=restart_PT(NX_XINIT,2)
		    do jj=1,nip
			CALL S_calc_thermo_one_PT(Pre=pold(1),Temp=Told(1),Tab_PT=TB_PT(jj),Tab_rT=TB_rT(jj),prop='all',val_out=valeurs)
			vi_j(jj)=1./valeurs(i_rho)
			ei_j(jj)=valeurs(i_ene)
		    enddo
		    vloc=0.;eloc=0.;
		    do jj=1,nip
			vloc=vloc+Fra_j(jj)*vi_j(jj)
			eloc=eloc+Fra_j(jj)*ei_j(jj)
		    enddo
		    ! conversion en volume massique 
		    ! on rescale les volumes et les energies. 
		    ! La méthode proposée par Cranfill suppose un volume et une énergie constante
		    CALL S_rescall(vol_ref,vloc,vi_j)
		    CALL S_rescall(ene_ref,eloc,ei_j)
		    dprec_P(1)=HUGE(dprec_P)
		    dprec_T(1)=HUGE(dprec_T)
		    ! Mise a jour du compteur d'essai en PT
		    NX_XINIT=NX_XINIT+1
		    ! Mise a jour du compteur du newton
		    kk=0
		    if (NX_XINIT.GT.nb_restart) then 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... 		On ne reussi pas a converger il faut choisir plus de 
! ..... 		point de départ (dans NPAR_PHI_PINIT et NPAR_PHI_TINIT)
! ..... 		ou changer de méthode 
! ..... 		G. Robert 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Prop_no(ii)%iconv=-60
			kk=NPAR_CONV_NEWTON%ITMAX+10
		    endif
		endif
		if (kk.gt.NPAR_CONV_NEWTON%ITMAX) exit
		do jj=1,nip
		! ......appel de l'interpolateur en ve pour les nip EE
		! ..... prop='irho'    : calcul de la position en rho       -> val_out(1)
		! ..... prop='jtemp'   : calcul de la position en temp      -> val_out(2)
		! ..... prop='rho'     : calcul de rho                      -> val_out(3)
		! ..... prop='Temp'    : calcul de Temp                     -> val_out(4)
		! ..... prop='Ptot'    : calcul de Ptot, dP/dv|T et dP/dT|v -> val_out(5:7)
		! ..... prop='etot'    : calcul de etot, de/dv|T et de/dT|v -> val_out(8:10)
		! ..... prop='stot'    : calcul de stot, ds/dv|T et ds/dT|v -> val_out(11:13)
		! ..... prop='Gibbs'   : calcul de G, dG/dv|T et dG/dT|v    -> val_out(14:16)
		! ..... autre          : + dv/dP|T,dv/dT|P,de/dP|T,de/dT|P -> val_out(17:20)
		! ..... autre          : + ds/dP|T,ds/dT|P,dG/dP|T,dG/dT|P -> val_out(21:24)
		    CALL S_calc_thermo_one_ve(vol=vi_j(jj),ene=ei_j(jj),Tab_rE=TB_re(jj),Tab_rT=TB_rT(jj),prop='all',val_out=valeurs)
		! if (ii.eq.1) write(*,*) " valeur ", valeurs
		    Temp_j(jj) =valeurs(4)
		    Pre_j(jj)  =valeurs(5)
		    dPdv_T_j(jj)=valeurs(6)
		    dPdT_v_j(jj)=valeurs(7)

		    dedv_T_j(jj)=valeurs(9)
		    dedT_v_j(jj)=valeurs(10)

		    Gibbs_j(jj)=valeurs(14)

		    dvdP_T_j(jj)=valeurs(17)
		    dvdT_P_j(jj)=valeurs(18)
		    dedP_T_j(jj)=valeurs(19)
		    dedT_P_j(jj)=valeurs(20)


		    IF (L_ERR_ISNAN_HUGE(dvdP_T_j(jj),ii=jj,xx=vol_ref,yy=ene_ref,msg='Erreur ISNAN dvdP_T ijwl v e valeur')) THEN
		      Prop_no(ii)%iconv=-21
		    ENDIF
		    IF (L_ERR_ISNAN_HUGE(dedP_T_j(jj),ii=jj,xx=vol_ref,yy=ene_ref,msg='Erreur ISNAN dedP_T ijwl v e valeur')) THEN
		      Prop_no(ii)%iconv=-22
		    ENDIF
		    IF (L_ERR_ISNAN_HUGE(dvdT_P_j(jj),ii=jj,xx=vol_ref,yy=ene_ref,msg='Erreur ISNAN dvdT_P ijwl v e valeur')) THEN
		      Prop_no(ii)%iconv=-23
		    ENDIF
		    IF (L_ERR_ISNAN_HUGE(dedT_P_j(jj),ii=jj,xx=vol_ref,yy=ene_ref,msg='Erreur ISNAN dedT_P ijwl v e valeur')) THEN
		      Prop_no(ii)%iconv=-24
		    ENDIF
		enddo
    ! ......	Calcul des constantes pour le newton 
		Av=0.;Bv=0.;Cv=0.;Ae=0.;Be=0.;Ce=0.;
		do jj=1,nip
			Av=Av+Fra_j(jj)*dvdP_T_j(jj)
			Bv=Bv+Fra_j(jj)*dvdT_P_j(jj)
			Cv=Cv+Fra_j(jj)*(Pre_j(jj)*dvdP_T_j(jj)+Temp_j(jj)*dvdT_P_j(jj))
			Ae=Ae+Fra_j(jj)*dedP_T_j(jj)
			Be=Be+Fra_j(jj)*dedT_P_j(jj)
			Ce=Ce+Fra_j(jj)*(Pre_j(jj)*dedP_T_j(jj)+Temp_j(jj)*dedT_P_j(jj))
		enddo


		ploc(1)=(Cv*Be-Bv*Ce)/(Av*Be-Bv*Ae)
		tloc(1)=(Av*Ce-Cv*Ae)/(Av*Be-Bv*Ae)
		!if (ii.eq.1) write(*,*) " temp et Pres ", Temp_j(1), "  ", Pre_j(1) 
		!if (ii.eq.1) write(*,*) " temp et Pres ", Temp_j(2), "  ", Pre_j(2) 
		!if (ii.eq.1) write(*,*) " temp et Pres ", Temp_j(3), "  ", Pre_j(3) 
		!if (ii.eq.1) write(*,*) " iteration ", kk
		!if (ii.eq.1) write(*,*) " coeffs", Av, " ", Bv, " ", Cv, " "   
		!if (ii.eq.1) write(*,*) " coeffs", Ae, " ", Be, " ", Ce, " "   
		
		dprec_P(1)=abs(ploc(1)-pold(1))
		dprec_T(1)=abs(Tloc(1)-Told(1))
		!if (ii.eq.1) write(*,*) " dprec " , dprec_P(1) , " " , dprec_T(1)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! ..... 	Quelques vérifications 
    ! ..... 	si il y a eu un prb iconv<0 et on repart d'un point en PT 
    ! ..... 	si T<Tmin ou T>Tmax
    ! ..... 	G. Robert 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if (L_ERR_ISNAN_HUGE(tloc(1),ii,tloc(1),ploc(1),msg='Erreur T stop'))	Prop_no(ii)%iconv=-30
		if (L_ERR_ISNAN_HUGE(ploc(1),ii,tloc(1),ploc(1),msg='Erreur P stop'))	Prop_no(ii)%iconv=-31
		if (tloc(1).lt.NPAR_CONV_NEWTON%TMIN) Prop_no(ii)%iconv=-41

		! on ne calcul que si Prop_no(ii)%iconv>=0
		if (Prop_no(ii)%iconv.GE.0) then
		    pold(1)=ploc(1)
		    Told(1)=Tloc(1)
! ......		! La convergence n'est pas atteinte 
			if ((dprec_P(1).ge.NPAR_CONV_NEWTON%dP).OR.(dprec_T(1).ge.NPAR_CONV_NEWTON%dT)) then 
! ......			calcul du dvol et dene de chaque constituant 
			    do jj=1,nip
				dv_loc(jj)=((ploc(1)-pre_j(jj))*dvdP_T_j(jj)+(tloc(1)-temp_j(jj))*dvdT_P_j(jj))
				de_loc(jj)=((ploc(1)-pre_j(jj))*dedP_T_j(jj)+(tloc(1)-temp_j(jj))*dedT_P_j(jj))
			    enddo    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ......			mu_v et mu_e sont des limiteurs du newton 
! ......			Dans les faits, ces limiteurs sont relativement 
! ......			peux efficaces et il vaut mieux changer le point de départ des 
! ......			volumes et énergies partielles via PT par exemple  
! ..... 			! G. Robert 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			    mu_v(:)=1.
			    mu_e(:)=1.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ......			Calcul des nouvelles fractions et du volume et energie moyennes
			    do jj=1,nip
! ......			application du delta
				vi_j(jj)=vi_j(jj)+dv_loc(jj)*mu_v(jj)
				ei_j(jj)=ei_j(jj)+de_loc(jj)*mu_e(jj)
! ......			application de limites physiques en V 
				if (vi_j(jj).lt.LIMIT_TABLE(jj)%VMIN) vi_j(jj)=LIMIT_TABLE(jj)%VMIN ! volume minimal sur chaque fraction 
				if (vi_j(jj).gt.LIMIT_TABLE(jj)%VMAX) vi_j(jj)=LIMIT_TABLE(jj)%VMAX ! volume minimal sur chaque fraction 
			    enddo
			    ! conservation de v et e 
			    vloc=0.;eloc=0.;
			    do jj=1,nip
				vloc=vloc+Fra_j(jj)*vi_j(jj)
				eloc=eloc+Fra_j(jj)*ei_j(jj)
			    enddo
			    ! on rescale les volumes et les energies. 
			    ! La méthode proposée par Cranfill suppose un volume et une énergie constante
			    CALL S_rescall(vol_ref,vloc,vi_j)
			    CALL S_rescall(ene_ref,eloc,ei_j)
			    dprec_v(1)=abs(vol_ref-vloc) ! /vol_ref     ! erreur relative modif 
			    dprec_e(1)=abs(ene_ref-eloc) ! /ene_ref      ! erreur relative modif 
			endif
		endif
		kk=kk+1
	    enddo
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! ..... Flag de non convergence tjs <0 
    ! ..... G. Robert 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    if (Prop_no(ii)%iconv.GE.0) then
		dprec_v(1)=abs(vol_ref-vloc)  ! /vol_ref   ! erreur relative modif 
		dprec_e(1)=abs(ene_ref-eloc)  ! /ene_ref   ! erreur relative modif

		if (dprec_v(1).GT.NPAR_CONV_NEWTON%dv) then
			Prop_no(ii)%iconv=-5
		elseif(dprec_e(1).GT.NPAR_CONV_NEWTON%de) then
			Prop_no(ii)%iconv=-6
			write(*,*) "Probleme de convergence : Maille " , ii 
			write(*,*) dprec_e(1) , " plus grand que ", NPAR_CONV_NEWTON%de
			Prop_no(ii)%iconv=0
		elseif (tloc(1).LT.NPAR_CONV_NEWTON%TMIN) then
			Prop_no(ii)%iconv=-7
		elseif (isnan(tloc(1))) then
			Prop_no(ii)%iconv=-8
		elseif (kk.ge.NPAR_CONV_NEWTON%ITMAX) then 
			Prop_no(ii)%iconv=-4  ! on dit que l'on a convergé 0
		else
			Prop_no(ii)%iconv=0
		endif

		Prop_no(ii)%Temp=Tloc(1)
		Prop_no(ii)%Ptot=Ploc(1)

		!if (ii.eq.1) then
		!write(*,*) " iteration ", kk
		!write(*,*) ii, 'Temperature locale ', Tloc(ii), 'Pression locale ', Ploc(ii)
		!write(*,*) ii, 'Volume massique locale ', vloc, 'Energie locale ', eloc
		!endif
		
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! ..... sorties additionnelles
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		Prop_no(ii)%rho=1./vloc
		Prop_no(ii)%ene=eloc
		Prop_no(ii)%v_j(1:nip)=vi_j(1:nip)
		Prop_no(ii)%e_j(1:nip)=ei_j(1:nip)

		Prop_no(ii)%dedP_T_j(1:nip)=dedP_T_j(1:nip)
		Prop_no(ii)%dvdP_T_j(1:nip)=dvdP_T_j(1:nip)
		Prop_no(ii)%dedT_P_j(1:nip)=dedT_P_j(1:nip)
		Prop_no(ii)%dvdT_P_j(1:nip)=dvdT_P_j(1:nip)
		
		Prop_no(ii)%Gibbs(1:nip)=Gibbs_j(1:nip)
	    else
	      !! il y a un probleme on renvoit la flag tel que 
	    endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else 
	  !deja calculé on fait rien 
	    
	  endif ! if (icalc(ii).ne.0) then  
	enddo  ! do ii=1,nb
    !!write(*,*) ' NEWTON iter ', kk, 'Pression finale ', Prop_no(1)%Ptot , ' Temperature finale ', Prop_no(1)%Temp
    RETURN

END SUBROUTINE S_SEARCH_NEWTON_ISOPT


END MODULE M_calc_equi_Newton

