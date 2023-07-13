MODULE M_calc_equi_PT

CONTAINS 


SUBROUTINE S_CHECK_FRACTION(nip,F_in)
	IMPLICIT NONE 
	INTEGER              :: nip
	REAL*8                 :: F_in(nip)
	REAL*8                 :: F_out(nip)
	INTEGER              :: jj
	REAL                 :: norm
	F_out(1:nip)=F_in(1:nip)
	norm=0.
	do jj=1,nip
	    if (F_out(jj).lt.0.) F_out(jj)=0.
	    if (F_out(jj).gt.1.) F_out(jj)=1.
		norm=norm+F_out(jj)
	enddo
	F_in(1:nip)=F_out(1:nip)/norm
	RETURN 
END SUBROUTINE S_CHECK_FRACTION



SUBROUTINE S_CALC_CINE_PT(nbmail,nip,Prop_pas_ni,Prop_pas_nf,Frac_ni,icalc)

USE M_calc_equi_Newton, ONLY : S_SEARCH_NEWTON_ISOPT

USE M_CST_TABLE_EE,     ONLY : Thermo_Newton ! type pour les propriétés thermo de sortie 
USE M_CST_TABLE_EE,     ONLY : TAB_APPEL_EE  ! type pour les propriétés thermo d'appel 



USE M_CST,              ONLY : CST_INIT_P,CST_INIT_T
USE M_CST,              ONLY : CST_INIT_P0,CST_INIT_T0

IMPLICIT NONE 



INTEGER       :: nbmail
INTEGER       :: nip
REAL*8                :: Frac_ni(nbmail,nip) ! Fraction au pas n0
TYPE(TAB_APPEL_EE)  :: Prop_pas_ni(nbmail) ! données entrées rho,ene,P,T
TYPE(Thermo_Newton) :: Prop_pas_nf(nbmail) ! données de sorties P,T, v_j etc etc 
INTEGER             :: icalc(nbmail)       ! si 0 on calcul rien 
INTEGER        :: ii


! renormalisation des fractions 
DO ii=1,nbmail
  CALL S_CHECK_FRACTION(nip,Frac_ni(ii,1:nip))
ENDDO
! initialisation des P_guess et T_guess
DO ii=1,nbmail
  if (Prop_pas_ni(ii)%Temp.lt.CST_INIT_T) then 
  ! initialisation a voir si on prend une valeur de la table .coeff
    Prop_pas_ni(ii)%Ptot=CST_INIT_P0
    Prop_pas_ni(ii)%Temp=CST_INIT_T0
  endif
ENDDO

CALL S_SEARCH_NEWTON_ISOPT(nb=nbmail,nip=nip,&
                           Prop_ni=Prop_pas_ni(1:nbmail),&
	  	           Prop_no=Prop_pas_nf(1:nbmail),&
			   Fra_in=Frac_ni(1:nbmail,1:nip),&
			   icalc=icalc(1:nbmail)&
			   )

RETURN 
END SUBROUTINE S_CALC_CINE_PT


END MODULE M_calc_equi_PT
