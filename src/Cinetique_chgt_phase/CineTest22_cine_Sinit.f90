MODULE M_cine_Sinit
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Routine d'initialisation appelée par le fichier .st
! ..... Module à supprimer et NCONST à mettre à 1 sans passage si modèle
! ..... composite
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE CineTest22_cine_Sinit(fcoeff,Nconst)
  USE M_CST_TABLE_CINE, ONLY : NPAR_Cine,nip_max_cine
  USE M_CST_TABLE_CINE, ONLY : Trans_max_up,Trans_max_dn
  USE M_CST_TABLE_CINE, ONLY : nip_cine
  USE M_CST_TABLE_EE,   ONLY : nip_ee

  IMPLICIT NONE


! ..... Buffers : un reel
  REAL(KIND=8)        ::  fcoeff(*)
  INTEGER     ::  Nconst
! ..... Compteurs 
  INTEGER     :: ii,jj

  write(*,*) 'Rentré dans CineTest22_cine_Sinit'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Initialisation par les jeux de coefficients inclus dans le buffer
! ..... Nb de phase 
  nip_cine=int(fcoeff(Nconst));Nconst=Nconst+1
  write(*,*) 'nip_cine', nip_cine

! ..... Verif nb phases max 
  IF (nip_cine.gt.nip_max_cine) THEN
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
  
! ..... Lecture du taux de transition maximal 
  NPAR_Cine%trans_orig=fcoeff(Nconst)
  NPAR_Cine%trans_orig=MAX(MIN(NPAR_Cine%trans_orig,Trans_max_up),Trans_max_dn)
  Nconst=Nconst+1

! ..... Lecture du type de cinetique
! ..... par defaut cinetique de Hayes
  NPAR_Cine%type_cinetique=0
  NPAR_Cine%type_cinetique=int(fcoeff(Nconst))
  write(*,*) 'NPAR_Cine%type_cinetique =',  NPAR_Cine%type_cinetique , "  ", fcoeff(Nconst)
  Nconst=Nconst+1

! ..... Decalage des energies libres Gii
  NPAR_Cine%shift_Gii(:)=0
  DO ii=1,nip_cine
     NPAR_Cine%shift_Gii(ii)=fcoeff(Nconst)
     Nconst=Nconst+1
  ENDDO
  

! ..... Lecture des propriétés de la cinétique 
  DO ii=1,nip_cine
    DO jj=1,nip_cine
      NPAR_Cine%B_ij(ii,jj)=fcoeff(Nconst)
      Nconst=Nconst+1
    ENDDO
  ENDDO
  DO ii=1,nip_cine
    DO jj=1,nip_cine
      NPAR_Cine%nu_ij(ii,jj)=fcoeff(Nconst)
      Nconst=Nconst+1
    ENDDO
  ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Les paramètres sont lus  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


END SUBROUTINE CineTest22_cine_Sinit

END MODULE M_cine_Sinit
