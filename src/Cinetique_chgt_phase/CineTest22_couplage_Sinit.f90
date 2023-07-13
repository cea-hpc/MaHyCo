MODULE SINIT
CONTAINS 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ..... Routine d'initialisation appelée par le fichier .st
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE CineTest22_couplage_Sinit(rho_std,ene_std)
USE buffer
USE M_cine_Sinit, ONLY : CineTest22_cine_Sinit
USE M_ee_Sinit,   ONLY : CineTest22_ee_Sinit
IMPLICIT NONE 
! ..... Buffers : reel
! REAL*8(KIND=8)(KIND=8)        ::  fcoeff(*)
  REAL*8        ::  rho_std,ene_std
! ..... Compteurs 
  INTEGER     ::  Nconst,ii


Nconst=1
write(*,'(2a)') ' lecture des coefficients cinetique '
! ..... Lecture des coefficients cinetique 
CALL CineTest22_cine_Sinit(fcoeff,Nconst)
write(*,'(2a)') ' lecture des coefficients EE '
! ..... Lecture des coefficients EE
CALL CineTest22_ee_Sinit(fcoeff,rho_std,ene_std,Nconst)
write(*,'(2a)') ' fin des lectures '

write(*,*) 'Valeur standard', ene_std, rho_std
!=======================================================================


END SUBROUTINE CineTest22_couplage_Sinit

END MODULE SINIT
