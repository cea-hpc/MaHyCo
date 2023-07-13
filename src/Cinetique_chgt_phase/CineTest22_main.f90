MODULE CHEMIN
CONTAINS


subroutine S_hugo_one(r0,e0,rh,eh,dt,Pres,Temp,Frac,ene,Frin)

USE  SCALC

implicit none
INTEGER, PARAMETER :: nip_max=6
INTEGER, PARAMETER :: nbmail=1
INTEGER, PARAMETER :: ii=1

REAL*8               :: dt

REAL*8  :: rho(nbmail)
REAL*8  :: ene(nbmail)
REAL*8  :: dtime(nbmail)

REAL*8  :: Frin(nbmail,nip_max)


REAL*8  :: Pres_i(nbmail)
REAL*8  :: Temp_i(nbmail)
REAL*8  :: Frac_i(nbmail,nip_max)

REAL*8  :: Pres(nbmail)
REAL*8  :: Temp(nbmail)
REAL*8  :: Frac(nbmail,nip_max)
REAL*8  :: dpde(nbmail)
REAL*8  :: cs2(nbmail)
REAL*8  :: conv(nbmail)

REAL*8 r0,e0,rh
REAL*8 DHug,dV
REAL*8 de,eh
integer kkl

Temp(:)=0.
Pres(:)=0.
Frac(:,:)=0.

Temp_i(:)=0.
Pres_i(:)=0.
Frac_i(:,:)=0.


Frac_i(ii,:)=Frin(ii,:)


rho(ii)=r0
ene(ii)=e0

dV=(1./r0-1./rh)*0.5

dtime(:)=dt
! calcul du pole 
CALL S_CALC_CINE_VE(nbmail,dtime(1:nbmail),rho(1:nbmail),ene(1:nbmail),&
				  Pres_i(1:nbmail),Temp_i(1:nbmail),&
				 Frac_i(1:nbmail,1),Frac_i(1:nbmail,2),Frac_i(1:nbmail,3),&
				 Frac_i(1:nbmail,4),Frac_i(1:nbmail,5),Frac_i(1:nbmail,6),&
				 dpde(1:nbmail),cs2(1:nbmail),conv(1:nbmail))

de=0.01

ene(ii)=eh-de
do while (de.gt.1.e-5)
kkl=0
DHug=5.
do while (DHug.gt.0.) 
if (kkl.gt.20) then
stop
endif


Pres(1:nbmail)=Pres_i(1:nbmail)
Temp(1:nbmail)=Temp_i(1:nbmail)
Frac(1:nbmail,1)=Frac_i(1:nbmail,1)
Frac(1:nbmail,2)=Frac_i(1:nbmail,2)
Frac(1:nbmail,3)=Frac_i(1:nbmail,3)
Frac(1:nbmail,4)=Frac_i(1:nbmail,4)
Frac(1:nbmail,5)=Frac_i(1:nbmail,5)
Frac(1:nbmail,6)=Frac_i(1:nbmail,6)
rho(ii)=rh
ene(ii)=ene(ii)+de


CALL S_CALC_CINE_VE(nbmail,dtime(1:nbmail),rho(1:nbmail),ene(1:nbmail),&
				  Pres(1:nbmail),Temp(1:nbmail),&
				 Frac(1:nbmail,1),Frac(1:nbmail,2),Frac(1:nbmail,3),&
				 Frac(1:nbmail,4),Frac(1:nbmail,5),Frac(1:nbmail,6),&
				 dpde(1:nbmail),cs2(1:nbmail),conv(1:nbmail))




DHug=dV*Pres(ii)-(ene(ii)-e0)
if (conv(1).ne.0.) DHug=5.
kkl=kkl+1
enddo
ene(ii)=ene(ii)-de
de=de/10.
enddo


RETURN
end subroutine S_hugo_one


subroutine S_hugo(r0,e0,rh,dt,nr,frin,frou,eh,ph,th)

INTEGER, PARAMETER :: nip_max=6
INTEGER, PARAMETER :: nbmail=1
INTEGER, PARAMETER :: ii=1

REAL*8  :: Pres(nbmail)
REAL*8  :: Temp(nbmail)
REAL*8  :: Frac(nbmail,nip_max)
REAL*8  :: ene(nbmail)

REAL*8  :: Frin(nbmail,nip_max)
REAL*8  :: Frou(nbmail,nip_max)

REAL*8    r0,e0,rh,dt
integer nr,jj
REAL*8 drho,rhi,eh,ph,th

drho=(rh-r0)/(1.*nr)

eh=e0

do jj=1,nr
rhi=r0+jj*drho


call S_hugo_one(r0,e0,rhi,eh,dt,Pres,Temp,Frac,ene,frin)
eh=ene(ii)
!    write(501,'(20d)') rhi,ene(ii)
!    write(601,'(20d)') Pres(ii),Frac(ii,1)
!    write(602,'(20d)') Pres(ii),Frac(ii,2)
!    write(603,'(20d)') Pres(ii),Frac(ii,3)
!    write(701,'(20d)') Temp(ii),Frac(ii,1)
!    write(702,'(20d)') Temp(ii),Frac(ii,2)
!    write(703,'(20d)') Temp(ii),Frac(ii,3)
!    write(801,'(20d)') Pres(ii),Temp(ii)
!    write(802,'(20d)') Pres(ii),rhi
!    write(803,'(20d)') rhi,Temp(ii)
ph=Pres(ii)
th=Temp(ii)
enddo
    write(501,'(a)') ''
    write(601,'(a)') ''
    write(602,'(a)') ''
    write(603,'(a)') ''
    write(701,'(a)') ''
    write(702,'(a)') ''
    write(703,'(a)') ''
    write(801,'(a)') ''
    write(802,'(a)') ''
    write(803,'(a)') ''
Frou(ii,:)=Frac(ii,:)
end subroutine S_hugo







subroutine S_isentrope(r0,e0,dr,rfin,dt,frin,frou,pin,tin,isens)
USE  SCALC
implicit none

INTEGER, PARAMETER :: nip_max=6
INTEGER, PARAMETER :: nbmail=1
INTEGER, PARAMETER :: ii=1

REAL*8  :: Frin(nbmail,nip_max)
REAL*8  :: Frou(nbmail,nip_max)
REAL*8  :: Pin,Tin

REAL*8  :: rho(nbmail)
REAL*8  :: ene(nbmail)
REAL*8  :: Pres(nbmail)
REAL*8  :: Temp(nbmail)
REAL*8  :: Frac(nbmail,nip_max)
REAL*8  :: dtime(nbmail)
REAL*8  :: dpde(nbmail)
REAL*8  :: cs2(nbmail)
REAL*8  :: conv(nbmail)
integer iloop
real r0,e0,dr,rfin,dt,drho
INTEGER :: isens

Temp(ii)=tin
Pres(ii)=pin
Frac(ii,:)=Frin(ii,:)
iloop=1

rho(ii)=r0
ene(ii)=e0
dtime(:)=dt
drho=dr
do while(iloop.gt.0)
iloop=iloop+1

CALL S_CALC_CINE_VE(nbmail,dtime(1:nbmail),rho(1:nbmail),ene(1:nbmail),&
				  Pres(1:nbmail),Temp(1:nbmail),&
				 Frac(1:nbmail,1),Frac(1:nbmail,2),Frac(1:nbmail,3),&
				 Frac(1:nbmail,4),Frac(1:nbmail,5),Frac(1:nbmail,6),&
				 dpde(1:nbmail),cs2(1:nbmail),conv(1:nbmail))

print*, Frac(1,:)
print*, rho,ene,Pres,Temp

if (conv(ii).ne.0) then
print*, 'iiii',conv(ii)

endif

if (iloop.eq.2000) stop

if (isens.gt.0) then 
  ene(ii)=ene(ii)-Pres(ii)*(1./(rho(ii)+drho)-1./rho(ii))
  rho(ii)=rho(ii)+drho
  if (rho(ii).gt.rfin) exit
else
  ene(ii)=ene(ii)+Pres(ii)*(1./(rho(ii)+drho)-1./rho(ii))
  rho(ii)=rho(ii)-drho
  if (rho(ii).lt.rfin) exit
endif

!    write(501,'(20d)') rho(ii),ene(ii)
!    write(502,'(20d)') rho(ii),ene(ii)
!    write(601,'(20d)') Pres(ii),Frac(ii,1)
!    write(602,'(20d)') Pres(ii),Frac(ii,2)
!    write(603,'(20d)') Pres(ii),Frac(ii,3)
!    write(701,'(20d)') Frac(ii,1),Temp(ii)
!    write(702,'(20d)') Frac(ii,2),Temp(ii)
!    write(703,'(20d)') Frac(ii,3),Temp(ii)
!    write(801,'(20d)') Pres(ii),Temp(ii)
!    write(802,'(20d)') rho(ii),Pres(ii)
!    write(803,'(20d)') rho(ii),Temp(ii)


enddo



    write(501,'(a)') ''
    write(601,'(a)') ''
    write(602,'(a)') ''
    write(603,'(a)') ''
    write(701,'(a)') ''
    write(702,'(a)') ''
    write(703,'(a)') ''
    write(801,'(a)') ''
    write(802,'(a)') ''
    write(803,'(a)') ''


Frou(ii,:)=Frac(ii,:)

end subroutine S_isentrope





subroutine S_isochore(r0,e0,de,efin,dt,Frin,Frou)

USE  SCALC
implicit none
INTEGER, PARAMETER :: nip_max=6
INTEGER, PARAMETER :: nbmail=1
INTEGER, PARAMETER :: ii=1

REAL*8  :: Frin(nbmail,nip_max)
REAL*8  :: Frou(nbmail,nip_max)

REAL*8  :: rho(nbmail)
REAL*8  :: ene(nbmail)
REAL*8  :: Pres(nbmail)
REAL*8  :: Temp(nbmail)
REAL*8  :: Frac(nbmail,nip_max)
REAL*8  :: dtime(nbmail)
REAL*8  :: dpde(nbmail)
REAL*8  :: cs2(nbmail)
REAL*8  :: conv(nbmail)

real r0,e0,de,efin,dt 

Temp(:)=0.
Pres(:)=0.
Frac(ii,:)=Frin(ii,:)


rho(ii)=r0
ene(ii)=e0
dtime(:)=dt
do while(ene(ii).lt.efin)
CALL S_CALC_CINE_VE(nbmail,dtime(1:nbmail),rho(1:nbmail),ene(1:nbmail),&
				  Pres(1:nbmail),Temp(1:nbmail),&
				 Frac(1:nbmail,1),Frac(1:nbmail,2),Frac(1:nbmail,3),&
				 Frac(1:nbmail,4),Frac(1:nbmail,5),Frac(1:nbmail,6),&
				 dpde(1:nbmail),cs2(1:nbmail),conv(1:nbmail))
ene(ii)=ene(ii)+de
!    write(601,'(20d)') Pres(ii),Frac(ii,1)
!    write(602,'(20d)') Pres(ii),Frac(ii,2)
!    write(603,'(20d)') Pres(ii),Frac(ii,3)
!    write(701,'(20d)') Frac(ii,1),Temp(ii)
!    write(702,'(20d)') Frac(ii,2),Temp(ii)
!    write(703,'(20d)') Frac(ii,3),Temp(ii)
!    write(801,'(20d)') Pres(ii),Temp(ii)
!    write(802,'(20d)') rho(ii),Pres(ii)
!    write(803,'(20d)') rho(ii),Temp(ii)
enddo
    write(601,'(a)') ''
    write(602,'(a)') ''
    write(603,'(a)') ''
    write(701,'(a)') ''
    write(702,'(a)') ''
    write(703,'(a)') ''
    write(801,'(a)') ''
    write(802,'(a)') ''
    write(803,'(a)') ''

Frou(ii,:)=Frac(ii,:)

end subroutine S_isochore



END MODULE CHEMIN 




PROGRAM TEST_CINE22

USE  SCALC
USE  SINIT
USE M_lec
USE buffer

!USE M_COUPLAGE_CINE, ONLY : S_CALC_CINE_VE

USE CHEMIN

IMPLICIT NONE

character*2000, fich

INTEGER,PARAMETER      :: nbmail_max=10
INTEGER,PARAMETER      :: nip_max=6


REAL*8  :: rho(nbmail_max)
REAL*8  :: ene(nbmail_max)
REAL*8  :: Pres(nbmail_max)
REAL*8  :: Temp(nbmail_max)
REAL*8  :: Frac(nbmail_max,nip_max)
REAL*8  :: dtime(nbmail_max)

REAL*8  :: Frin(1,nip_max)
REAL*8  :: Frou(1,nip_max)

REAL*8  :: dpde(nbmail_max)
REAL*8  :: cs2(nbmail_max)
REAL*8  :: conv(nbmail_max)


REAL*8  :: rho_std,ene_std


INTEGER        :: nbmail


INTEGER        :: ii,jj,kk
real dt
real ef,eh,rf,ph,th
real rho_init(4)
real drho

!=======================================================================
! definition du jeu de coeff de la cinetique et de l'EE
!=======================================================================
fich='ee.CineTest22#.Sn.00#.coeff'
!lecture coeff
CALL S_lec_coeff(fich)

CALL CineTest22_couplage_Sinit(rho_std,ene_std)

! /!\ ATTENTION PAS DE VERIFICATION DE LA COHERENCE SUR NIP ENTRE EE ET CINETIQUE 
print*, 'Etat Std'
print*, '--------'
print*, rho_std,ene_std
print*, '--------'

rho_init(1)=rho_std
rho_init(2)=7.5
rho_init(3)=8.
rho_init(4)=9.

! TEST POINT A POINT
nbmail=1
! TEST POINT A POINT

dt=1.e-10
!dt=1.e-10
!dt=1.e-6
!dt=1.e-1

Frin(1,:)=0.
Frin(1,1)=1.
Frin(1,2)=0.
Frin(1,3)=0.



!CALL S_hugo(rho_std,ene_std,9.5,dt*1.e3,80,Frin,Frou,eh,ph,th)
stop
!Frin(1,:)=0.
!Frin(1,1)=1.
!Frin(1,2)=0.
!Frin(1,3)=0.
!rf=8.7
!CALL S_hugo(rho_std,ene_std,rf,dt*1.e3,80,Frin,Frou,eh,ph,th)
!Frin(1,:)=Frou(1,:)
!CALL S_isentrope(rf,eh,0.01,7.0,dt*1.e3,Frin,Frou,ph,th,-1)
!stop

!Frin(1,:)=0.
!Frin(1,1)=1.
!Frin(1,2)=0.
!Frin(1,3)=0.
!ef=0.03
!rf=8.7
!CALL S_hugo(rho_std,ef,rf,dt*1.e3,80,Frin,Frou,eh,ph,th)


!Frin(1,:)=Frou(1,:)
!CALL S_isentrope(rf,eh,0.01,7.0,dt*1.e3,Frin,Frou,ph,th,-1)

Frin(1,:)=0.
Frin(1,1)=1.
Frin(1,2)=0.
Frin(1,3)=0.
ef=0.08
rf=8.7
!CALL S_hugo(rho_std,ef,rf,dt*1.e3,80,Frin,Frou,eh,ph,th)
Frin(1,:)=Frou(1,:)
!CALL S_isentrope(rf,eh,0.01,7.0,dt*1.e3,Frin,Frou,ph,th,-1)
stop
Frin(1,:)=0.
Frin(1,1)=1.
Frin(1,2)=0.
Frin(1,3)=0.
ef=0.15
rf=8.7
!CALL S_hugo(rho_std,ef,rf,dt*1.e3,80,Frin,Frou,eh,ph,th)
Frin(1,:)=Frou(1,:)
!CALL S_isentrope(rf,eh,0.01,7.0,dt*1.e3,Frin,Frou,ph,th,-1)


stop
! hugoniot 
Frin(1,:)=0.
Frin(1,1)=1.
Frin(1,2)=0.
Frin(1,3)=0.
ef=0.03
rf=8.5
!CALL S_hugo(rho_std,ene_std,rf,dt*1.e3,80,Frin,Frou,eh,ph,th)
Frin=Frou
!CALL S_isentrope(rf,eh,0.01,7.1,dt*1.e3,Frin,Frou,ph,th,-1)


!CALL S_hugo(rho_std,ef,8.5,dt*1.e3,80,Frin,Frou)

stop

!CALL S_isochore(7.5,ene_std,0.01,1.5,dt*1.e3)
!CALL S_isochore(7.8,ene_std,0.01,1.5,dt*1.e3)
!CALL S_isochore(8.0,ene_std,0.01,1.5,dt*1.e3)

Frin(1,:)=0.
Frin(1,1)=1.
Frin(1,2)=0.
Frin(1,3)=0.
!CALL S_isentrope(rho_std,ene_std,0.01,9.5,dt*1.e3,Frin,Frou,0.,0.,1)

Frin(1,:)=0.
Frin(1,1)=1.
Frin(1,2)=0.
Frin(1,3)=0.
ef=0.08
!CALL S_isochore(rho_std,ene_std,0.001,ef,dt*1.e3,Frin,Frou)
Frin=Frou
!CALL S_isentrope(rho_std,ef,0.01,9.5,dt*1.e3,Frin,Frou,0.,0.,1)

Frin(1,:)=0.
Frin(1,1)=1.
Frin(1,2)=0.
Frin(1,3)=0.
ef=0.15
!CALL S_isochore(rho_std,ene_std,0.001,ef,dt*1.e3,Frin,Frou)
Frin=Frou
!CALL S_isentrope(rho_std,ef,0.01,9.5,dt*1.e3,Frin,Frou,0.,0.,1)


print*, 'OK'


END 


! passage et convergence : 
! --> CST_NEWTON%ITMAX = 100 
! --> DANS S_CALC_CINE_VE_ONE frac_max=0.01 + petit mieux c'est 
! --> les evolutions Pn-1 et Tn-1 sont primordiaux permet d'avoir un bon point de départ et  frac_max plus grand (plus stable)


