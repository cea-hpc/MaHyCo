MODULE M_polygone
CONTAINS


LOGICAL FUNCTION F_in_polygone(n_xy,Data_xy,x,y)
! renvoi TRUE si dans le polygone
IMPLICIT NONE
REAL*8,    PARAMETER   :: NPAR_prec_machine_poly=abs(1.-tan(atan(1.)))*100.         !  precision machine + large pour eos

INTEGER n_xy

REAL*8     :: Data_xy(n_xy,2),x,y
INTEGER  :: nb_crois,ishift

REAL*8 Poly_x(n_xy+1),Poly_y(n_xy+1)
REAL*8 Segment_x(2),Segment_y(2),Dx(2)
REAL*8 a,b,xmin,xmax
INTEGER ni_segment
REAL*8   x_i(2*n_xy,2)

INTEGER Nb_Poly,icount,jcount

Nb_Poly=n_xy+1
if (Nb_Poly.lt.3) then
! le polygone ne peut pas être fermé 
  F_in_polygone =.TRUE.
  RETURN
endif

! affectation et fermeture du polygone
Poly_x(1:Nb_Poly-1)=Data_xy(1:Nb_Poly-1,1)
Poly_y(1:Nb_Poly-1)=Data_xy(1:Nb_Poly-1,2)
Poly_x(Nb_Poly)=Poly_x(1)
Poly_y(Nb_Poly)=Poly_y(1)
DO icount=1,Nb_Poly
  IF (abs(Poly_x(icount)-x).le.NPAR_prec_machine_poly) THEN 
    IF (abs(Poly_y(icount)-y).le.NPAR_prec_machine_poly) THEN 
      F_in_polygone =.TRUE.
      RETURN
    ENDIF 
  ENDIF
ENDDO

ni_segment=0
F_in_polygone=.FALSE.
DO icount=1,Nb_Poly-1
    Segment_x(1)=Poly_x(icount);Segment_x(2)=Poly_x(icount+1);
    Segment_y(1)=Poly_y(icount);Segment_y(2)=Poly_y(icount+1);
    IF (y.ge.MINVAL(Segment_y,DIM=1).and.(y.le.MAXVAL(Segment_y,DIM=1))) ni_segment=ni_segment+1

    IF (abs(Segment_y(1)-Segment_y(2)).le.NPAR_prec_machine_poly) THEN 
    ni_segment=ni_segment+1
    
      IF (abs(Segment_y(1)-y).le.NPAR_prec_machine_poly) THEN
	IF (x.gt.MINVAL(Segment_x,DIM=1).and.(x.lt.MAXVAL(Segment_x,DIM=1))) THEN
	  F_in_polygone =.TRUE.
	  RETURN
	ENDIF
      ENDIF
    ENDIF
ENDDO

DO icount=1,Nb_Poly
IF (Poly_y(icount).eq.y) Poly_y(icount)=Poly_y(icount)+NPAR_prec_machine_poly
ENDDO

ni_segment=0
DO icount=1,Nb_Poly-1
Segment_x(1)=Poly_x(icount);Segment_x(2)=Poly_x(icount+1);
Segment_y(1)=Poly_y(icount);Segment_y(2)=Poly_y(icount+1);
IF (y.ge.MINVAL(Segment_y,DIM=1).and.(y.le.MAXVAL(Segment_y,DIM=1))) THEN 
  IF (abs(Segment_y(1)-Segment_y(2)).le.NPAR_prec_machine_poly) THEN
  ni_segment=ni_segment+1;
   x_i(ni_segment,1)=Segment_x(1);x_i(ni_segment,2)=1.
  ni_segment=ni_segment+1;
   x_i(ni_segment,1)=Segment_x(2);x_i(ni_segment,2)=1.
  ELSE
  ni_segment=ni_segment+1
  a=((Segment_x(1)-Segment_x(2))/(Segment_y(1)-Segment_y(2)))
  b=Segment_x(1)-a*Segment_y(1)
  x_i(ni_segment,1)=a*y+b;;x_i(ni_segment,2)=0.;
  ENDIF

  IF (abs(Segment_x(1)-Segment_x(2)).le.NPAR_prec_machine_poly) THEN
! cas particulier de ligne verticale
   if (abs(x-Segment_x(2)).le.NPAR_prec_machine_poly) THEN
   F_in_polygone =.TRUE.
   RETURN
   ENDIF
  ENDIF
ENDIF
ENDDO
!limite min et max
xmin=huge(xmin)
xmax=-huge(xmax)
DO icount=1,ni_segment
if (x_i(icount,1).le.xmin) xmin=x_i(icount,1)
if (x_i(icount,1).ge.xmax) xmax=x_i(icount,1)
ENDDO

if (x.lt.xmin) then
F_in_polygone =.FALSE.
RETURN
endif
if (x.gt.xmax) then
F_in_polygone =.FALSE.
RETURN
endif



! cas particulier de ligne horizontale
DO icount=1,ni_segment-1
  IF (x_i(icount,2).gt.0.) THEN 
  Dx(1)=x_i(icount,1);Dx(2)=x_i(icount+1,1);
  IF (x.gt.MINVAL(Dx,DIM=1).and.(x.lt.MAXVAL(Dx,DIM=1))) THEN
  F_in_polygone=.TRUE.
RETURN
  ENDIF
ENDIF
ENDDO


! supression des doubles points
DO icount=1,ni_segment-1
 DO jcount=icount+1,ni_segment
    IF (abs(x_i(icount,1)-x_i(jcount,1)).lt.1.e-10) x_i(jcount,2)=2.
ENDDO
ENDDO
nb_crois=0;ishift=0
DO icount=1,ni_segment
IF (x_i(icount,2).lt.0.1) then
ishift=ishift+1
IF (x.ge.x_i(icount,1)) nb_crois=nb_crois+1
ENDIF
ENDDO


ishift=mod(ishift,2)
nb_crois=nb_crois+ishift
if (mod(nb_crois,2).eq.1) F_in_polygone=.TRUE.

RETURN
END FUNCTION F_in_polygone
END MODULE M_polygone
