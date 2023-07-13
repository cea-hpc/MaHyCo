!! modele interpolation DE KERLEY UTILISE PAR EOS ET LOCATE DE NUMREC
MODULE M_interpol_EOS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE indexx(n,arr,indx)
      IMPLICIT NONE
      INTEGER n,indx(n),M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL*8 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
!!        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END SUBROUTINE indexx

      SUBROUTINE sort_m(n,m,im,ra) ! tableau n,m et dont la colonne im est trie
      IMPLICIT NONE
      INTEGER n,m,im,iwksp(n)
      REAL*8 ra(n,m),wksp(n)
      INTEGER j,k
      call indexx(n,ra(1:n,im),iwksp)
      do  k=1,m
        do  j=1,n;wksp(j)=ra(j,k);enddo;
        do  j=1,n;ra(j,k)=wksp(iwksp(j));enddo
      enddo
      return
      END SUBROUTINE  sort_m
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE sort_n(n,m,in,ra) ! tableau n,m et dont la ligne in est trie
      IMPLICIT NONE
      INTEGER n,m,in,iwksp(m)
      REAL*8 ra(n,m),wksp(m)
      INTEGER j,k
      call indexx(m,ra(in,1:m),iwksp)

      do  k=1,n
        do  j=1,m;wksp(j)=ra(k,j);enddo;
        do  j=1,m;ra(k,j)=wksp(iwksp(j));enddo
      enddo


      return
      END SUBROUTINE  sort_n
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE sort1(n,ra)
      IMPLICIT NONE
      INTEGER n,iwksp(n)
      REAL*8 ra(n),wksp(n)
      INTEGER j
      call indexx(n,ra,iwksp)

      do  j=1,n;wksp(j)=ra(j);enddo; !sauv 
      do  j=1,n;ra(j)=wksp(iwksp(j));enddo; ! recopie
      
      return
      END SUBROUTINE  sort1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SUBROUTINE KERLEY EN 1 DIMENSION 
SUBROUTINE S_INTERPOLE_RATF1_eos (XT,FT,N,X,FRF1,DF)												 
!-----------------------------------------------------------------------									 
!																		 
!  PURPOSE:        INTERPOLATE FOR A FUNCTION F(X) FROM TABLES OF										 
!                  F AND X.  USES RATIONAL FUNCTION METHOD WITH 										 
!                  QUADRATIC ESTIMATE OF DERIVATIVES AT END POINTS.										 
!																		 
!  ARGUMENTS:      X  (INP) - INDEPENDENT VARIABLE												 
!                  XT (INP) - TABLE OF INDEPENDENT VARIABLE											 
!                  FT (INP) - TABLE OF DEPENDENT VARIABLE											 
!                  N  (INP) - LENGTH OF ARRAYS XT AND FT											 
!                  DF (OUT) - DERIVATIVE OF FUNCTION												 
!                  FRF1 IS THE VALUE OF THE FUNCTION AT X											 
!																		 
!  PROGRAMMER:     G. I. KERLEY, T-4.														 
IMPLICIT NONE																	 
      INTEGER N 																 
      REAL*8 XT(N),FT(N),X,FRF1,DF														 
      INTEGER I,J,JP																 
      REAL*8 Q,D,R,S,SP,C1,C2,C3,C4,DM,SM 													 
      																		 
!  SEARCH FOR INDEX																 
																		 
      I = 1																	 
      J = N																	 
      !verification que xt(1)<xt(n)														 
      IF (XT(I).GT.XT(J)) THEN															 
      print*, XT(I),XT(J)															 
      WRITE(*,*) "RATF1 doit etre croissant en X";STOP; 											 
      ENDIF																	 
      																		 
 1    IF(J-I.EQ.1) GO TO 3															 
      JP = .5*(J+I)																 
      IF(X.LT.XT(JP))GO TO 2															 
      I = JP																	 
      GO TO 1																	 
 2    J = JP																	 
      GO TO 1																	 
!  COMPUTE INTERPOLATION FUNCTION														 
																		 
 3    Q = X-XT(I)																 
      D = XT(I+1)-XT(I) 															 
      IF (D.EQ.0.) THEN
        FRF1=HUGE(FRF1)  
        DF=HUGE(DF)
	RETURN
      ENDIF																 
      R = D-Q																	 
      S = (FT(I+1)-FT(I))/D															 
    																		 
      IF(I.GT.1) GO TO 4
      
      IF ((XT(I+2)-XT(I+1)).EQ.0.) THEN
        FRF1=HUGE(FRF1)  
        DF=HUGE(DF)
	RETURN
      ENDIF																 
      SP = (FT(I+2)-FT(I+1))/(XT(I+2)-XT(I+1))													 

      IF ((XT(I+2)-XT(I)).EQ.0.) THEN
        FRF1=HUGE(FRF1)  
        DF=HUGE(DF)
	RETURN
      ENDIF																 

      C2 = (SP-S)/(XT(I+2)-XT(I))														 
      IF(S*(S-D*C2).LE.0.) C2=S/D														 
      FRF1 = FT(I)+Q*(S-R*C2)															 
      DF = S+(Q-R)*C2																 
      RETURN																	 
 4    DM = XT(I)-XT(I-1)
      IF (DM.EQ.0.) THEN
        FRF1=HUGE(FRF1)  
        DF=HUGE(DF)
	RETURN
      ENDIF																 
      SM = (FT(I)-FT(I-1))/DM															 
      C1 = (S-SM)/(D+DM)															 
      IF(I.LT.N-1) GO TO 5															 
      FRF1 = FT(I)+Q*(S-R*C1)															 
      DF = S+(Q-R)*C1																 
      RETURN																	 
 5    IF(I.GT.2) GO TO 6															 
      IF(SM*(SM-DM*C1).LE.0.) C1=(S-SM-SM)/D													 
 6    SP = 1. 
      IF ((XT(I+2)-XT(I+1)).EQ.0.) THEN
        FRF1=HUGE(FRF1)  
        DF=HUGE(DF)
	RETURN
      ENDIF																 
      SP = (FT(I+2)-FT(I+1))/(XT(I+2)-XT(I+1))													 
      C2 = (SP-S)/(XT(I+2)-XT(I))														 
      C3 = ABS(C2*R)																 
      C4 = C3+ABS(C1*Q) 
      IF(C4.GT.0.) C3=C3/C4															 
      C4 = C2+C3*(C1-C2)															 
      FRF1 = FT(I)+Q*(S-R*C4)															 
      DF = S+(Q-R)*C4+D*(C4-C2)*(1.-C3) 													 
      
      RETURN																	 
END SUBROUTINE S_INTERPOLE_RATF1_eos														 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
																		 
																		 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! location dans une grille xx															
      SUBROUTINE locate_eos(xx,n,x,j)	!ok													
      INTEGER j,n																
      REAL*8 x,xx(n)																
      INTEGER jl,jm,ju																
      jl=0																	
      ju=n+1																	
      do while (ju-jl.gt.1)															
        jm=(ju+jl)/2																
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then												
          jl=jm 																
        else																	
          ju=jm 																
        endif																	
      enddo																	
																		
      if(x.eq.xx(1))then															
        j=1																	
      else if(x.eq.xx(n))then															
        j=n-1																	
      else																	
        j=jl																	
        !!if (xx(1).ge.xx(n)) j=j+1 
      endif																	
      ! borne min atteinte															
      if (j.eq.0) j=1																
      																		
      return																	
      END SUBROUTINE locate_eos														
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine EXTRA2(nxt, nyt,xtbls, ytbls, f1tbls,xvalv, &
      yvalv, f1v, f1dxv, f1dyv,ixv, iyv)

!
!-----Interpolation d une fonction de 2 variables independantes
!     par une fonction lineaire log-log
!     (16/12/96 A.B.)
!
      IMPLICIT NONE
      INTEGER nxt,nyt

      REAL*8    xtbls(nxt), ytbls(nyt)
      REAL*8    f1tbls(nxt,nyt)
      REAL*8    xvalv, yvalv
      REAL*8    f1v, f1dxv, f1dyv
      INTEGER ixv, iyv,i,j
      REAL*8    xi,xip1,dx1,dx2
      REAL*8    dqx,dqy,dqzj,dqzjp1
      REAL*8    zij,zip1j,zijp1,zip1jp1
      REAL*8    f1vp1,deltat
      REAL*8    cstmin,valmin
      REAL*8, PARAMETER::     par_cstmin=1.

      dqy = 0.
!
!.....calcul indices et fonctions aux 4 coins du domaine
!
      i = ixv
      j = iyv

      if ((xvalv.lt.xtbls(i)).and.(i.gt.1)) i=i-1 
      if ((xvalv.gt.xtbls(i+1)).and.(i+1.lt.nxt)) i=i+1 
      if ((yvalv.lt.ytbls(j)).and.(j.gt.1)) j=j-1 
      if ((yvalv.gt.ytbls(j+1)).and.(j+1.lt.nyt)) j=j+1 

      xi   = xtbls(i)
      xip1 = xtbls(i+1)
      dx1  = xvalv  -xi
      dx2  = xip1-xvalv

! definition de cstmin pour le log 
      valmin = MINVAL(f1tbls(i:i+1,j:j+1))
      if (valmin.gt.par_cstmin) then
        cstmin=0.
      else
        cstmin=-valmin+par_cstmin
      endif


!
      zij     = f1tbls(i,j)    +cstmin
      zip1j   = f1tbls(i+1,j)  +cstmin
      zijp1   = f1tbls(i,j+1)  +cstmin
      zip1jp1 = f1tbls(i+1,j+1)+cstmin

!.....calcul de la Fonction et de ses 2 Derivees partielles
!
!.....calcul des coefficients de ponderation
!
      dqx  = log(xip1 / xi)
      dqzj  = log(zip1j / zij)
      dqzjp1  = log(zip1jp1 / zijp1)
!
!.....calcul de Z(x,y)
!
      f1v = zip1j  * exp(log(xvalv/xip1)*dqzj/dqx)
      f1vp1 = zip1jp1  * exp(log(xvalv/xip1)*dqzjp1/dqx)
!
!.....calcul de dZ/dx(x,y)
!
!     f1dxv = f1v * dqy / (dqx*xvalv)
      deltaT = (ytbls(j+1) - yvalv) / (ytbls(j+1)-ytbls(j))
      
      
      
      f1dxv =  deltaT * f1v *  dqzj / (dqx*xvalv)+ (1-deltaT) * &
      f1vp1 *  dqzjp1 / (dqx*xvalv)
!
!.....calcul de dZ/dy(x,y)
!
      f1dyv = (f1vp1-f1v)/(ytbls(j+1)-ytbls(j))
!     *    + dqy*(zz2  -zz1
!     *    + aqx*(zij  -zijp1  )
!     *    +  qx*(zip1j-zip1jp1))

      f1v = f1v + f1dyv * (yvalv - ytbls(j))
      
!     maj en fonction de   cstmin
      f1v = f1v -cstmin
      
      

      return
      end subroutine EXTRA2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine S_INTERPOLL2(nxt, nyt, nxtbl, nytbl,xtbls, ytbls, &
      f1tbls,xvalv, yvalv, f1v, f1dxv,f1dyv,ixv,iyv)

!
!-----Interpolation d une fonction de 2 variables independantes
!     par la methode de KERLEY (fractions rationnelles)
!     (06/05/92 G.D.)
!

      INTEGER nxt, nyt,nxtbl,nytbl
      REAL*8    xtbls(nxt), ytbls(nyt)
      REAL*8    f1tbls(nxt,nyt)
      REAL*8    xvalv, yvalv
      REAL*8    f1v, f1dxv, f1dyv
      INTEGER ixv, iyv,i,j
      REAL*8    xi,xip1,dx1,dx2
      REAL*8    d0,d1,d2,dqx,d2pd1,d1pd0
      REAL*8    s0,s1,s2,c1,c2
      REAL*8    u1,u2,tn,q,qx,qy
      REAL*8    aqx,aqy
      REAL*8    yj,yjp1,dy1,dy2,dqy
      REAL*8    zij,zip1j,zijp1,zip1jp1
      REAL*8    zdd,zz1,zd1,zz2,zd2,zz3,zd3,zz4,zd4
!
!.....calcul indices et fonctions aux 4 coins du domaine
!
!.....calcul de Z(x,yj+1) et Z(x,yj)
!
      i = ixv
      j = iyv

      if ((xvalv.lt.xtbls(i)).and.(i.gt.1)) i=i-1 
      if ((xvalv.gt.xtbls(i+1)).and.(i+1.lt.nxt)) i=i+1 
      if ((yvalv.lt.ytbls(j)).and.(j.gt.1)) j=j-1 
      if ((yvalv.gt.ytbls(j+1)).and.(j+1.lt.nyt)) j=j+1 

      xi   = xtbls(i)
      xip1 = xtbls(i+1)
      dx1  = xvalv  -xi
      dx2  = xip1-xvalv
!
      zij = f1tbls(i,j)
      zip1j = f1tbls(i+1,j)
      zijp1 = f1tbls(i,j+1)
      zip1jp1 = f1tbls(i+1,j+1)

      if (i.eq.1) then
         d1 = xip1    -xi
         d2 = xtbls(i+2)-xip1
         dqx= d1
         d2pd1= d2+d1
         s1 = (   zip1j   -zij  )/d1
         s2 = (f1tbls(i+2,j)-zip1j)/d2
         c2 = (s2-s1)/(d2pd1)
         if (s1*(s1-d1*c2).le.0.) c2 = s1/d1
         zdd =          (s1-c2*dx2)
         zz1 = zij + dx1*zdd
         zd1 =           zdd + dx1*c2
         s1 = (   zip1jp1   -zijp1  )/d1
         s2 = (f1tbls(i+2,j+1)-zip1jp1)/d2
         c2 = (s2-s1)/(d2pd1)
         if (s1*(s1-d1*c2).le.0.) c2 = s1/d1
         zdd =          (s1-c2*dx2)
         zz2 = zijp1 + dx1*zdd
         zd2 =             zdd + dx1*c2
         goto 10
      endif
      if (i.eq.nxtbl) then
         d1 = xi  -xtbls(i-1)
         d2 = xip1-xi  
         dqx= d2
         d2pd1= d2+d1
         s1 = (zij  -f1tbls(i-1,j))/d1
         s2 = (zip1j-   zij   )/d2
         c1 = (s2-s1)/(d2pd1)
         if (s2*(s2-c1*(dx2-dx1)).lt.0.) c1 = s2/d2
         zdd =          (s2-c1*dx2)
         zz1 = zij + dx1*zdd
         zd1 =           zdd + dx1*c1
         s1 = (zijp1  -f1tbls(i-1,j+1))/d1
         s2 = (zip1jp1-   zijp1   )/d2
         c1 = (s2-s1)/(d2pd1)
         if (s2*(s2-c1*(dx2-dx1)).lt.0.) c1 = s2/d2
         zdd =          (s2-c1*dx2)
         zz2 = zijp1 + dx1*zdd
         zd2 =             zdd + dx1*c1
         goto 10
      endif
      d0 = xi      -xtbls(i-1)
      d1 = xip1    -xi  
      d2 = xtbls(i+2)-xip1
      dqx= d1
      d1pd0= d1+d0
      d2pd1= d2+d1
      s0 = (   zij     -f1tbls(i-1,j))/d0
      s1 = (   zip1j   -   zij   )/d1
      s2 = (f1tbls(i+2,j)-   zip1j )/d2
      c1 = (s1-s0)/(d1pd0)
      c2 = (s2-s1)/(d2pd1)
      u1 = abs(c2*dx2)
      u2 = abs(c1*dx1)
      tn = c1*u1 + c2*u2
      if (tn.eq.0.) then
         zz1 = zij + dx1*s1
         zd1 =           s1
      else
         q   = tn/(u1+u2)
         zdd =          (s1-q*dx2)
         zz1 = zij + dx1*zdd
         zd1 =           zdd + dx1*q
      endif
      s0 = (   zijp1     -f1tbls(i-1,j+1))/d0
      s1 = (   zip1jp1   -   zijp1   )/d1
      s2 = (f1tbls(i+2,j+1)-   zip1jp1 )/d2
      c1 = (s1-s0)/(d1pd0)
      c2 = (s2-s1)/(d2pd1)
      u1 = abs(c2*dx2)
      u2 = abs(c1*dx1)
      tn = c1*u1 + c2*u2
      if (tn.eq.0.) then
         zz2 = zijp1 + dx1*s1
         zd2 =             s1
      else
         q   = tn/(u1+u2)
         zdd =          (s1-q*dx2)
         zz2 = zijp1 + dx1*zdd
         zd2 =             zdd + dx1*q
      endif
 10   continue
!
!.....calcul de Z(xi,y) et de Z(xi+1,y)
!
      yj   = ytbls(j)
      yjp1 = ytbls(j+1)
      dy1  = yvalv   -yj
      dy2  = yjp1-yvalv
      if (j.eq.1) then
         d1 = yjp1    -yj  
         d2 = ytbls(j+2)-yjp1
         dqy= d1
         d2pd1= d2+d1
         s1 = (   zijp1     -zij  )/d1
	 ! CORRECTION f1tbls(i,j+2) au lieu de f1tbls(i+2,j+1)
         s2 = (f1tbls(i,j+2)-zijp1)/d2
         c2 = (s2-s1)/(d2pd1)
         if (s1*(s1-d1*c2).le.0.) c2 = s1/d1
         zdd =          (s1-c2*dy2)
         zz3 = zij + dy1*zdd
         zd3 =           zdd + dy1*c2
         s1 = (   zip1jp1     -zip1j  )/d1
	 ! CORRECTION f1tbls(i+1,j+2) au lieu de f1tbls(i+2,j+1)
         s2 = (f1tbls(i+1,j+2)-zip1jp1)/d2
         c2 = (s2-s1)/(d2pd1)
         if (s1*(s1-d1*c2).le.0.) c2 = s1/d1
         zdd =          (s1-c2*dy2)
         zz4 = zip1j + dy1*zdd
         zd4 =             zdd + dy1*c2
         goto 30
      endif
      if (j.eq.nytbl) then
         d1 = yj - ytbls(j-1)
         d2 = yjp1 - yj  
         dqy= d2
         d2pd1= d2+d1
         s1 = (zij -f1tbls(i,j-1))/d1
         s2 = (zijp1- zij )/d2
         c1 = (s2-s1)/(d2pd1)
         if (s2*(s2-c1*(dy2-dy1)).lt.0.) c1 = s2/d2
         zdd =          (s2-c1*dy2)
         zz3 = zij + dy1*zdd
         zd3 =           zdd + dy1*c1
         s1 = (zip1j  -f1tbls(i+1,j-1))/d1
         s2 = (zip1jp1-   zip1j     )/d2
         c1 = (s2-s1)/(d2pd1)
         if (s2*(s2-c1*(dy2-dy1)).lt.0.) c1 = s2/d2
         zdd =          (s2-c1*dy2)
         zz4 = zip1j + dy1*zdd
         zd4 =             zdd + dy1*c1
         goto 30
      endif
      d0 = yj      -ytbls(j-1)
      d1 = yjp1    -yj  
      d2 = ytbls(j+2)-yjp1
      dqy= d1
      d1pd0= d1+d0
      d2pd1= d2+d1
      s0 = (   zij       -f1tbls(i,j-1))/d0
      s1 = (   zijp1     -   zij     )/d1
      s2 = (f1tbls(i,j+2)-   zijp1   )/d2
      c1 = (s1-s0)/(d1pd0)
      c2 = (s2-s1)/(d2pd1)
      u1 = abs(c2*dy2)
      u2 = abs(c1*dy1)
      tn = c1*u1 + c2*u2
      if (tn.eq.0.) then
         zz3 = zij + dy1*s1
         zd3 =           s1
      else
         q   = tn/(u1+u2)
         zdd =          (s1-q*dy2)
         zz3 = zij + dy1*zdd
         zd3 =           zdd + dy1*q
      endif
      s0 = (   zip1j       -f1tbls(i+1,j-1))/d0
      s1 = (   zip1jp1     -   zip1j     )/d1
      s2 = (f1tbls(i+1,j+2)-   zip1jp1   )/d2
      c1 = (s1-s0)/(d1pd0)
      c2 = (s2-s1)/(d2pd1)
      u1 = abs(c2*dy2)
      u2 = abs(c1*dy1)
      tn = c1*u1 + c2*u2
      if (tn.eq.0.) then
         zz4 = zip1j + dy1*s1
         zd4 =             s1
      else
         q   = tn/(u1+u2)
         zdd =          (s1-q*dy2)
         zz4 = zip1j + dy1*zdd
         zd4 =             zdd + dy1*q
      endif
 30   continue
!
!.....calcul de la Fonction et de ses 2 Derivees partielles
!
!.....calcul des coefficients de ponderation
!
      dqx  = 1./dqx
      dqy  = 1./dqy
       qx  = dx1*dqx
       qy  = dy1*dqy
      aqx  = 1.-qx
      aqy  = 1.-qy
!
!.....calcul de Z(x,y)
!
      f1v = zz1*aqy + zz2*qy + zz3*aqx + zz4*qx- zij	 *aqx*aqy- &
      zijp1  *aqx* qy- zip1j  * qx*aqy- zip1jp1* qx* qy
!
!.....calcul de dZ/dx(x,y)
!
      f1dxv = aqy*zd1 + qy*zd2+ dqx*(zz4  -zz3+ aqy*(zij  -zip1j  )+&
        qy*(zijp1-zip1jp1))
!
!.....calcul de dZ/dy(x,y)
!
      f1dyv = aqx*zd3 + qx*zd4 + dqy*(zz2  -zz1 + aqx*(zij  -zijp1  )+&
        qx*(zip1j-zip1jp1))

      return
      end subroutine S_INTERPOLL2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     externals -- srchrv,  min, max, abs.
!
!     purpose   -- performs vector bi-rational-function (12-point)
!                  search-interpolation on eos function f(x,y).
!                  uses linear extrapolation.
!     arguments
!	nxtbl  : (IN) nb de tabulations selon x
!	nytbl  : (IN) nb de tabulations selon y
!	xtbls  : (IN) tabulations selon x
!	ytbls  : (IN) tabulations selon y
!	f1tbls : (IN) tabulations de la fonction
!       nzons  : (IN) nb de mailles
!	xvalv  : (IN) valeurs de x
!	yvalv  : (IN) valeurs de y
!	f1v    : (OUT) valeurs interpol\E9es de la fonction
!	f1dxv  : (OUT) valeurs interpol\E9es de la d\E9riv\E9e en x
!	f1dyv  : (OUT) valeurs interpol\E9es de la d\E9riv\E9e en y
!       ixv    : (IN) positions des x dans la grille
!	iyv    : (IN) positions des y dans la grille
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine S_INTERPOLV2(nxtbl,nytbl,xtbls,ytbls,f1tbls,nzons,xvalv&
      ,yvalv,f1v,f1dxv,f1dyv,ixv,iyv)
      INTEGER nxtbl,nytbl
      INTEGER nzons
      REAL*8    zero,un
      REAL*8    xtbls(nxtbl), ytbls(nytbl)
      REAL*8    f1tbls(nxtbl,nytbl)
      REAL*8    xvalv(nzons), yvalv(nzons)
      REAL*8    f1v(nzons), f1dxv(nzons), f1dyv(nzons)
      INTEGER ixv(nzons), iyv(nzons)
      INTEGER i,jx,jy,jv
      REAL*8    rx3,rx4,ry3,ry4,syv,tyv
      REAL*8    f1x1,f1xx1,g1xx1,f1x2,f1xx2,g1xx2
      REAL*8    f1y1,f1yy1,g1yy1,f1y2,f1yy2,g1yy2
      REAL*8    w1x1,w1y1,w1y2,g1xy1,g1yx1,f1xy1
      REAL*8    sxv,txv,w1x2
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      parameter (zero=0.,un=1.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL*8 dxg,dyg,rxv,ryv
      dxg(i,jv) = xtbls(jx + (i-1)) - xtbls(jx + (i-2))
      dyg(i,jv) = ytbls(jy + (i-1)) - ytbls(jy + (i-2))
      rxv(jv)   = ( xvalv(jv)-xtbls(jx) )/dxg(2,jx)
      ryv(jv)   = ( yvalv(jv)-ytbls(jy) )/dyg(2,jy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!.... perform vector table searches to load search vectors containing
!.... indices and spacings of nearest data table x,y values to x,yvalv.
!
!.... interpolate values of eos function f(x,y) and derivatives.
!
      do 30 jv = 1, nzons
      
         jx = ixv(jv) !+1
         jy = iyv(jv) !+1
! loc en T max T >= Tmax-1
         if(yvalv(jv).ge.ytbls(nytbl-1)) then  
           call S_interpoll2(nxtbl, nytbl,nxtbl-1,nytbl-1,xtbls,ytbls,&
	   f1tbls,xvalv(jv),yvalv(jv),f1v(jv),f1dxv(jv),f1dyv(jv),ixv(jv),iyv(jv)+1)
! loc en x x >= xmax-1 et x <= xmax
       else if(xvalv(jv).ge.xtbls(nxtbl-1).and.xvalv(jv).lt.xtbls(nxtbl)) then              
           call S_interpoll2(nxtbl, nytbl,nxtbl-1,nytbl-1,xtbls,ytbls,f1tbls,&
	   xvalv(jv),yvalv(jv),f1v(jv),f1dxv(jv),f1dyv(jv),ixv(jv)+1,iyv(jv)+1)
! loc en x x >= xmax
         else if(xvalv(jv).ge.xtbls(nxtbl)) then
       	    call extra2(nxtbl, nytbl, xtbls, ytbls,f1tbls, xvalv(jv), yvalv(jv),&
	    f1v(jv), f1dxv(jv), f1dyv(jv),ixv(jv)+1,iyv(jv))
	else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.... attention, l appel a searchv a partir de ro(1) et t(1)
!..............induit le +1 dans les ixv et iyv
!............a besoin de jx,jy -1, 0, 1 , 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rx3  = -dxg(1,jv)/dxg(2,jv)
        rx4  =  dxg(3,jv)/dxg(2,jv) + 1.
        sxv  =  min(rx4, max(rx3, rxv(jv)))
        txv  =  min( un, max( zero, rxv(jv)))
        ry3  = -dyg(1,jv)/dyg(2,jv)
        ry4  =  dyg(3,jv)/dyg(2,jv) + 1.
        syv  =  min(ry4, max(ry3, ryv(jv)))
        tyv  =  min( un, max( zero, ryv(jv)))
! -- calcul de : f1v et f1dyv
        f1x1 =  f1tbls(jx+1,jy  )-f1tbls(jx  ,jy  )
        f1xx1= (f1x1-(f1tbls(jx-1,jy  )-f1tbls(jx  ,jy  ))/rx3)/(1.-rx3)
        g1xx1= (f1x1-(f1tbls(jx+2,jy  )-f1tbls(jx  ,jy  ))/rx4)/(1.-rx4)- f1xx1
        f1x2 =  f1tbls(jx+1,jy+1)-f1tbls(jx  ,jy+1)
        f1xx2= (f1x2-(f1tbls(jx-1,jy+1)-f1tbls(jx  ,jy+1))/rx3)/(1.-rx3)
        g1xx2= (f1x2-(f1tbls(jx+2,jy+1)-f1tbls(jx  ,jy+1))/rx4)/(1.-rx4)- f1xx2
        f1y1 =  f1tbls(jx  ,jy+1) - f1tbls(jx  ,jy  )
        f1yy1= (f1y1-(f1tbls(jx  ,jy-1)-f1tbls(jx  ,jy  ))/ry3)/(1.-ry3)
        g1yy1= (f1y1-(f1tbls(jx  ,jy+2)-f1tbls(jx  ,jy  ))/ry4)/(1.-ry4)- f1yy1
        f1y2 =  f1tbls(jx+1,jy+1) - f1tbls(jx+1,jy  )
        f1yy2= (f1y2-(f1tbls(jx+1,jy-1)-f1tbls(jx+1,jy  ))/ry3)/(1.-ry3)
        g1yy2= (f1y2-(f1tbls(jx+1,jy+2)-f1tbls(jx+1,jy  ))/ry4)/(1.-ry4)- f1yy2
        w1x1 =  txv*abs(f1xx1)
        w1x1 =  w1x1/(w1x1+(1.-txv)*abs(f1xx1+g1xx1)+1.e-37)
        w1x2 =  txv*abs(f1xx2)
        w1x2 =  w1x2/(w1x2+(1.-txv)*abs(f1xx2+g1xx2)+1.e-37)
        w1y1 =  tyv*abs(f1yy1)
        w1y1 =  w1y1/(w1y1+(1.-tyv)*abs(f1yy1+g1yy1)+1.e-37)
        w1y2 =  tyv*abs(f1yy2)
        w1y2 =  w1y2/(w1y2+(1.-tyv)*abs(f1yy2+g1yy2)+1.e-37)
        f1xx1=  f1xx1 + w1x1*g1xx1
        g1xy1=  f1xx2 + w1x2*g1xx2 - f1xx1
        f1yy1=  f1yy1 + w1y1*g1yy1
        g1yx1=  f1yy2 + w1y2*g1yy2 - f1yy1
        f1xy1=  f1x2 - f1x1 - (1.-sxv-sxv)*g1xy1 - (1.-syv-syv)*g1yx1
        f1x1 =  f1x1 - (1.-sxv-sxv)*f1xx1 - syv*syv*g1yx1
        f1y1 =  f1y1 - (1.-syv-syv)*f1yy1 - sxv*sxv*g1xy1
        f1v(jv) =  f1tbls(jx  ,jy  ) + rxv(jv)*f1x1 +ryv(jv)*(f1y1 + rxv(jv)*f1xy1)&
	 -(sxv*sxv*f1xx1 + syv*syv*f1yy1)
        f1dxv(jv) = (f1x1 + ryv(jv)*(f1xy1 - w1x2*(1.-w1x2)*g1xx2) -(1.-ryv(jv))*&
	w1x1*(1.-w1x1)*g1xx1)/dxg(2,jv)
        f1dyv(jv) = (f1y1 + rxv(jv)*(f1xy1 - w1y2*(1.-w1y2)*g1yy2) -(1.-rxv(jv))*&
	w1y1*(1.-w1y1)*g1yy1)/dyg(2,jv)
        endif
   30   continue
        
	


      end subroutine S_INTERPOLV2





END MODULE M_interpol_EOS
