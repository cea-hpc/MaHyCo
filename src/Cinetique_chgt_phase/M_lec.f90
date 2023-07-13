MODULE buffer


REAL(KIND=8), ALLOCATABLE  ::  fcoeff(:)
integer ibuffer

END MODULE buffer


MODULE M_lec
CONTAINS 

SUBROUTINE S_lec_coeff(fich)
USE buffer
IMPLICIT NONE
INTEGER ios(5)
INTEGER ii,jj
integer itab_on
character*200, tmpc
character*200, tmpc2
character*60, fich
real xval(100),xx 
write(*,'(2a)') 'ouverture du fichier ',trim(fich)

open(unit=40,file=trim(fich))
ios(:)=0
itab_on=0
ibuffer=0
do while (ios(1).eq.0) 
   read(40,'(a)',IOSTAT=ios(1)) tmpc
   read(tmpc,*,IOSTAT=ios(2)) tmpc2,tmpc2,xx

    if ((ios(2).eq.0).and.(itab_on.ne.1)) then 
      if (tmpc(ii:ii).ne.'#') then 
      ibuffer=ibuffer+1 
      endif
    else 
      do ii=1,LEN_trim(tmpc)
        if (tmpc(ii:ii).eq.'#') then 
        !commentaire 
          exit
        elseif (tmpc(ii:ii+2).eq.'= [') then 
         itab_on=1
         exit
        elseif (tmpc(ii:ii).eq.']') then 
         itab_on=0
         exit
        endif 
      enddo
    endif

   if ((itab_on.eq.1)) then 
      xval(:)=HUGE(xx )
      read(tmpc,*,IOSTAT=ios(3)) xval
      
      if ((xval(1).eq.0).and.(xval(2).gt.0.9*HUGE(xx))) then
      ! on considere que l'on a un i,n, etc et la leture conduit à un 0
      xval(1)=xval(2)
      endif
      
      jj=1
      do while (xval(jj).lt.0.9*HUGE(xx))
        ibuffer=ibuffer+1 
        jj=jj+1
      enddo
   endif
enddo
! nb de mots ibuffer
close(40)

IF (ALLOCATED(fcoeff))	 DEALLOCATE(fcoeff);
ALLOCATE(fcoeff(ibuffer))

open(unit=40,file=trim(fich))
ios(:)=0
itab_on=0
ibuffer=0
do while (ios(1).eq.0) 
  read(40,'(a)',IOSTAT=ios(1)) tmpc
    read(tmpc,*,IOSTAT=ios(2)) tmpc2,tmpc2,xx

    if ((ios(2).eq.0).and.(itab_on.ne.1)) then 
      if (tmpc(ii:ii).ne.'#') then 
      ibuffer=ibuffer+1 
      fcoeff(ibuffer)=xx
      endif
    else 
      do ii=1,LEN_trim(tmpc)
        if (tmpc(ii:ii).eq.'#') then 
        !commentaire 
          exit
        elseif (tmpc(ii:ii+2).eq.'= [') then 
         itab_on=1
         exit
        elseif (tmpc(ii:ii).eq.']') then 
         itab_on=0
         exit
        endif 
      enddo
    endif


   if ((itab_on.eq.1)) then 
      xval(:)=HUGE(xx )
      read(tmpc,*,IOSTAT=ios(3)) xval

      if ((xval(1).eq.0).and.(xval(2).gt.0.9*HUGE(xx))) then
      ! on considere que l'on a un i,n, etc et la leture conduit à un 0
      xval(1)=xval(2)
      endif

      jj=1
      do while (xval(jj).lt.0.9*HUGE(xx))
        ibuffer=ibuffer+1 
        fcoeff(ibuffer)=xval(jj)
        jj=jj+1
      enddo
   endif
enddo

write(*,'(2a)') ' fin routine ok'
! nb de mots ibuffer

close(40)
END SUBROUTINE S_lec_coeff

END MODULE M_lec


