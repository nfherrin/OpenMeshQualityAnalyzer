!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to analyze mesh quality indicators
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE mesh_analyze
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calcvols
CONTAINS

  !calculate volumes of each tet
  SUBROUTINE calcvols()
    REAL(8) :: totalvol1
    REAL(8) :: a(3),b(3),c(3),d(3)
    REAL(8), ALLOCATABLE :: tetvol(:),regvol(:)
    INTEGER, ALLOCATABLE :: tets_in_reg(:)
    INTEGER :: i,minreg,maxreg

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'

    minreg=MINVAL(el_tag(:))
    maxreg=MAXVAL(el_tag(:))
    ALLOCATE(tetvol(num_tets),regvol(minreg:maxreg),tets_in_reg(minreg:maxreg))
    tets_in_reg=0
    tetvol=0
    totalvol1=0
    regvol=0
    prog=0
    !compute tet volumes and add to both total volumes and region volumes
    DO i=1,num_tets
      a(:)=vertex(element(i,1),:)
      b(:)=vertex(element(i,2),:)
      c(:)=vertex(element(i,3),:)
      d(:)=vertex(element(i,4),:)
      tetvol(i)=ABS((-c(2)*d(1)+b(2)*(-c(1)+d(1))+b(1)*(c(2)-d(2))+c(1)*d(2))*(a(3)-d(3))+(a(1)-d(1)) &
        *(-c(3)*d(2)+b(3)*(-c(2)+d(2))+b(2)*(c(3)-d(3))+c(2)*d(3))+(a(2)-d(2))*(b(3)*(c(1)-d(1)) &
        +c(3)*d(1)-c(1)*d(3)+b(1)*(-c(3)+d(3))))/6
      regvol(el_tag(i))=regvol(el_tag(i))+tetvol(i)
      tets_in_reg(el_tag(i))=tets_in_reg(el_tag(i))+1
      totalvol1=totalvol1+tetvol(i)
      IF(MOD(i,CEILING(num_tets*1.0/(max_prog-1.0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)

    DO i=minreg,maxreg
      WRITE(*,'(A,I0,A,I0)')'Region ',i,' tets: ',tets_in_reg(i)
      WRITE(*,'(A,I0,A,ES24.16)')'Region ',i,' volume: ',regvol(i)
      WRITE(*,'(A,I0,A,ES24.16)')'Region ',i,' equivalent radius: ',(3.0/4.0/pi*regvol(i))**(1.0/3.0)
    ENDDO
    WRITE(*,'(A,I0)')'Total number of tets: ',SUM(tets_in_reg)
    WRITE(*,'(A,ES24.16)')'Total system volume: ',totalvol1
    WRITE(*,'(A,ES24.16)')'Equivalent radius: ',(3.0/4.0/pi*totalvol1)**(1.0/3.0)
  ENDSUBROUTINE calcvols
END MODULE mesh_analyze
