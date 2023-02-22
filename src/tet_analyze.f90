!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to analyze tet quality indicators
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE tet_analyze
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_tet_vols,comp_tet_skew,comp_tet_ar
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate volumes of each tet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calc_tet_vols()
    TYPE(vertex_type) :: a,b,c,d
    INTEGER :: i

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=0

    !compute tet volumes
    DO i=1,tot_tets
      a=tet(i)%corner(1)%p
      b=tet(i)%corner(2)%p
      c=tet(i)%corner(3)%p
      d=tet(i)%corner(4)%p
      tet(i)%vol=ABS((-c%y*d%x+b%y*(-c%x+d%x)+b%x*(c%y-d%y)+c%x*d%y)*(a%z-d%z)+(a%x-d%x) &
          *(-c%z*d%y+b%z*(-c%y+d%y)+b%y*(c%z-d%z)+c%y*d%z)+(a%y-d%y)*(b%z*(c%x-d%x) &
          +c%z*d%x-c%x*d%z+b%x*(-c%z+d%z)))/6

      IF(MOD(i,CEILING(tot_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
  ENDSUBROUTINE calc_tet_vols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate skewness of the mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE comp_tet_skew()
    !counting indeces
    INTEGER :: i
    !side length of regular tet in the circumsphere
    REAL(8) :: a_side
    !volume of the regular tet in the circumsphere
    REAL(8) :: vol_reg

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=0

    !compute tet skews
    DO i=1,tot_tets
      !compute the side of the regular tet in the circumsphere of the tet
      a_side=4.0D0*tet(i)%sphere_rad()/SQRT(6.0D0)
      !compute the volume of the regular tet in the circumsphere of the tet
      vol_reg=a_side**3/(6.0D0*SQRT(2.0D0))
      !compute the cell skew
      tet(i)%skew=(vol_reg-tet(i)%vol)/vol_reg

      !cell skew average
      tot_mesh%skew_avg=tot_mesh%skew_avg+tet(i)%skew
      reg_mesh(tet(i)%reg)%skew_avg=reg_mesh(tet(i)%reg)%skew_avg+tet(i)%skew
      IF(MOD(i,CEILING(tot_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
    tot_mesh%skew_avg=tot_mesh%skew_avg/(tot_tets*1.0D0)
    DO i=minreg,maxreg
      IF(reg_mesh(i)%num_el .GT. 0)reg_mesh(i)%skew_avg=reg_mesh(i)%skew_avg/ &
          (reg_mesh(i)%num_el*1.0D0)
    ENDDO

    !compute skew standard deviation
    DO i=1,tot_tets
      tot_mesh%skew_sd=tot_mesh%skew_sd+(tet(i)%skew-tot_mesh%skew_avg)**2
      reg_mesh(tet(i)%reg)%skew_sd=reg_mesh(tet(i)%reg)%skew_sd+(tet(i)%skew &
          -reg_mesh(tet(i)%reg)%skew_avg)**2
    ENDDO
    tot_mesh%skew_sd=SQRT(tot_mesh%skew_sd/(tot_tets-1.0D0))
    DO i=minreg,maxreg
      IF(reg_mesh(i)%num_el .GT. 0)reg_mesh(i)%skew_sd=SQRT(reg_mesh(i)%skew_sd/ &
          (reg_mesh(i)%num_el-1.0D0))
    ENDDO

    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)
  ENDSUBROUTINE comp_tet_skew

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !calculate aspect ratio of the mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE comp_tet_ar()
      !counting indeces
      INTEGER :: i
      !length of the sides of a given tet
      REAL(8) :: llen(6)

      WRITE(*,'(A)',ADVANCE='NO')'Progress:'
      prog=0

      !compute tet ars
      DO i=1,tot_tets
        !compute the length of each of the six sides
        llen=0
        llen(1)=tet(i)%corner(1)%p%distance(tet(i)%corner(2)%p)
        llen(2)=tet(i)%corner(1)%p%distance(tet(i)%corner(3)%p)
        llen(3)=tet(i)%corner(1)%p%distance(tet(i)%corner(4)%p)
        llen(4)=tet(i)%corner(2)%p%distance(tet(i)%corner(3)%p)
        llen(5)=tet(i)%corner(2)%p%distance(tet(i)%corner(4)%p)
        llen(6)=tet(i)%corner(3)%p%distance(tet(i)%corner(4)%p)
        !compute the cell ar
        tet(i)%aspect_ratio=MAXVAL(llen)/MINVAL(llen)

        !cell ar average
        tot_mesh%ar_avg=tot_mesh%ar_avg+tet(i)%aspect_ratio
        reg_mesh(tet(i)%reg)%ar_avg=reg_mesh(tet(i)%reg)%ar_avg+tet(i)%aspect_ratio
        IF(MOD(i,CEILING(tot_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
          WRITE(*,'(A)',ADVANCE='NO')'*'
          prog=prog+1
        ENDIF
      ENDDO
      tot_mesh%ar_avg=tot_mesh%ar_avg/(tot_tets*1.0D0)
      DO i=minreg,maxreg
        IF(reg_mesh(i)%num_el .GT. 0)reg_mesh(i)%ar_avg=reg_mesh(i)%ar_avg/ &
            (reg_mesh(i)%num_el*1.0D0)
      ENDDO

      !compute ar standard deviation
      DO i=1,tot_tets
        tot_mesh%ar_sd=tot_mesh%ar_sd+(tet(i)%aspect_ratio-tot_mesh%ar_avg)**2
        reg_mesh(tet(i)%reg)%ar_sd=reg_mesh(tet(i)%reg)%ar_sd+(tet(i)%aspect_ratio- &
            reg_mesh(tet(i)%reg)%ar_avg)**2
      ENDDO
      tot_mesh%ar_sd=SQRT(tot_mesh%ar_sd/(tot_tets-1.0D0))
      DO i=minreg,maxreg
        IF(reg_mesh(i)%num_el .GT. 0)reg_mesh(i)%ar_sd=SQRT(reg_mesh(i)%ar_sd/ &
            (reg_mesh(i)%num_el-1.0D0))
      ENDDO

      DO i=prog,max_prog
        WRITE(*,'(A)',ADVANCE='NO')'*'
      ENDDO
      WRITE(*,*)
    ENDSUBROUTINE comp_tet_ar
END MODULE tet_analyze
