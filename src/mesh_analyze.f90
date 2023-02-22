!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to analyze mesh quality indicators
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE mesh_analyze
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calcvols,comp_skew,comp_ar
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate volumes of each tet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calcvols()
    TYPE(vertex_type) :: a,b,c,d
    INTEGER :: i

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'

    minreg=MINVAL(tet(:)%reg)
    maxreg=MAXVAL(tet(:)%reg)
    ALLOCATE(reg_vol(minreg:maxreg),tets_in_reg(minreg:maxreg),&
        reg_vol_sd(minreg:maxreg))
    tets_in_reg=0
    tot_vol=0
    reg_vol=0
    prog=0

    !compute tet volumes and add to both total volumes and region volumes
    DO i=1,tot_tets
      a=tet(i)%corner(1)%p
      b=tet(i)%corner(2)%p
      c=tet(i)%corner(3)%p
      d=tet(i)%corner(4)%p
      tet(i)%vol=ABS((-c%y*d%x+b%y*(-c%x+d%x)+b%x*(c%y-d%y)+c%x*d%y)*(a%z-d%z)+(a%x-d%x) &
          *(-c%z*d%y+b%z*(-c%y+d%y)+b%y*(c%z-d%z)+c%y*d%z)+(a%y-d%y)*(b%z*(c%x-d%x) &
          +c%z*d%x-c%x*d%z+b%x*(-c%z+d%z)))/6
      reg_vol(tet(i)%reg)=reg_vol(tet(i)%reg)+tet(i)%vol
      tets_in_reg(tet(i)%reg)=tets_in_reg(tet(i)%reg)+1
      tot_vol=tot_vol+tet(i)%vol
      IF(MOD(i,CEILING(tot_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO

    !compute volume standard deviations
    tot_vol_sd=0.0
    reg_vol_sd=0.0
    DO i=1,tot_tets
      tot_vol_sd=tot_vol_sd+(tet(i)%vol-tot_vol/(tot_tets*1.0D0))**2
      IF(tets_in_reg(tet(i)%reg) .NE. 0)THEN
        reg_vol_sd(tet(i)%reg)=reg_vol_sd(tet(i)%reg)+ &
            (tet(i)%vol-reg_vol(tet(i)%reg)/(tets_in_reg(tet(i)%reg)*1.0D0))**2
      ENDIF
    ENDDO
    tot_vol_sd=SQRT(tot_vol_sd/(tot_tets-1.0D0))
    DO i=minreg,maxreg
      IF(tets_in_reg(i) .NE. 0)reg_vol_sd(i)=SQRT(reg_vol_sd(i)/(tets_in_reg(i)-1.0D0))
    ENDDO

    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)
  ENDSUBROUTINE calcvols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate skewness of the mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE comp_skew()
    !counting indeces
    INTEGER :: i
    !side length of regular tet in the circumsphere
    REAL(8) :: a_side
    !volume of the regular tet in the circumsphere
    REAL(8) :: vol_reg

    ALLOCATE(reg_avg_skew(minreg:maxreg),reg_sd_skew(minreg:maxreg))

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=0

    reg_avg_skew=0.0
    reg_sd_skew=0.0
    !compute tet skews
    DO i=1,tot_tets
      !compute the side of the regular tet in the circumsphere of the tet
      a_side=4.0D0*tet(i)%sphere_rad()/SQRT(6.0D0)
      !compute the volume of the regular tet in the circumsphere of the tet
      vol_reg=a_side**3/(6.0D0*SQRT(2.0D0))
      !compute the cell skew
      tet(i)%skew=(vol_reg-tet(i)%vol)/vol_reg

      !cell skew average
      tot_avg_skew=tot_avg_skew+tet(i)%skew
      reg_avg_skew(tet(i)%reg)=reg_avg_skew(tet(i)%reg)+tet(i)%skew
      IF(MOD(i,CEILING(tot_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
    tot_avg_skew=tot_avg_skew/(tot_tets*1.0D0)
    DO i=minreg,maxreg
      IF(tets_in_reg(i) .GT. 0)reg_avg_skew(i)=reg_avg_skew(i)/(tets_in_reg(i)*1.0D0)
    ENDDO

    !compute skew standard deviation
    DO i=1,tot_tets
      tot_sd_skew=tot_sd_skew+(tet(i)%skew-tot_avg_skew)**2
      reg_sd_skew(tet(i)%reg)=reg_sd_skew(tet(i)%reg)+(tet(i)%skew-reg_avg_skew(tet(i)%reg))**2
    ENDDO
    tot_sd_skew=SQRT(tot_sd_skew/(tot_tets-1.0D0))
    DO i=minreg,maxreg
      IF(tets_in_reg(i) .GT. 0)reg_sd_skew(i)=SQRT(reg_sd_skew(i)/(tets_in_reg(i)-1.0D0))
    ENDDO

    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)
  ENDSUBROUTINE comp_skew

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !calculate aspect ratio of the mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE comp_ar()
      !counting indeces
      INTEGER :: i
      !length of the sides of a given tet
      REAL(8) :: llen(6)

      ALLOCATE(reg_avg_ar(minreg:maxreg),reg_sd_ar(minreg:maxreg))

      WRITE(*,'(A)',ADVANCE='NO')'Progress:'
      prog=0

      reg_avg_ar=0.0
      reg_sd_ar=0.0
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
        tot_avg_ar=tot_avg_ar+tet(i)%aspect_ratio
        reg_avg_ar(tet(i)%reg)=reg_avg_ar(tet(i)%reg)+tet(i)%aspect_ratio
        IF(MOD(i,CEILING(tot_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
          WRITE(*,'(A)',ADVANCE='NO')'*'
          prog=prog+1
        ENDIF
      ENDDO
      tot_avg_ar=tot_avg_ar/(tot_tets*1.0D0)
      DO i=minreg,maxreg
        IF(tets_in_reg(i) .GT. 0)reg_avg_ar(i)=reg_avg_ar(i)/(tets_in_reg(i)*1.0D0)
      ENDDO

      !compute ar standard deviation
      DO i=1,tot_tets
        tot_sd_ar=tot_sd_ar+(tet(i)%aspect_ratio-tot_avg_ar)**2
        reg_sd_ar(tet(i)%reg)=reg_sd_ar(tet(i)%reg)+(tet(i)%aspect_ratio-reg_avg_ar(tet(i)%reg))**2
      ENDDO
      tot_sd_ar=SQRT(tot_sd_ar/(tot_tets-1.0D0))
      DO i=minreg,maxreg
        IF(tets_in_reg(i) .GT. 0)reg_sd_ar(i)=SQRT(reg_sd_ar(i)/(tets_in_reg(i)-1.0D0))
      ENDDO

      DO i=prog,max_prog
        WRITE(*,'(A)',ADVANCE='NO')'*'
      ENDDO
      WRITE(*,*)
    ENDSUBROUTINE comp_ar
END MODULE mesh_analyze
