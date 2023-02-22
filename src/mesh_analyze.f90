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

  !tet volumes
  REAL(8), ALLOCATABLE :: tetvol(:)
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
    write(*,*)'minreg,maxreg',minreg,maxreg
    ALLOCATE(tetvol(num_tets),reg_vol(minreg:maxreg),tets_in_reg(minreg:maxreg),&
        reg_vol_sd(minreg:maxreg))
    tets_in_reg=0
    tetvol=0
    tot_vol=0
    reg_vol=0
    prog=0

    !compute tet volumes and add to both total volumes and region volumes
    DO i=1,num_tets
      a=tet(i)%corner(1)%p
      b=tet(i)%corner(2)%p
      c=tet(i)%corner(3)%p
      d=tet(i)%corner(4)%p
      tetvol(i)=ABS((-c%y*d%x+b%y*(-c%x+d%x)+b%x*(c%y-d%y)+c%x*d%y)*(a%z-d%z)+(a%x-d%x) &
          *(-c%z*d%y+b%z*(-c%y+d%y)+b%y*(c%z-d%z)+c%y*d%z)+(a%y-d%y)*(b%z*(c%x-d%x) &
          +c%z*d%x-c%x*d%z+b%x*(-c%z+d%z)))/6
      reg_vol(tet(i)%reg)=reg_vol(tet(i)%reg)+tetvol(i)
      tets_in_reg(tet(i)%reg)=tets_in_reg(tet(i)%reg)+1
      tot_vol=tot_vol+tetvol(i)
      IF(MOD(i,CEILING(num_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO

    !compute volume standard deviations
    tot_vol_sd=0.0
    reg_vol_sd=0.0
    DO i=1,num_tets
      tot_vol_sd=tot_vol_sd+(tetvol(i)-tot_vol/(num_tets*1.0D0))**2
      IF(tets_in_reg(tet(i)%reg) .NE. 0)THEN
        reg_vol_sd(tet(i)%reg)=reg_vol_sd(tet(i)%reg)+ &
            (tetvol(i)-reg_vol(tet(i)%reg)/(tets_in_reg(tet(i)%reg)*1.0D0))**2
      ENDIF
    ENDDO
    tot_vol_sd=SQRT(tot_vol_sd/(num_tets-1.0D0))
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
    !how skewed each tet is
    REAL(8) :: cell_skew(num_tets)
    !counting indeces
    INTEGER :: i
    !side length of regular tet in the circumsphere
    REAL(8) :: a_side
    !volume of the regular tet in the circumsphere
    REAL(8) :: vol_reg

    ALLOCATE(reg_avg_skew(minreg:maxreg),reg_sd_skew(minreg:maxreg))

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=0

    cell_skew=0.0
    reg_avg_skew=0.0
    reg_sd_skew=0.0
    !compute tet skews
    DO i=1,num_tets
      !compute the side of the regular tet in the circumsphere of the tet
      a_side=4.0D0*sphere_rad(tet(i)%corner(1)%p,tet(i)%corner(2)%p,tet(i)%corner(3)%p, &
          tet(i)%corner(4)%p)/SQRT(6.0D0)
      !compute the volume of the regular tet in the circumsphere of the tet
      vol_reg=a_side**3/(6.0D0*SQRT(2.0D0))
      !compute the cell skew
      cell_skew(i)=(vol_reg-tetvol(i))/vol_reg

      !cell skew average
      tot_avg_skew=tot_avg_skew+cell_skew(i)
      reg_avg_skew(tet(i)%reg)=reg_avg_skew(tet(i)%reg)+cell_skew(i)
      IF(MOD(i,CEILING(num_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
    tot_avg_skew=tot_avg_skew/(num_tets*1.0D0)
    DO i=minreg,maxreg
      IF(tets_in_reg(i) .GT. 0)reg_avg_skew(i)=reg_avg_skew(i)/(tets_in_reg(i)*1.0D0)
    ENDDO

    !compute skew standard deviation
    DO i=1,num_tets
      tot_sd_skew=tot_sd_skew+(cell_skew(i)-tot_avg_skew)**2
      reg_sd_skew(tet(i)%reg)=reg_sd_skew(tet(i)%reg)+(cell_skew(i)-reg_avg_skew(tet(i)%reg))**2
    ENDDO
    tot_sd_skew=SQRT(tot_sd_skew/(num_tets-1.0D0))
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
      !aspect ratio of each tet
      REAL(8) :: cell_ar(num_tets)
      !counting indeces
      INTEGER :: i
      !length of the sides of a given tet
      REAL(8) :: llen(6)

      ALLOCATE(reg_avg_ar(minreg:maxreg),reg_sd_ar(minreg:maxreg))

      WRITE(*,'(A)',ADVANCE='NO')'Progress:'
      prog=0

      cell_ar=0.0
      reg_avg_ar=0.0
      reg_sd_ar=0.0
      !compute tet ars
      DO i=1,num_tets
        !compute the length of each of the six sides
        llen=0
        llen(1)=line_length(tet(i)%corner(1)%p,tet(i)%corner(2)%p)
        llen(2)=line_length(tet(i)%corner(1)%p,tet(i)%corner(3)%p)
        llen(3)=line_length(tet(i)%corner(1)%p,tet(i)%corner(4)%p)
        llen(4)=line_length(tet(i)%corner(2)%p,tet(i)%corner(3)%p)
        llen(5)=line_length(tet(i)%corner(2)%p,tet(i)%corner(4)%p)
        llen(6)=line_length(tet(i)%corner(3)%p,tet(i)%corner(4)%p)
        !compute the cell ar
        cell_ar(i)=MAXVAL(llen)/MINVAL(llen)

        !cell ar average
        tot_avg_ar=tot_avg_ar+cell_ar(i)
        reg_avg_ar(tet(i)%reg)=reg_avg_ar(tet(i)%reg)+cell_ar(i)
        IF(MOD(i,CEILING(num_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
          WRITE(*,'(A)',ADVANCE='NO')'*'
          prog=prog+1
        ENDIF
      ENDDO
      tot_avg_ar=tot_avg_ar/(num_tets*1.0D0)
      DO i=minreg,maxreg
        IF(tets_in_reg(i) .GT. 0)reg_avg_ar(i)=reg_avg_ar(i)/(tets_in_reg(i)*1.0D0)
      ENDDO

      !compute ar standard deviation
      DO i=1,num_tets
        tot_sd_ar=tot_sd_ar+(cell_ar(i)-tot_avg_ar)**2
        reg_sd_ar(tet(i)%reg)=reg_sd_ar(tet(i)%reg)+(cell_ar(i)-reg_avg_ar(tet(i)%reg))**2
      ENDDO
      tot_sd_ar=SQRT(tot_sd_ar/(num_tets-1.0D0))
      DO i=minreg,maxreg
        IF(tets_in_reg(i) .GT. 0)reg_sd_ar(i)=SQRT(reg_sd_ar(i)/(tets_in_reg(i)-1.0D0))
      ENDDO

      DO i=prog,max_prog
        WRITE(*,'(A)',ADVANCE='NO')'*'
      ENDDO
      WRITE(*,*)
    ENDSUBROUTINE comp_ar

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !computes line length between two points
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(8) FUNCTION sphere_rad(a,b,c,d)
      TYPE(vertex_type), INTENT(IN) :: a,c,b,d

      REAL(8) :: ba(3),ca(3),da(3)
      REAL(8) :: len_ba,len_ca,len_da
      REAL(8) :: cross_cd(3),cross_db(3),cross_bc(3)
      REAL(8) :: denominator
      REAL(8) :: circ(3)

      ba(1)=b%x-a%x
      ba(2)=b%y-a%y
      ba(3)=b%z-a%z

      ca(1)=c%x-a%x
      ca(2)=c%y-a%y
      ca(3)=c%z-a%z

      da(1)=d%x-a%x
      da(2)=d%y-a%y
      da(3)=d%z-a%z

      len_ba=ba(1)*ba(1)+ba(2)*ba(2)+ba(3)*ba(3);
      len_ca=ca(1)*ca(1)+ca(2)*ca(2)+ca(3)*ca(3);
      len_da=da(1)*da(1)+da(2)*da(2)+da(3)*da(3);

      cross_cd(1)=ca(2)*da(3)-da(2)*ca(3);
      cross_cd(2)=ca(3)*da(1)-da(3)*ca(1);
      cross_cd(3)=ca(1)*da(2)-da(1)*ca(2);

      cross_db(1)=da(2)*ba(3)-ba(2)*da(3);
      cross_db(2)=da(3)*ba(1)-ba(3)*da(1);
      cross_db(3)=da(1)*ba(2)-ba(1)*da(2);

      cross_bc(1)=ba(2)*ca(3)-ca(2)*ba(3);
      cross_bc(2)=ba(3)*ca(1)-ca(3)*ba(1);
      cross_bc(3)=ba(1)*ca(2)-ca(1)*ba(2);

      denominator=0.5D0/(ba(1)*cross_cd(1)+ba(2)*cross_cd(2)+ba(3)*cross_cd(3));

      circ(1)=(len_ba*cross_cd(1)+len_ca*cross_db(1)+len_da*cross_bc(1))*denominator
      circ(2)=(len_ba*cross_cd(2)+len_ca*cross_db(2)+len_da*cross_bc(2))*denominator
      circ(3)=(len_ba*cross_cd(3)+len_ca*cross_db(3)+len_da*cross_bc(3))*denominator

      sphere_rad=SQRT(circ(1)**2+circ(2)**2+circ(3)**2)
    ENDFUNCTION sphere_rad

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !computes the length of the line between two vertices
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL(8) FUNCTION line_length(a,b)
        TYPE(vertex_type), INTENT(IN) :: a,b

        line_length=SQRT((a%x-b%x)**2+(a%y-b%y)**2+(a%z-b%z)**2)
      ENDFUNCTION line_length
END MODULE mesh_analyze
