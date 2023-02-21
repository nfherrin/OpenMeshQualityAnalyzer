!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to analyze mesh quality indicators
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE mesh_analyze
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calcvols,comp_skew

  !tet volumes
  REAL(8), ALLOCATABLE :: tetvol(:)
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate volumes of each tet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calcvols()
    REAL(8) :: a(3),b(3),c(3),d(3)
    INTEGER :: i

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'

    minreg=MINVAL(el_tag(:))
    maxreg=MAXVAL(el_tag(:))
    ALLOCATE(tetvol(num_tets),reg_vol(minreg:maxreg),tets_in_reg(minreg:maxreg),&
        reg_vol_sd(minreg:maxreg))
    tets_in_reg=0
    tetvol=0
    tot_vol=0
    reg_vol=0
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
      reg_vol(el_tag(i))=reg_vol(el_tag(i))+tetvol(i)
      tets_in_reg(el_tag(i))=tets_in_reg(el_tag(i))+1
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
      IF(tets_in_reg(el_tag(i)) .NE. 0)THEN
        reg_vol_sd(el_tag(i))=reg_vol_sd(el_tag(i))+ &
            (tetvol(i)-reg_vol(el_tag(i))/(tets_in_reg(el_tag(i))*1.0D0))**2
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
      a_side=4.0D0*sphere_rad(vertex(element(i,1),:),vertex(element(i,2),:),vertex(element(i,3),:), &
          vertex(element(i,4),:))/SQRT(6.0D0)
      !compute the volume of the regular tet in the circumsphere of the tet
      vol_reg=a_side**3/(6.0D0*SQRT(2.0D0))
      !compute the cell skew
      cell_skew(i)=(vol_reg-tetvol(i))/vol_reg

      !cell skew average
      tot_avg_skew=tot_avg_skew+cell_skew(i)
      reg_avg_skew(el_tag(i))=reg_avg_skew(el_tag(i))+cell_skew(i)
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
      reg_sd_skew(el_tag(i))=reg_sd_skew(el_tag(i))+(cell_skew(i)-reg_avg_skew(el_tag(i)))**2
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
  !computes the circumsphere radius for a given tet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(8) FUNCTION sphere_rad(a,b,c,d)
    REAL(8), INTENT(IN) :: a(3),c(3),b(3),d(3)

    REAL(8) :: ba(3),ca(3),da(3)
    REAL(8) :: len_ba,len_ca,len_da
    REAL(8) :: cross_cd(3),cross_db(3),cross_bc(3)
    REAL(8) :: denominator
    REAL(8) :: circ(3)

    ba(1)=b(1)-a(1)
    ba(2)=b(2)-a(2)
    ba(3)=b(3)-a(3)

    ca(1)=c(1)-a(1)
    ca(2)=c(2)-a(2)
    ca(3)=c(3)-a(3)

    da(1)=d(1)-a(1)
    da(2)=d(2)-a(2)
    da(3)=d(3)-a(3)

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
END MODULE mesh_analyze
