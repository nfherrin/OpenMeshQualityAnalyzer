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
    REAL(8) :: totalvol1
    REAL(8) :: a(3),b(3),c(3),d(3)
    REAL(8), ALLOCATABLE :: regvol(:)
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
      IF(MOD(i,CEILING(num_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
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
      WRITE(*,'(A,I0,A,ES24.16)')'Region ',i,' equivalent radius: ', &
          (3.0/4.0/pi*regvol(i))**(1.0/3.0)
    ENDDO
    WRITE(*,'(A,I0)')'Total number of tets: ',SUM(tets_in_reg)
    WRITE(*,'(A,ES24.16)')'Total system volume: ',totalvol1
    WRITE(*,'(A,ES24.16)')'Equivalent radius: ',(3.0/4.0/pi*totalvol1)**(1.0/3.0)
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
    !Total mesh average cell skew
    REAL(8) :: tot_avg_skew=0.0D0
    !Total mesh skew standard deviation
    REAL(8) :: tot_sd_kew=0.0D0

    cell_skew=0.0
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
    ENDDO
    tot_avg_skew=tot_avg_skew/(num_tets*1.0D0)
    !compute skew standard deviation
    DO i=1,num_tets
      tot_sd_kew=tot_sd_kew+(cell_skew(i)-tot_avg_skew)**2
    ENDDO
    tot_sd_kew=SQRT(tot_sd_kew/(num_tets-1.0D0))
    WRITE(*,'(A,ES16.8)')'Mesh min skew: ',MINVAL(cell_skew)
    WRITE(*,'(A,ES16.8)')'Mesh max skew: ',MAXVAL(cell_skew)
    WRITE(*,'(A,ES16.8)')'Mesh average skew: ',tot_avg_skew
    WRITE(*,'(A,ES16.8)')'Mesh skew standard deviation: ',tot_sd_kew
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
