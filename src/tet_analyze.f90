!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to analyze tet quality indicators
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE tet_analyze
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_tet_vols,comp_tet_skew,comp_tet_ar,comp_tet_smooth
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
      IF(MOD(i,CEILING(tot_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
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

      IF(MOD(i,CEILING(tot_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
  ENDSUBROUTINE comp_tet_ar

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate the smoothness of a given tet (volume compared to surrounding tets)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE comp_tet_smooth()
    INTEGER :: i,j,adj_i
    REAL(8) :: num_adj_tets

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=0

    DO i=1,tot_tets
      num_adj_tets=0
      tet(i)%smoothness=0
      !compute summed smoothness related to adjacent cells
      DO j=1,4
        adj_i=tet(i)%adj_id(j)
        IF(adj_i .NE. 0)THEN
          num_adj_tets=num_adj_tets+1.0
          !add to smoothness sum, which is larger volume divided by smaller volume
          IF(tet(adj_i)%vol .GE. tet(i)%vol)THEN
            tet(i)%smoothness=tet(i)%smoothness+tet(adj_i)%vol/tet(i)%vol
          ELSE
            tet(i)%smoothness=tet(i)%smoothness+tet(i)%vol/tet(adj_i)%vol
          ENDIF
        ENDIF
      ENDDO
      IF(tet(i)%smoothness .LE. 0)STOP 'tet cannot be solo!'
      !average smoothness
      tet(i)%smoothness=tet(i)%smoothness/num_adj_tets
      IF(MOD(i,CEILING(tot_tets*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
  ENDSUBROUTINE comp_tet_smooth
END MODULE tet_analyze
