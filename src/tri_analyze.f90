!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to analyze tri quality indicators
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE tri_analyze
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_tri_areas, comp_tri_skew, comp_tri_ar, comp_tri_smooth
CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate areas of each tri
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calc_tri_areas()
    INTEGER :: i
    REAL(8) :: ll(3),la,lb,lc

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=0
    DO i=1,tot_bcf

      ll(1)=tri(i)%corner(1)%p%distance(tri(i)%corner(2)%p)
      ll(2)=tri(i)%corner(1)%p%distance(tri(i)%corner(3)%p)
      ll(3)=tri(i)%corner(2)%p%distance(tri(i)%corner(3)%p)
      !bubble sort to make sure our orders are good for numerical stability
      CALL bubble_sort(ll(:))
      la=ll(3)
      lb=ll(2)
      lc=ll(1)
      tri(i)%area=SQRT((la+(lb+lc))*(lc-(la-lb))*(lc+(la-lb))*(la+(lb-lc)))/4.0D0

      IF(MOD(i,CEILING(tot_bcf*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
  ENDSUBROUTINE calc_tri_areas

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate skew of each tri
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE comp_tri_skew()
    !counting indeces
    INTEGER :: i
    !side length of regular tri in the circumcircle
    REAL(8) :: a_side
    !area of the regular tri in the circumcircle
    REAL(8) :: area_reg

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=0

    !compute tri skews
    DO i=1,tot_bcf
      !compute the side of the regular tri in the circumcircle of the tri
      a_side=tri(i)%circ_rad()*SQRT(3.0D0)
      !compute the area of the regular tri in the circumcircle of the tri
      area_reg=SQRT(3.0D0)*a_side**2/4.0D0
      !compute the cell skew
      tri(i)%skew=(area_reg-tri(i)%area)/area_reg

      !cell skew average
      IF(MOD(i,CEILING(tot_bcf*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
  ENDSUBROUTINE comp_tri_skew

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !calculate aspect ratio of the mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE comp_tri_ar()
    !counting indeces
    INTEGER :: i
    !length of the sides of a given tri
    REAL(8) :: llen(3)

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=0

    !compute tri ars
    DO i=1,tot_bcf
      !compute the length of each of the six sides
      llen=0
      llen(1)=tri(i)%corner(1)%p%distance(tri(i)%corner(2)%p)
      llen(2)=tri(i)%corner(1)%p%distance(tri(i)%corner(3)%p)
      llen(3)=tri(i)%corner(2)%p%distance(tri(i)%corner(3)%p)
      !compute the cell ar
      tri(i)%aspect_ratio=MAXVAL(llen)/MINVAL(llen)

      IF(MOD(i,CEILING(tot_bcf*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
  ENDSUBROUTINE comp_tri_ar

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !calculate the smoothness of a given tri (area compared to surrounding tris)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE comp_tri_smooth()
    INTEGER :: i,j,adj_i
    REAL(8) :: num_adj_tris

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=0

    DO i=1,tot_bcf
      num_adj_tris=0
      tri(i)%smoothness=0
      !compute summed smoothness related to adjacent cells
      DO j=1,3
        adj_i=tri(i)%adj_id(j)
        IF(adj_i .NE. 0)THEN
          num_adj_tris=num_adj_tris+1.0
          !add to smoothness sum, which is larger area divided by smaller area
          IF(tri(adj_i)%area .GE. tri(i)%area)THEN
            tri(i)%smoothness=tri(i)%smoothness+tri(adj_i)%area/tri(i)%area
          ELSE
            tri(i)%smoothness=tri(i)%smoothness+tri(i)%area/tri(adj_i)%area
          ENDIF
        ENDIF
      ENDDO
      IF(tri(i)%smoothness .LE. 0)STOP 'tri cannot be solo!'
      !average smoothness
      tri(i)%smoothness=tri(i)%smoothness/num_adj_tris
      IF(MOD(i,CEILING(tot_bcf*1.0D0/(max_prog-1.0D0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
  ENDSUBROUTINE comp_tri_smooth
END MODULE tri_analyze
