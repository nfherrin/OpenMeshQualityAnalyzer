!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to output mesh statistics
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE out_stats
  USE globals
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: output_statistics
CONTAINS

  !outputs mesh statistics
  SUBROUTINE output_statistics()
    INTEGER :: i,j

    OPEN(UNIT=30,FILE=TRIM(ADJUSTL(mesh_infile))//'_stats.csv',ACTION='WRITE',STATUS='REPLACE')

    WRITE(30,'(A)')'3D Mesh Data,'
    WRITE(30,'(A)')"Region, Vol, Tets, Avg Tet Vol, Tet Vol SD, Avg Skew, Skew SD, Avg AR, AR SD, &
        &Avg Smooth, Smooth SD,"

    DO i=minreg,maxreg
      IF(reg_mesh(i)%num_el .GT. 0)THEN
        WRITE(30,'(I0,A)',ADVANCE='NO')i,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%vol,', '
        WRITE(30,'(I0,A)',ADVANCE='NO')reg_mesh(i)%num_el,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%vol_avg,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%vol_sd,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%skew_avg,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%skew_sd,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%ar_avg,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%ar_sd,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%smooth_avg,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%smooth_sd,', '
        WRITE(30,*)
      ENDIF
    ENDDO
    WRITE(30,'(A,A)',ADVANCE='NO')'total, '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_mesh%vol,', '
    WRITE(30,'(I0,A)',ADVANCE='NO')tot_mesh%num_el,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_mesh%vol_avg,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_mesh%vol_sd,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_mesh%skew_avg,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_mesh%skew_sd,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_mesh%ar_avg,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_mesh%ar_sd,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_mesh%smooth_avg,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_mesh%smooth_sd,', '
    WRITE(30,*)
    WRITE(30,'(A)')', '
    WRITE(30,'(A)')', '
    WRITE(30,'(A)')', '
    WRITE(30,'(A)')'2D Mesh Data,'
    DO i=1,6
      IF(side_flat(i))THEN
        SELECTCASE(i)
          CASE(1)
            WRITE(30,'(A)')'-x Boundary Mesh,'
          CASE(2)
            WRITE(30,'(A)')'+x boundary,'
          CASE(3)
            WRITE(30,'(A)')'-y Boundary Mesh,'
          CASE(4)
            WRITE(30,'(A)')'+y Boundary Mesh,'
          CASE(5)
            WRITE(30,'(A)')'-z Boundary Mesh,'
          CASE(6)
            WRITE(30,'(A)')'+z Boundary Mesh,'
        ENDSELECT
        WRITE(30,'(A)')"Region, Area, Tris, Avg Tri Area, Tri Area SD, Avg Skew, Skew SD, Avg AR, &
            &AR SD, Avg Smooth, Smooth SD,"
        DO j=minreg,maxreg
          IF(reg_side_mesh(i,j)%num_el .GT. 0)THEN
            WRITE(30,'(I0,A)',ADVANCE='NO')j,', '
            WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_side_mesh(i,j)%area,', '
            WRITE(30,'(I0,A)',ADVANCE='NO')reg_side_mesh(i,j)%num_el,', '
            WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_side_mesh(i,j)%area_avg,', '
            WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_side_mesh(i,j)%area_sd,', '
            WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_side_mesh(i,j)%skew_avg,', '
            WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_side_mesh(i,j)%skew_sd,', '
            WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_side_mesh(i,j)%ar_avg,', '
            WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_side_mesh(i,j)%ar_sd,', '
            WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_side_mesh(i,j)%smooth_avg,', '
            WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_side_mesh(i,j)%smooth_sd,', '
            WRITE(30,*)
          ENDIF
        ENDDO
        WRITE(30,'(A,A)',ADVANCE='NO')'total, '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_side_mesh(i)%area,', '
        WRITE(30,'(I0,A)',ADVANCE='NO')tot_side_mesh(i)%num_el,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_side_mesh(i)%area_avg,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_side_mesh(i)%area_sd,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_side_mesh(i)%skew_avg,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_side_mesh(i)%skew_sd,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_side_mesh(i)%ar_avg,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_side_mesh(i)%ar_sd,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_side_mesh(i)%smooth_avg,', '
        WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_side_mesh(i)%smooth_sd,', '
        WRITE(30,*)
      ENDIF
    ENDDO
    CLOSE(30)
  ENDSUBROUTINE output_statistics
END MODULE out_stats
