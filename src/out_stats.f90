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
    INTEGER :: i

    OPEN(UNIT=30,FILE=TRIM(ADJUSTL(mesh_infile))//'_stats.csv',ACTION='WRITE',STATUS='REPLACE')

    WRITE(30,'(A)')"Region, Vol, Tets, Avg Tet Vol, Tet Vol SD, Avg Skew, Skew SD, Avg AR, AR SD,"

    DO i=minreg,maxreg
      WRITE(30,'(I0,A)',ADVANCE='NO')i,', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%vol,', '
      WRITE(30,'(I0,A)',ADVANCE='NO')reg_mesh(i)%num_el,', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%vol_avg,', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%vol_sd,', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%skew_avg,', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%skew_sd,', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%ar_avg,', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_mesh(i)%ar_sd,', '
      WRITE(30,*)
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
    CLOSE(30)
  ENDSUBROUTINE output_statistics
END MODULE out_stats
