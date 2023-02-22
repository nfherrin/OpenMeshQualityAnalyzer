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
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_vol(i),', '
      WRITE(30,'(I0,A)',ADVANCE='NO')tets_in_reg(i),', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_vol(i)/(tets_in_reg(i)*1.0D0),', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_vol_sd(i),', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_avg_skew(i),', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_sd_skew(i),', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_avg_ar(i),', '
      WRITE(30,'(ES14.8,A)',ADVANCE='NO')reg_sd_ar(i),', '
      WRITE(30,*)
    ENDDO
    WRITE(30,'(A,A)',ADVANCE='NO')'total, '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_vol,', '
    WRITE(30,'(I0,A)',ADVANCE='NO')num_tets,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_vol/(num_tets*1.0D0),', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_vol_sd,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_avg_skew,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_sd_skew,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_avg_ar,', '
    WRITE(30,'(ES14.8,A)',ADVANCE='NO')tot_sd_ar,', '
    CLOSE(30)
  ENDSUBROUTINE output_statistics
END MODULE out_stats
