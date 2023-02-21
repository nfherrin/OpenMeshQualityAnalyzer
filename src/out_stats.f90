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

    WRITE(30,'(A)')"Region, Vol, Tets, Avg Tet Vol, Min Skew, Max Skew, Avg Skew, Skew SD"

    DO i=minreg,maxreg
      WRITE(30,*)i
    ENDDO
    CLOSE(30)
  ENDSUBROUTINE output_statistics
END MODULE out_stats
