!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to ingest a file in the
!! THOR mesh format (.thrm)
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE read_thrm
  USE globals
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_thrm_file
CONTAINS

  !reads in a gmsh file
  SUBROUTINE read_thrm_file()
    STOP 'read_thrm_file not yet supported'
  ENDSUBROUTINE read_thrm_file
END MODULE read_thrm
