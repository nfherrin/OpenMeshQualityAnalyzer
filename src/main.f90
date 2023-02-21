!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
! OpenMeshQualityAnalyzer to go from gmsh to thrm
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
PROGRAM openmeshqualityanalyzer
  USE globals
  USE read_gmsh
  USE mesh_analyze
  USE boundary_conditions
  IMPLICIT NONE

  !> String for temporarily holding input data for validation prior to storing
  !! in globals
  CHARACTER(200) :: temp_string = ""

  !> The number of provided command line arguments
  INTEGER :: arg_count,i,j

  !Get number of command line args and check for proper invocation
  arg_count = COMMAND_ARGUMENT_COUNT()

  IF (arg_count .NE. 1)STOP 'must give mesh file and mesh file ONLY'

  !get mesh in file name
  CALL GET_COMMAND_ARGUMENT(1, mesh_infile)

  WRITE(*,'(A)')'----------------------- Reading in mesh:'
  CALL read_gmsh_file()

  WRITE(*,'(A)')'----------------------- Calculating adjacencies:'
  CALL adjacency_calc()

  WRITE(*,'(A)')'----------------------- Calculating volumes:'
  CALL calcvols()

  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'---------------------- OpenMeshQualityAnalyzer successful ----------------------'
ENDPROGRAM openmeshqualityanalyzer
