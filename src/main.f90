!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
! OpenMeshQualityAnalyzer to go from gmsh to thrm
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
PROGRAM openmeshqualityanalyzer
  USE globals
  USE read_gmsh
  USE read_thrm
  USE mesh_analyze
  USE boundary_conditions
  IMPLICIT NONE

  !> The number of provided command line arguments
  INTEGER :: arg_count

  !type of mesh format
  CHARACTER(20) :: mesh_format="",tchar1,tchar2

  !Get number of command line args and check for proper invocation
  arg_count = COMMAND_ARGUMENT_COUNT()
  IF (arg_count .NE. 1)STOP 'must give mesh file and mesh file ONLY'

  !get mesh in file name
  CALL GET_COMMAND_ARGUMENT(1, mesh_infile)

  !check if file type is gmsh or thrm
  OPEN(unit=20,FILE=mesh_infile,ACTION='READ',STATUS='OLD')

  READ(20,*)tchar1
  READ(20,*)tchar2
  IF(tchar1 .EQ. "$MeshFormat" .AND. tchar2 .EQ. "4.1")mesh_format="gmsh"
  IF(mesh_format .EQ. "")THEN
    READ(20,*)tchar1
    READ(20,*)tchar2
    IF(tchar1 .EQ. "1" .AND. tchar2 .EQ. "1")mesh_format="thrm"
  ENDIF

  IF(mesh_format .EQ. "") STOP "Mesh format not recognized, currently on gmsh and thrm supported"

  CLOSE(20)

  !read in the mesh and calculate adjacencies if needed
  WRITE(*,'(A)')'----------------------- Reading in mesh:'
  SELECTCASE(mesh_format)
    CASE("gmsh")
      CALL read_gmsh_file()

      WRITE(*,'(A)')'----------------------- Calculating adjacencies:'
      CALL adjacency_calc()
    CASE("thrm")
      CALL read_thrm_file()
    CASE DEFAULT
      STOP "Mesh format not yet supported"
  ENDSELECT

  WRITE(*,'(A)')'----------------------- Calculating volumes:'
  CALL calcvols()

  WRITE(*,'(A)')'----------------------- Computing mesh skewness:'
  CALL comp_skew()

  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'---------------------- OpenMeshQualityAnalyzer successful ----------------------'
ENDPROGRAM openmeshqualityanalyzer
