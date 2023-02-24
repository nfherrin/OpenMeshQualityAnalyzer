!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
! OpenMeshQualityAnalyzer to go from gmsh to thrm
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
PROGRAM openmeshqualityanalyzer
  USE globals
  USE read_gmsh
  USE out_stats
  USE read_thrm
  USE mesh_analyze_3d
  USE mesh_analyze_2d
  USE tet_analyze
  USE tri_analyze
  USE boundary_conditions
  USE mesh_2d_construct
  IMPLICIT NONE

  !> The number of provided command line arguments
  INTEGER :: arg_count,i,j

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

  IF(mesh_format .EQ. "") STOP "Mesh format not recognized, currently only gmsh and thrm supported"

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
  CALL compute_bc_sides()

  !assign tets to each volume
  CALL assign_tets()

  WRITE(*,'(A)')'----------------------- Calculating 3D volumes:'
  CALL calc_tet_vols()
  !for each mesh structure
  CALL mesh_vol_analysis(tot_mesh)
  DO i=minreg,maxreg
    IF(reg_mesh(i)%num_el .GT. 0)CALL mesh_vol_analysis(reg_mesh(i))
  ENDDO
  !finish up progress bar
  DO i=prog,max_prog
    WRITE(*,'(A)',ADVANCE='NO')'*'
  ENDDO
  WRITE(*,*)

  WRITE(*,'(A)')'----------------------- Computing 3D mesh skewness:'
  CALL comp_tet_skew()
  !for each mesh structure
  CALL mesh_skew_analysis(tot_mesh)
  DO i=minreg,maxreg
    IF(reg_mesh(i)%num_el .GT. 0)CALL mesh_skew_analysis(reg_mesh(i))
  ENDDO
  !finish up progress bar
  DO i=prog,max_prog
    WRITE(*,'(A)',ADVANCE='NO')'*'
  ENDDO
  WRITE(*,*)

  WRITE(*,'(A)')'----------------------- Computing 3D mesh aspect ratio:'
  CALL comp_tet_ar()
  !for each mesh structure
  CALL mesh_ar_analysis(tot_mesh)
  DO i=minreg,maxreg
    IF(reg_mesh(i)%num_el .GT. 0)CALL mesh_ar_analysis(reg_mesh(i))
  ENDDO
  !finish up progress bar
  DO i=prog,max_prog
    WRITE(*,'(A)',ADVANCE='NO')'*'
  ENDDO
  WRITE(*,*)

  WRITE(*,'(A)')'----------------------- Computing 3D mesh smoothness:'
  CALL comp_tet_smooth()
  !for each mesh structure
  CALL mesh_smooth_analysis(tot_mesh)
  DO i=minreg,maxreg
    IF(reg_mesh(i)%num_el .GT. 0)CALL mesh_smooth_analysis(reg_mesh(i))
  ENDDO
  !finish up progress bar
  DO i=prog,max_prog
    WRITE(*,'(A)',ADVANCE='NO')'*'
  ENDDO
  WRITE(*,*)

  WRITE(*,'(A)')'----------------------- Constructing 2D Meshes:'
  CALL generate_bc_triangles()
  ALLOCATE(reg_side_mesh(6,minreg:maxreg))
  DO i=1,6
    CALL construct_bc_mesh(tot_side_mesh(i),reg_side_mesh(i,:),i)
  ENDDO

  WRITE(*,'(A)')'----------------------- Calculating 2D areas:'
  CALL calc_tri_areas()
  !for each mesh structure
  DO i=1,6
    CALL mesh_2d_area_analysis(tot_side_mesh(i))
    DO j=minreg,maxreg
      IF(reg_side_mesh(i,j)%num_el .GT. 0)CALL mesh_2d_area_analysis(reg_side_mesh(i,j))
    ENDDO
  ENDDO
  !finish up progress bar
  DO i=prog,max_prog
    WRITE(*,'(A)',ADVANCE='NO')'*'
  ENDDO
  WRITE(*,*)

  WRITE(*,'(A)')'----------------------- Calculating 2D mesh skewness:'
  CALL comp_tri_skew()
  !for each mesh structure
  DO i=1,6
    CALL mesh_2d_skew_analysis(tot_side_mesh(i))
    DO j=minreg,maxreg
      IF(reg_side_mesh(i,j)%num_el .GT. 0)CALL mesh_2d_skew_analysis(reg_side_mesh(i,j))
    ENDDO
  ENDDO
  !finish up progress bar
  DO i=prog,max_prog
    WRITE(*,'(A)',ADVANCE='NO')'*'
  ENDDO
  WRITE(*,*)

  WRITE(*,'(A)')'----------------------- Calculating 2D mesh aspect ratio:'
  CALL comp_tri_ar()
  !for each mesh structure
  DO i=1,6
    CALL mesh_2d_ar_analysis(tot_side_mesh(i))
    DO j=minreg,maxreg
      IF(reg_side_mesh(i,j)%num_el .GT. 0)CALL mesh_2d_ar_analysis(reg_side_mesh(i,j))
    ENDDO
  ENDDO
  !finish up progress bar
  DO i=prog,max_prog
    WRITE(*,'(A)',ADVANCE='NO')'*'
  ENDDO
  WRITE(*,*)

  WRITE(*,'(A)')'----------------------- Calculating 2D mesh aspect smoothness:'
  CALL comp_tri_smooth()
  !for each mesh structure
  DO i=1,6
    CALL mesh_2d_smooth_analysis(tot_side_mesh(i))
    DO j=minreg,maxreg
      IF(reg_side_mesh(i,j)%num_el .GT. 0)CALL mesh_2d_smooth_analysis(reg_side_mesh(i,j))
    ENDDO
  ENDDO
  !finish up progress bar
  DO i=prog,max_prog
    WRITE(*,'(A)',ADVANCE='NO')'*'
  ENDDO
  WRITE(*,*)

  WRITE(*,'(A)')'---------------------- Outputting results to '//TRIM(ADJUSTL(mesh_infile))//&
      '_stats.csv'
  CALL output_statistics()
  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'---------------------- OpenMeshQualityAnalyzer successful ----------------------'
ENDPROGRAM openmeshqualityanalyzer
