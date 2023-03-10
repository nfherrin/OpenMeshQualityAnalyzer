!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!>    This module contains variables and functions used to store and manipulate
!!    data common to both the input mesh and output mesh.
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE globals
  USE mesh_types
  IMPLICIT NONE

  !input mesh filename
  CHARACTER(200) :: mesh_infile

  !number of vertices
  INTEGER :: tot_verts

  !number of tets
  INTEGER :: tot_tets

  !number of boundary condition faces
  INTEGER :: tot_bcf

  !progress bar counting
  INTEGER :: prog

  !progress bar counting
  INTEGER, PARAMETER :: max_prog=70

  !vertex data (tot_verts) locations of each vertex
  TYPE(vertex_type), TARGET, ALLOCATABLE :: vertex(:)

  !element data (tot_tets) indices of vertices of element
  TYPE(element_type_3d), TARGET, ALLOCATABLE :: tet(:)

  !triangles on the surfaces, even the non-flat ones
  TYPE(element_type_2d), TARGET, ALLOCATABLE :: tri(:)

  !boundary conditions data. (num_bcf,3) first data column is element id,
  !second data column is face id,
  !and third data column is direction 1 to 6 ordered -x +x -y +y -z +z
  INTEGER, ALLOCATABLE :: bc_data(:,:)

  !side flatness ordered -x +x -y +y -z +z
  LOGICAL :: side_flat(6)=.TRUE.

  !minimum and maximum region bounds
  INTEGER :: minreg=0,maxreg=0

  !total mesh
  TYPE(mesh_type_3d) :: tot_mesh

  !region based mesh
  TYPE(mesh_type_3d), ALLOCATABLE :: reg_mesh(:)

  !total side 2d mesh ordered as -x +x -y +y -z +z (so west, east, south, north, bot, top)
  TYPE(mesh_type_2d) :: tot_side_mesh(6)

  !region based side 2d mesh, first dimension is side second dimension is region
  TYPE(mesh_type_2d), ALLOCATABLE :: reg_side_mesh(:,:)

  !side bc locations (how far each side goes)
  REAL(8) :: bc_locs(6)=0

  !pi
  REAL(8),PARAMETER :: PI=4.D0*DATAN(1.D0)
CONTAINS

  !all the sorts are very small so bubble sort is fine
  SUBROUTINE bubble_sort(vec)
    REAL(8),INTENT(INOUT) :: vec(:)
    INTEGER :: m,i,changes
    REAL(8) :: tempr

    m=SIZE(vec)
    DO
      changes=0
      DO i=1,m-1
        IF(vec(i) .GT. vec(i+1))THEN
          tempr=vec(i)
          vec(i)=vec(i+1)
          vec(i+1)=tempr
          changes=changes+1
        ENDIF
      ENDDO
      IF(changes .EQ. 0)EXIT
    ENDDO
  ENDSUBROUTINE bubble_sort

END MODULE globals
