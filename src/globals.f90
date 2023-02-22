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

  !boundary conditions data. (num_bcf,2) first data column is element id and second data column is face id
  INTEGER, ALLOCATABLE :: bc_data(:,:)

  !side flatnesses, ordered same as BCs.
  LOGICAL :: side_flat(6)=.FALSE.

  !minimum and maximum region bounds
  INTEGER :: minreg=0,maxreg=0

  !total mesh
  TYPE(mesh_type_3d) :: tot_mesh

  !region based mesh
  TYPE(mesh_type_3d), ALLOCATABLE :: reg_mesh(:)

  !pi
  REAL(8),PARAMETER :: PI=4.D0*DATAN(1.D0)
CONTAINS

END MODULE globals
