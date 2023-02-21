!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!>    This module contains variables and functions used to store and manipulate
!!    data common to both the input mesh and output mesh.
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE globals
  IMPLICIT NONE

  !input mesh filename
  CHARACTER(200) :: mesh_infile

  !number of vertices
  INTEGER :: num_verts

  !number of tets
  INTEGER :: num_tets

  !number of boundary condition faces
  INTEGER :: num_bcf

  !progress bar counting
  INTEGER :: prog

  !progress bar counting
  INTEGER, PARAMETER :: max_prog=70

  !vertex data
  REAL(8), ALLOCATABLE :: vertex(:,:)

  !element data
  INTEGER, ALLOCATABLE :: element(:,:)

  !element region tags
  INTEGER, ALLOCATABLE :: el_tag(:)

  !adjacency list
  INTEGER, ALLOCATABLE :: adj_list(:,:)

  !boundary conditions data. num_bcf,data first data column is element id and second data column is face id
  INTEGER, ALLOCATABLE :: bc_data(:,:)

  !side flatnesses, ordered same as BCs.
  LOGICAL :: side_flat(6)=.FALSE.

  !pi
  REAL(8),PARAMETER :: PI=4.D0*DATAN(1.D0)
CONTAINS

END MODULE globals
