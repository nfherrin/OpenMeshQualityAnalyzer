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

  !vertex data (num_verts,3) locations of each vertex
  REAL(8), ALLOCATABLE :: vertex(:,:)

  !element data (num_tets,4) indices of vertices of element
  INTEGER, ALLOCATABLE :: element(:,:)

  !element region tags (num_tets)
  INTEGER, ALLOCATABLE :: el_tag(:)

  !adjacency list (num_tets*4,4) first column is one element id, second column is that element's face, third column is other element id, fourth column is that other element's face
  !last two columns are zero for boundaries
  INTEGER, ALLOCATABLE :: adj_list(:,:)

  !boundary conditions data. (num_bcf,2) first data column is element id and second data column is face id
  INTEGER, ALLOCATABLE :: bc_data(:,:)

  !side flatnesses, ordered same as BCs.
  LOGICAL :: side_flat(6)=.FALSE.

  !minimum and maximum region bounds
  INTEGER :: minreg=0,maxreg=0

  !region number of tets (minreg:maxreg)
  INTEGER,ALLOCATABLE :: tets_in_reg(:)

  !region and total volumes (minreg:maxreg)
  REAL(8), ALLOCATABLE :: reg_vol(:)

  !total volume
  REAL(8) :: tot_vol

  !region maximum, minimum, average and standard deviation skew (minreg:maxreg)
  REAL(8), ALLOCATABLE :: reg_avg_skew(:),reg_sd_skew(:)

  !total maximum, minimum, average and standard deviation skew
  REAL(8) :: tot_max_skew,tot_min_skew,tot_avg_skew,tot_sd_skew

  !pi
  REAL(8),PARAMETER :: PI=4.D0*DATAN(1.D0)
CONTAINS

END MODULE globals
