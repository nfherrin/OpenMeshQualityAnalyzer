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
  INTEGER :: num_verts

  !number of tets
  INTEGER :: num_tets

  !number of boundary condition faces
  INTEGER :: num_bcf

  !progress bar counting
  INTEGER :: prog

  !progress bar counting
  INTEGER, PARAMETER :: max_prog=70

  !vertex data (num_verts) locations of each vertex
  TYPE(vertex_type), TARGET, ALLOCATABLE :: vertex(:)

  !element data (num_tets,4) indices of vertices of element
  TYPE(element_type_3d), ALLOCATABLE :: tet(:)

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

  !region volumes and standard deviation (minreg:maxreg)
  REAL(8), ALLOCATABLE :: reg_vol(:),reg_vol_sd(:)

  !total volume and standard deviation
  REAL(8) :: tot_vol,tot_vol_sd

  !region average and standard deviation skew (minreg:maxreg)
  REAL(8), ALLOCATABLE :: reg_avg_skew(:),reg_sd_skew(:)

  !total average and standard deviation skew
  REAL(8) :: tot_avg_skew,tot_sd_skew

  !region average and standard deviation aspect ratio (minreg:maxreg)
  REAL(8), ALLOCATABLE :: reg_avg_ar(:),reg_sd_ar(:)

  !total average and standard deviation aspec ratio
  REAL(8) :: tot_avg_ar,tot_sd_ar

  !pi
  REAL(8),PARAMETER :: PI=4.D0*DATAN(1.D0)
CONTAINS

END MODULE globals
