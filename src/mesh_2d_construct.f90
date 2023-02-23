!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to construct a 2d mesh from a side
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE mesh_2d_construct
  USE globals
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: construct_bc_mesh
CONTAINS

  !reads in a gmsh file
  SUBROUTINE construct_bc_mesh(this_mesh,side_id)
    TYPE(mesh_type_2d), INTENT(INOUT) :: this_mesh
    INTEGER, INTENT(IN) :: side_id

    IF(is_flat(side_id))THEN
      !WRITE(*,*)'it is flat'
    ELSE
      !WRITE(*,*)'it is not flat'
    ENDIF
    this_mesh%num_el=side_id
  ENDSUBROUTINE construct_bc_mesh

  LOGICAL FUNCTION is_flat(side_id)
    INTEGER, INTENT(IN) :: side_id
    is_flat=.TRUE.
    IF(MOD(side_id,2) == 0)is_flat=.TRUE.
  ENDFUNCTION is_flat

END MODULE mesh_2d_construct
