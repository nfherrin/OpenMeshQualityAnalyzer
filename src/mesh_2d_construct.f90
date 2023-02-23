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

  !constructs a 2D bc mesh for a given side
  SUBROUTINE construct_bc_mesh(this_mesh,side_id)
    TYPE(mesh_type_2d), INTENT(INOUT) :: this_mesh
    INTEGER, INTENT(IN) :: side_id

    IF(side_flat(side_id))THEN
      WRITE(*,*)'it is flat'
    ELSE
      !do nothing, the side is not flat and so won't be analyzed
    ENDIF
    write(*,*)this_mesh%num_el,side_id
  ENDSUBROUTINE construct_bc_mesh

END MODULE mesh_2d_construct
