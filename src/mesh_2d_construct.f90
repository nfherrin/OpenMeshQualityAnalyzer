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
  SUBROUTINE construct_bc_mesh(this_mesh,this_mesh_reg,side_id)
    TYPE(mesh_type_2d), INTENT(INOUT) :: this_mesh
    TYPE(mesh_type_2d), INTENT(INOUT) :: this_mesh_reg(minreg:maxreg)
    INTEGER, INTENT(IN) :: side_id
    INTEGER :: i=0,el_id=0
    INTEGER, ALLOCATABLE :: reg_tri_id(:)

    ALLOCATE(reg_tri_id(minreg:maxreg))
    reg_tri_id=0

    !only construct a side mesh if it's flat
    IF(side_flat(side_id))THEN
      DO i=1,tot_bcf
        IF(bc_data(i,3) .EQ. side_id)this_mesh%num_el=this_mesh%num_el+1
      ENDDO
      ALLOCATE(this_mesh%tri(this_mesh%num_el))
      el_id=0
      !assign total mesh elements and count elements for each region
      DO i=1,tot_bcf
        IF(bc_data(i,3) .EQ. side_id)THEN
          el_id=el_id+1
          this_mesh%tri(el_id)%p => tri(i)
          this_mesh_reg(tri(i)%reg)%num_el=this_mesh_reg(tri(i)%reg)%num_el+1
        ENDIF
      ENDDO
      !allocate and assign tris in the region
      DO i=minreg,maxreg
        ALLOCATE(this_mesh_reg(i)%tri(this_mesh_reg(i)%num_el))
      ENDDO
      DO i=1,tot_bcf
        IF(bc_data(i,3) .EQ. side_id)THEN
          reg_tri_id(tri(i)%reg)=reg_tri_id(tri(i)%reg)+1
          this_mesh_reg(tri(i)%reg)%tri(reg_tri_id(tri(i)%reg))%p => tri(i)
        ENDIF
      ENDDO
    ELSE
      !do nothing, the side is not flat and so won't be analyzed
    ENDIF
  ENDSUBROUTINE construct_bc_mesh

END MODULE mesh_2d_construct
