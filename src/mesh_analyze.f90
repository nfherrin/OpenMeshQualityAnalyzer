!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to analyze mesh quality indicators
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE mesh_analyze
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: assign_tets,mesh_vol_analysis
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !assigns tets to each region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE assign_tets()
    INTEGER :: i
    INTEGER, ALLOCATABLE :: reg_tet_id(:)

    !allocate the region mesh structures
    minreg=MINVAL(tet(:)%reg)
    maxreg=MAXVAL(tet(:)%reg)
    ALLOCATE(reg_mesh(minreg:maxreg),reg_tet_id(minreg:maxreg))
    reg_tet_id=0
    tot_mesh%num_el=tot_tets
    ALLOCATE(tot_mesh%tet(tot_mesh%num_el))

    !count number of elements in each region
    DO i=1,tot_tets
      reg_mesh(tet(i)%reg)%num_el=reg_mesh(tet(i)%reg)%num_el+1
      tot_mesh%tet(i)%p => tet(i)
    ENDDO
    !allocate tets in the region
    DO i=minreg,maxreg
      ALLOCATE(reg_mesh(i)%tet(reg_mesh(i)%num_el))
    ENDDO
    DO i=1,tot_tets
      reg_tet_id(tet(i)%reg)=reg_tet_id(tet(i)%reg)+1
      reg_mesh(tet(i)%reg)%tet(reg_tet_id(tet(i)%reg))%p => tet(i)
    ENDDO
  ENDSUBROUTINE assign_tets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate volumes of each tet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE mesh_vol_analysis(this_mesh)
    INTEGER :: i
    TYPE(mesh_type_3d), INTENT(INOUT) :: this_mesh
    this_mesh%vol=0
    DO i=1,this_mesh%num_el
      this_mesh%vol=this_mesh%vol+this_mesh%tet(i)%p%vol
    ENDDO
    this_mesh%vol_avg=this_mesh%vol/(this_mesh%num_el*1.0D0)

    this_mesh%vol_sd=0
    DO i=1,this_mesh%num_el
      this_mesh%vol_sd=this_mesh%vol_sd+(this_mesh%tet(i)%p%vol-this_mesh%vol_avg)**2
    ENDDO
    this_mesh%vol_sd=SQRT(this_mesh%vol_sd/(this_mesh%num_el-1.0D0))
  ENDSUBROUTINE mesh_vol_analysis
END MODULE mesh_analyze
