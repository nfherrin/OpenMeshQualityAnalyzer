!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to analyze mesh quality indicators
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE mesh_analyze
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: assign_tets,mesh_vol_analysis,mesh_skew_analysis,mesh_ar_analysis,mesh_smooth_analysis
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !assigns tets to each region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE assign_tets()
    INTEGER, ALLOCATABLE :: reg_tet_id(:)
    INTEGER :: i

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
  !volume analysis of a mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE mesh_vol_analysis(this_mesh)
    TYPE(mesh_type_3d), INTENT(INOUT) :: this_mesh
    INTEGER :: i
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !skew analysis of a mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE mesh_skew_analysis(this_mesh)
      TYPE(mesh_type_3d), INTENT(INOUT) :: this_mesh
      INTEGER :: i
      this_mesh%skew_avg=0
      DO i=1,this_mesh%num_el
        this_mesh%skew_avg=this_mesh%skew_avg+this_mesh%tet(i)%p%skew
      ENDDO
      this_mesh%skew_avg=this_mesh%skew_avg/(this_mesh%num_el*1.0D0)

      this_mesh%skew_sd=0
      DO i=1,this_mesh%num_el
        this_mesh%skew_sd=this_mesh%skew_sd+(this_mesh%tet(i)%p%skew-this_mesh%skew_avg)**2
      ENDDO
      this_mesh%skew_sd=SQRT(this_mesh%skew_sd/(this_mesh%num_el-1.0D0))
    ENDSUBROUTINE mesh_skew_analysis

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !aspect ratio analysis of a mesh
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE mesh_ar_analysis(this_mesh)
        TYPE(mesh_type_3d), INTENT(INOUT) :: this_mesh
        INTEGER :: i
        this_mesh%ar_avg=0
        DO i=1,this_mesh%num_el
          this_mesh%ar_avg=this_mesh%ar_avg+this_mesh%tet(i)%p%aspect_ratio
        ENDDO
        this_mesh%ar_avg=this_mesh%ar_avg/(this_mesh%num_el*1.0D0)

        this_mesh%ar_sd=0
        DO i=1,this_mesh%num_el
          this_mesh%ar_sd=this_mesh%ar_sd+(this_mesh%tet(i)%p%aspect_ratio-this_mesh%ar_avg)**2
        ENDDO
        this_mesh%ar_sd=SQRT(this_mesh%ar_sd/(this_mesh%num_el-1.0D0))
      ENDSUBROUTINE mesh_ar_analysis

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !smoothness analysis of a mesh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE mesh_smooth_analysis(this_mesh)
          TYPE(mesh_type_3d), INTENT(INOUT) :: this_mesh
          INTEGER :: i
          this_mesh%smooth_avg=0
          DO i=1,this_mesh%num_el
            this_mesh%smooth_avg=this_mesh%smooth_avg+this_mesh%tet(i)%p%smoothness
          ENDDO
          this_mesh%smooth_avg=this_mesh%smooth_avg/(this_mesh%num_el*1.0D0)

          this_mesh%smooth_sd=0
          DO i=1,this_mesh%num_el
            this_mesh%smooth_sd=this_mesh%smooth_sd+(this_mesh%tet(i)%p%smoothness&
                -this_mesh%smooth_avg)**2
          ENDDO
          this_mesh%smooth_sd=SQRT(this_mesh%smooth_sd/(this_mesh%num_el-1.0D0))
        ENDSUBROUTINE mesh_smooth_analysis
END MODULE mesh_analyze
