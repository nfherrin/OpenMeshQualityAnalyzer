!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to analyze 2d mesh quality indicators
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE mesh_analyze_2d
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: mesh_2d_area_analysis, mesh_2d_skew_analysis, mesh_2d_ar_analysis
  PUBLIC :: mesh_2d_smooth_analysis
CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !area analysis of a mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE mesh_2d_area_analysis(this_mesh)
    TYPE(mesh_type_2d), INTENT(INOUT) :: this_mesh
    INTEGER :: i

    !compute areas
    this_mesh%area=0
    DO i=1,this_mesh%num_el
      this_mesh%area=this_mesh%area+this_mesh%tri(i)%p%area
    ENDDO
    this_mesh%area_avg=this_mesh%area/(this_mesh%num_el*1.0D0)

    !compute area standar deviation
    this_mesh%area_sd=0
    DO i=1,this_mesh%num_el
      this_mesh%area_sd=this_mesh%area_sd+(this_mesh%tri(i)%p%area-this_mesh%area_avg)**2
    ENDDO
    this_mesh%area_sd=SQRT(this_mesh%area_sd/(this_mesh%num_el-1.0D0))
  ENDSUBROUTINE mesh_2d_area_analysis

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !skew analysis of a mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE mesh_2d_skew_analysis(this_mesh)
    TYPE(mesh_type_2d), INTENT(INOUT) :: this_mesh
    INTEGER :: i
    this_mesh%skew_avg=0
    DO i=1,this_mesh%num_el
      this_mesh%skew_avg=this_mesh%skew_avg+this_mesh%tri(i)%p%skew
    ENDDO
    this_mesh%skew_avg=this_mesh%skew_avg/(this_mesh%num_el*1.0D0)

    this_mesh%skew_sd=0
    DO i=1,this_mesh%num_el
      this_mesh%skew_sd=this_mesh%skew_sd+(this_mesh%tri(i)%p%skew-this_mesh%skew_avg)**2
    ENDDO
    this_mesh%skew_sd=SQRT(this_mesh%skew_sd/(this_mesh%num_el-1.0D0))
  ENDSUBROUTINE mesh_2d_skew_analysis

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !aspect ratio analysis of a mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE mesh_2d_ar_analysis(this_mesh)
    TYPE(mesh_type_2d), INTENT(INOUT) :: this_mesh
    INTEGER :: i
    this_mesh%ar_avg=0
    DO i=1,this_mesh%num_el
      this_mesh%ar_avg=this_mesh%ar_avg+this_mesh%tri(i)%p%aspect_ratio
    ENDDO
    this_mesh%ar_avg=this_mesh%ar_avg/(this_mesh%num_el*1.0D0)

    this_mesh%ar_sd=0
    DO i=1,this_mesh%num_el
      this_mesh%ar_sd=this_mesh%ar_sd+(this_mesh%tri(i)%p%aspect_ratio-this_mesh%ar_avg)**2
    ENDDO
    this_mesh%ar_sd=SQRT(this_mesh%ar_sd/(this_mesh%num_el-1.0D0))
  ENDSUBROUTINE mesh_2d_ar_analysis

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !smoothness analysis of a mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE mesh_2d_smooth_analysis(this_mesh)
    TYPE(mesh_type_2d), INTENT(INOUT) :: this_mesh
    INTEGER :: i
    this_mesh%smooth_avg=0
    DO i=1,this_mesh%num_el
      this_mesh%smooth_avg=this_mesh%smooth_avg+this_mesh%tri(i)%p%smoothness
    ENDDO
    this_mesh%smooth_avg=this_mesh%smooth_avg/(this_mesh%num_el*1.0D0)

    this_mesh%smooth_sd=0
    DO i=1,this_mesh%num_el
      this_mesh%smooth_sd=this_mesh%smooth_sd+(this_mesh%tri(i)%p%smoothness&
          -this_mesh%smooth_avg)**2
    ENDDO
    this_mesh%smooth_sd=SQRT(this_mesh%smooth_sd/(this_mesh%num_el-1.0D0))
  ENDSUBROUTINE mesh_2d_smooth_analysis
END MODULE mesh_analyze_2d
