!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!>    This module contains variables and functions used to store and manipulate
!!    data common to both the input mesh and output mesh.
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE mesh_types
  IMPLICIT NONE

  !the vertex type
  TYPE :: vertex_type
    !location for x,y,z
    REAL(8) :: x=0
    REAL(8) :: y=0
    REAL(8) :: z=0
    !vertex id/index
    INTEGER :: id=0
  ENDTYPE

  !the base element type
  TYPE, ABSTRACT :: base_element_type
    !element region
    INTEGER :: reg=0
    !element id/index
    INTEGER :: id=0
    !skewness of the element
    REAL(8) :: skew=0
    !aspect ratio of the element
    REAL(8) :: aspect_ratio=0
  ENDTYPE

  !the specific 2d element type
  TYPE, EXTENDS(base_element_type) :: element_type_2d
    !3 corners of the triangle
    TYPE(vertex_type), POINTER :: corner(:)
    !area of the triangle
    REAL(8) :: area=0
  ENDTYPE

  !the specific 3d element type
  TYPE, EXTENDS(base_element_type) :: element_type_3d
    !4 corners of the tet
    TYPE(vertex_type), POINTER :: corner(:)
    !volume of the tet
    REAL(8) :: vol=0
  ENDTYPE

  !the base mesh type
  TYPE, ABSTRACT :: base_mesh_type
    !mesh id (0 for total mesh)
    INTEGER :: id=0
  ENDTYPE

  !the 2d mesh type
  TYPE, EXTENDS(base_mesh_type) :: mesh_type_2d
    !triangle elements in the mesh
    TYPE(element_type_2d), POINTER :: tri(:)
  ENDTYPE

  !the 3d mesh type
  TYPE, EXTENDS(base_mesh_type) :: mesh_type_3d
    !tet elements in the mesh
    TYPE(element_type_2d), POINTER :: tet(:)
  ENDTYPE
CONTAINS

END MODULE mesh_types
