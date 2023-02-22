!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!>    This module contains variables and functions used to store and manipulate
!!    data common to both the input mesh and output mesh.
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE mesh_types
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: vertex_type, element_type_2d, element_type_3d, mesh_type_2d, mesh_type_3d

  !vertex pointer necessary to make array of pointers
  TYPE :: vert_ptr
    !actual pointer to the vertex type
    TYPE(vertex_type), POINTER :: p
  ENDTYPE

  !the vertex type
  TYPE :: vertex_type
    !location for x,y,z
    REAL(8) :: x=0
    REAL(8) :: y=0
    REAL(8) :: z=0
    !vertex id/index
    INTEGER :: id=0
    CONTAINS
      !compute distance to another point
      PROCEDURE :: distance
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
    TYPE(vert_ptr) :: corner(3)
    !area of the triangle
    REAL(8) :: area=0
  ENDTYPE

  !the specific 3d element type
  TYPE, EXTENDS(base_element_type) :: element_type_3d
    !4 corners of the tet
    TYPE(vert_ptr) :: corner(4)
    !volume of the tet
    REAL(8) :: vol=0
    !Adjacent tet id for faces 1 to 4
    INTEGER :: adj_id(4)=0
    !Adjacent tet face for faces 1 to 4
    INTEGER :: adj_face(4)=0
    CONTAINS
      !compute the circumsphere radius
      PROCEDURE :: sphere_rad
  ENDTYPE

  !the base mesh type
  TYPE, ABSTRACT :: base_mesh_type
    !mesh id (0 for total mesh)
    INTEGER :: id=0
    !number of elements
    INTEGER :: num_el=0
    !element average skew
    REAL(8) :: skew_avg=0
    !element skew standard deviation
    REAL(8) :: skew_sd=0
    !element average aspect ratio
    REAL(8) :: ar_avg=0
    !element aspect ratio standard deviation
    REAL(8) :: ar_sd=0
  ENDTYPE

  !tri pointer necessary to make array of pointers
  TYPE :: tri_ptr
    !actual pointer to the tri type
    TYPE(element_type_2d), POINTER :: p
  ENDTYPE

  !the 2d mesh type
  TYPE, EXTENDS(base_mesh_type) :: mesh_type_2d
    !triangle elements in the mesh
    TYPE(tri_ptr), ALLOCATABLE :: tri(:)
    !total area
    REAL(8) :: area=0
    !average tri area
    REAL(8) :: area_avg=0
    !tri area standard deviation
    REAL(8) :: area_sd=0
  ENDTYPE

  !tet pointer necessary to make array of pointers
  TYPE :: tet_ptr
    !actual pointer to the tet type
    TYPE(element_type_3d), POINTER :: p
  ENDTYPE

  !the 3d mesh type
  TYPE, EXTENDS(base_mesh_type) :: mesh_type_3d
    !tet elements in the mesh
    TYPE(tet_ptr), ALLOCATABLE :: tet(:)
    !total volume
    REAL(8) :: vol=0
    !average tet volume
    REAL(8) :: vol_avg=0
    !tet volume standard deviation
    REAL(8) :: vol_sd=0
  ENDTYPE
CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !computes radius of the circumsphere around a tet
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(8) FUNCTION sphere_rad(this_tet)
    CLASS(element_type_3d), INTENT(IN) :: this_tet

    REAL(8) :: ba(3),ca(3),da(3)
    REAL(8) :: len_ba,len_ca,len_da
    REAL(8) :: cross_cd(3),cross_db(3),cross_bc(3)
    REAL(8) :: denominator
    REAL(8) :: circ(3)

    ba(1)=this_tet%corner(2)%p%x-this_tet%corner(1)%p%x
    ba(2)=this_tet%corner(2)%p%y-this_tet%corner(1)%p%y
    ba(3)=this_tet%corner(2)%p%z-this_tet%corner(1)%p%z

    ca(1)=this_tet%corner(3)%p%x-this_tet%corner(1)%p%x
    ca(2)=this_tet%corner(3)%p%y-this_tet%corner(1)%p%y
    ca(3)=this_tet%corner(3)%p%z-this_tet%corner(1)%p%z

    da(1)=this_tet%corner(4)%p%x-this_tet%corner(1)%p%x
    da(2)=this_tet%corner(4)%p%y-this_tet%corner(1)%p%y
    da(3)=this_tet%corner(4)%p%z-this_tet%corner(1)%p%z

    len_ba=ba(1)*ba(1)+ba(2)*ba(2)+ba(3)*ba(3);
    len_ca=ca(1)*ca(1)+ca(2)*ca(2)+ca(3)*ca(3);
    len_da=da(1)*da(1)+da(2)*da(2)+da(3)*da(3);

    cross_cd(1)=ca(2)*da(3)-da(2)*ca(3);
    cross_cd(2)=ca(3)*da(1)-da(3)*ca(1);
    cross_cd(3)=ca(1)*da(2)-da(1)*ca(2);

    cross_db(1)=da(2)*ba(3)-ba(2)*da(3);
    cross_db(2)=da(3)*ba(1)-ba(3)*da(1);
    cross_db(3)=da(1)*ba(2)-ba(1)*da(2);

    cross_bc(1)=ba(2)*ca(3)-ca(2)*ba(3);
    cross_bc(2)=ba(3)*ca(1)-ca(3)*ba(1);
    cross_bc(3)=ba(1)*ca(2)-ca(1)*ba(2);

    denominator=0.5D0/(ba(1)*cross_cd(1)+ba(2)*cross_cd(2)+ba(3)*cross_cd(3));

    circ(1)=(len_ba*cross_cd(1)+len_ca*cross_db(1)+len_da*cross_bc(1))*denominator
    circ(2)=(len_ba*cross_cd(2)+len_ca*cross_db(2)+len_da*cross_bc(2))*denominator
    circ(3)=(len_ba*cross_cd(3)+len_ca*cross_db(3)+len_da*cross_bc(3))*denominator

    sphere_rad=SQRT(circ(1)**2+circ(2)**2+circ(3)**2)
  ENDFUNCTION sphere_rad

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !computes the distance to another vertex
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(8) FUNCTION distance(a,b)
      CLASS(vertex_type), INTENT(IN) :: a
      TYPE(vertex_type), INTENT(IN) :: b

      distance=SQRT((a%x-b%x)**2+(a%y-b%y)**2+(a%z-b%z)**2)
    ENDFUNCTION distance

END MODULE mesh_types
