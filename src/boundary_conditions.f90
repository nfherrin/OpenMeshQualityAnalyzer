!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to find boundary
!! conditions and adjacency info for a given set of elements
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE boundary_conditions
    USE globals
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: adjacency_calc, compute_bc_sides,generate_bc_triangles

    INTEGER, ALLOCATABLE :: tbound_cond(:,:)
CONTAINS

  !calculate adjacencies between tets
  SUBROUTINE adjacency_calc()
    INTEGER :: i,og_face(3),adj_idx

    DO i=1,tot_tets
      CALL orderverts(tet(i))
    ENDDO

    ALLOCATE(tbound_cond(tot_tets*4,2))
    tbound_cond=0
    !loop over all tets
    adj_idx=0
    tot_bcf=0
    prog=0
    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    DO i=1,tot_tets
      IF(MOD(i,CEILING(tot_tets*1.0/(max_prog-1.0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
      !first face
      og_face=(/tet(i)%corner(2)%p%id,tet(i)%corner(3)%p%id,tet(i)%corner(4)%p%id/)
      CALL find_adj(og_face,i,0,adj_idx)
      !second face
      og_face=(/tet(i)%corner(1)%p%id,tet(i)%corner(3)%p%id,tet(i)%corner(4)%p%id/)
      CALL find_adj(og_face,i,1,adj_idx)
      !third face
      og_face=(/tet(i)%corner(1)%p%id,tet(i)%corner(2)%p%id,tet(i)%corner(4)%p%id/)
      CALL find_adj(og_face,i,2,adj_idx)
      !fourth face
      og_face=(/tet(i)%corner(1)%p%id,tet(i)%corner(2)%p%id,tet(i)%corner(3)%p%id/)
      CALL find_adj(og_face,i,3,adj_idx)
    ENDDO
    ALLOCATE(bc_data(tot_bcf,3))
    bc_data=0
    DO i=1,tot_bcf
      bc_data(i,1:2)=tbound_cond(i,:)
    ENDDO

    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)
    DEALLOCATE(tbound_cond)
  ENDSUBROUTINE adjacency_calc

  !find adjacency for a given face
  SUBROUTINE find_adj(face,el_idx,faceid,adj_idx)
    INTEGER,INTENT(IN) :: face(3)
    INTEGER,INTENT(IN) :: el_idx
    INTEGER,INTENT(IN) :: faceid
    INTEGER,INTENT(INOUT) :: adj_idx
    INTEGER :: j,comp_face(3)
    LOGICAL :: match

    match=.FALSE.
    DO j=1,tot_tets
      !compare for first face
      comp_face=(/tet(j)%corner(2)%p%id,tet(j)%corner(3)%p%id,tet(j)%corner(4)%p%id/)
      CALL check_face(face,comp_face,el_idx,j,faceid,0,adj_idx,match)
      IF(match)EXIT
      !compare for second face
      comp_face=(/tet(j)%corner(1)%p%id,tet(j)%corner(3)%p%id,tet(j)%corner(4)%p%id/)
      CALL check_face(face,comp_face,el_idx,j,faceid,1,adj_idx,match)
      IF(match)EXIT
      !compare for third face
      comp_face=(/tet(j)%corner(1)%p%id,tet(j)%corner(2)%p%id,tet(j)%corner(4)%p%id/)
      CALL check_face(face,comp_face,el_idx,j,faceid,2,adj_idx,match)
      IF(match)EXIT
      !compare for fourth face
      comp_face=(/tet(j)%corner(1)%p%id,tet(j)%corner(2)%p%id,tet(j)%corner(3)%p%id/)
      CALL check_face(face,comp_face,el_idx,j,faceid,3,adj_idx,match)
      IF(match)EXIT
    ENDDO
    !if we go through the whole thing and don't exit
    IF(j .EQ. tot_tets+1)THEN
      !we didn't find a matching face so it's a boundary condition
      adj_idx=adj_idx+1
      tet(el_idx)%adj_id(faceid+1)=0
      tet(el_idx)%adj_face(faceid+1)=0
      tot_bcf=tot_bcf+1
      tbound_cond(tot_bcf,1)=el_idx
      tbound_cond(tot_bcf,2)=faceid+1
    ENDIF
  ENDSUBROUTINE find_adj

  !check to see if two faces match
  SUBROUTINE check_face(face1,face2,el_idx1,el_idx2,faceid1,faceid2,adj_idx,match)
    INTEGER,INTENT(IN) :: face1(3)
    INTEGER,INTENT(IN) :: face2(3)
    INTEGER,INTENT(IN) :: el_idx1
    INTEGER,INTENT(IN) :: el_idx2
    INTEGER,INTENT(IN) :: faceid1
    INTEGER,INTENT(IN) :: faceid2
    INTEGER,INTENT(INOUT) :: adj_idx
    LOGICAL,INTENT(INOUT) :: match
    match=.FALSE.
    !check to see if the faces match
    IF(face1(1) .EQ. face2(1) .AND. face1(2) .EQ. face2(2) &
        .AND. face1(3) .EQ. face2(3) .AND. el_idx1 .NE. el_idx2)THEN
      adj_idx=adj_idx+1
      tet(el_idx1)%adj_id(faceid1+1)=el_idx2
      tet(el_idx1)%adj_face(faceid1+1)=faceid2+1
      match=.TRUE.
    ENDIF
  ENDSUBROUTINE check_face

  !order the vertices by bubble sort for a given tet
  SUBROUTINE orderverts(this_tet)
    TYPE(element_type_3d), INTENT(INOUT) :: this_tet
    TYPE(vertex_type), POINTER :: temp_vert
    INTEGER :: i,changes

    !bubble sort algorithm, pretty cheap for only 4 points
    DO
      changes=0
      DO i=1,3
        IF(this_tet%corner(i)%p%id .GE. this_tet%corner(i+1)%p%id)THEN
          temp_vert => vertex(this_tet%corner(i)%p%id)
          this_tet%corner(i)%p => vertex(this_tet%corner(i+1)%p%id)
          this_tet%corner(i+1)%p => vertex(temp_vert%id)
          changes=changes+1
        ENDIF
      ENDDO
      IF(changes .EQ. 0)EXIT
    ENDDO
  ENDSUBROUTINE orderverts

  SUBROUTINE compute_bc_sides()
    INTEGER :: i,j,face_idx
    REAL(8) :: face_point(3,3),ext_point(3),norm_vec(3),lambda,offset
    bc_locs(1)=MINVAL(vertex(:)%x)
    bc_locs(2)=MAXVAL(vertex(:)%x)
    bc_locs(3)=MINVAL(vertex(:)%y)
    bc_locs(4)=MAXVAL(vertex(:)%y)
    bc_locs(5)=MINVAL(vertex(:)%z)
    bc_locs(6)=MAXVAL(vertex(:)%z)
    ext_point=0

    IF(MINVAL(bc_data(:,2)) .LE. 0 .OR. MAXVAL(bc_data(:,2)) .GE. 5)THEN
      WRITE(*,*)'error, bc faces werent increased',MINVAL(bc_data(:,2)),MAXVAL(bc_data(:,2))
      STOP 'here1'
    ENDIF
    !assume that sides are flat, figure out if they're not...
    side_flat=.TRUE.
    !go through all bc faces
    DO i=1,tot_bcf
      !assign the face and extruded points for the cross product
      face_idx=0
      DO j=1,4
        IF(bc_data(i,2) .EQ. j)THEN
          !if it equals the face index the it's the extruded point
          ext_point(:)=(/tet(bc_data(i,1))%corner(j)%p%x,tet(bc_data(i,1))%corner(j)%p%y,&
              tet(bc_data(i,1))%corner(j)%p%z/)
        ELSE
          !if doesn't equal the face index the it's one of the face points
          face_idx=face_idx+1
          face_point(face_idx,:)=(/tet(bc_data(i,1))%corner(j)%p%x, &
              tet(bc_data(i,1))%corner(j)%p%y,tet(bc_data(i,1))%corner(j)%p%z/)
        ENDIF
      ENDDO
      !get the outward going normal vector for the tet for this face
      norm_vec=cross(face_point(2,:)-face_point(1,:), face_point(3,:)-face_point(1,:))
      offset=face_point(1,1)*norm_vec(1)+face_point(1,2)*norm_vec(2)+face_point(1,3)*norm_vec(3)
      lambda=(offset-norm_vec(1)*ext_point(1)-norm_vec(2)*ext_point(2)-norm_vec(3)*ext_point(3)) &
          /(norm_vec(1)**2+norm_vec(2)**2+norm_vec(3)**2)
      norm_vec=norm_vec*lambda
      norm_vec=norm_vec/(SQRT(norm_vec(1)**2+norm_vec(2)**2+norm_vec(3)**2))

      !figure out which side the bc is on
      IF(ABS(MAXVAL(norm_vec)) .GT. ABS(MINVAL(norm_vec)))THEN
        !the primary direction is positive, so this is a +x, +y, or +z face
        IF(MAXLOC(norm_vec,1) .EQ. 1)THEN
          !+x
          bc_data(i,3)=2
        ELSEIF(MAXLOC(norm_vec,1) .EQ. 2)THEN
          !+y
          bc_data(i,3)=4
        ELSE
          !+z
          bc_data(i,3)=6
        ENDIF
      ELSE
        !the primary direction is negative, so this is a -x, -y, or -z face
        IF(MINLOC(norm_vec,1) .EQ. 1)THEN
          !-x
          bc_data(i,3)=1
        ELSEIF(MINLOC(norm_vec,1) .EQ. 2)THEN
          !-y
          bc_data(i,3)=3
        ELSE
          !-z
          bc_data(i,3)=5
        ENDIF
      ENDIF

      !check fot see if the side is flat, it only takes one to make it not flat
      IF(MAXVAL(ABS(norm_vec))-1.0D0 .LE. -1.0D-14)side_flat(bc_data(i,3))=.FALSE.
    ENDDO
  ENDSUBROUTINE compute_bc_sides

  SUBROUTINE generate_bc_triangles()
    INTEGER :: i,j,ii,jj,point_i
    ALLOCATE(tri(tot_bcf))

    !assign vertex pointers
    DO i=1,tot_bcf
      point_i=0
      tri(i)%id=i
      tri(i)%reg=tet(bc_data(i,1))%reg
      DO j=1,4
        !point to vertices if they are not the extruded (internal) point
        IF(j .NE. bc_data(i,2))THEN
          point_i=point_i+1
          tri(i)%corner(point_i)%p => tet(bc_data(i,1))%corner(j)%p
        ENDIF
      ENDDO
    ENDDO

    prog=0
    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    !compute the tri adjacencies
    DO i=1,tot_bcf
      IF(MOD(i,CEILING(tot_bcf*1.0/(max_prog-1.0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
      !loop over the sides of the tri
      DO j=1,3
        !only check if the adjacency hasn't already been found to avoid redundancy
        IF(tri(i)%adj_id(j) .EQ. 0)THEN
          !loop over all other tris
          DO ii=1,tot_bcf
            !loop over other tri's sides to see if they match only if it's on the same side to
            !check if any match
            IF(bc_data(i,3) .EQ. bc_data(ii,3))THEN
              DO jj=1,3
                !if the sides match, then assign the adjancencies
                IF(check_side(tri(i),j,tri(ii),jj))THEN
                  tri(i)%adj_id(j)=ii
                  tri(ii)%adj_id(jj)=i
                  tri(i)%adj_side(j)=jj
                  tri(ii)%adj_side(jj)=j
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)
  ENDSUBROUTINE generate_bc_triangles

  !check to see if two sides match
  LOGICAL FUNCTION check_side(tri1,side1,tri2,side2)
    TYPE(element_type_2d), INTENT(IN) :: tri1,tri2
    INTEGER, INTENT(IN) :: side1,side2
    INTEGER :: p1,p2,pp1,pp2
    check_side=.FALSE.
    !find first side points
    SELECTCASE(side1)
      CASE(1)
        p1=tri1%corner(2)%p%id
        p2=tri1%corner(3)%p%id
      CASE(2)
        p1=tri1%corner(1)%p%id
        p2=tri1%corner(3)%p%id
      CASE(3)
        p1=tri1%corner(1)%p%id
        p2=tri1%corner(2)%p%id
      CASE DEFAULT
        STOP 'side 1 must be 1 to 3'
    ENDSELECT
    !find second side points
    SELECTCASE(side2)
      CASE(1)
        pp1=tri2%corner(2)%p%id
        pp2=tri2%corner(3)%p%id
      CASE(2)
        pp1=tri2%corner(1)%p%id
        pp2=tri2%corner(3)%p%id
      CASE(3)
        pp1=tri2%corner(1)%p%id
        pp2=tri2%corner(2)%p%id
      CASE DEFAULT
        STOP 'side 2 must be 1 to 3'
    ENDSELECT
    IF(pp1 .EQ. p1 .AND. pp2 .EQ. p2)THEN
      check_side=.TRUE.
    ELSEIF(pp2 .EQ. p1 .AND. pp1 .EQ. p2)THEN
      check_side=.TRUE.
    ELSE
      check_side=.FALSE.
    ENDIF
  ENDFUNCTION

  !cross product
  FUNCTION cross(a, b)
    REAL(8) :: cross(3)
    REAL(8), INTENT(IN) :: a(3), b(3)

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  ENDFUNCTION cross
END MODULE boundary_conditions
