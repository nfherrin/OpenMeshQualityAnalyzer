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
    PUBLIC :: adjacency_calc

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
      tbound_cond(tot_bcf,2)=faceid
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

  !cross product
  FUNCTION cross(a, b)
    REAL(8) :: cross(3)
    REAL(8), INTENT(IN) :: a(3), b(3)

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  ENDFUNCTION cross

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
END MODULE boundary_conditions
