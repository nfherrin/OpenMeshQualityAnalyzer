!OpenMeshQualityAnalyzer is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to ingest a file in the
!! THOR mesh format (.thrm)
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE read_thrm
  USE globals
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_thrm_file
CONTAINS

  !reads in a gmsh file
  SUBROUTINE read_thrm_file()
    !counting variables
    INTEGER :: i,j
    !temporary variables
    INTEGER :: temp_int,el_vec(4),ti2

    OPEN(unit=20,FILE=mesh_infile,ACTION='READ',STATUS='OLD')

    !get number of tets and verts and get past placeholder values
    READ(20,*)tot_verts
    READ(20,*)tot_tets
    READ(20,*)
    READ(20,*)

    !read in the vertices
    ALLOCATE(vertex(tot_verts))
    DO i=1,tot_verts
      READ(20,*)temp_int,vertex(i)%x,vertex(i)%y,vertex(i)%z
      vertex(i)%id=temp_int
      IF(i .NE. temp_int)STOP 'thrm vertices disordered'
    ENDDO

    !read in the element tags
    ALLOCATE(tet(tot_tets))
    DO i=1,tot_tets
      READ(20,*)temp_int,tet(i)%reg
      tet(i)%id=temp_int
      IF(i .NE. temp_int)STOP 'thrm element tags disordered'
    ENDDO

    !read in the tet composition (data/vertices)
    DO i=1,tot_tets
      READ(20,*)temp_int,el_vec(:)
      DO j=1,4
        tet(i)%corner(j)%p => vertex(el_vec(j))
      ENDDO
      IF(i .NE. temp_int)STOP 'thrm element compositions disordered'
    ENDDO

    !read in the boundary face data
    READ(20,*)tot_bcf
    ALLOCATE(bc_data(tot_bcf,3))
    DO i=1,tot_bcf
      READ(20,*)bc_data(i,1:2)
    ENDDO

    !read in the adjancency data
    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=1
    READ(20,*)
    DO i=1,tot_tets
      !loop over for a single tet
      DO j=1,4
        READ(20,*)temp_int,ti2,tet(i)%adj_id(j),tet(i)%adj_face(j)
        IF(i .NE. temp_int)STOP 'thrm adjacencies disordered'
        IF(j-1 .NE. ti2)STOP 'thrm adjacencies disordered'
      ENDDO
      IF(MOD(i,CEILING(tot_tets*4*1.0/(max_prog-1.0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
    !adjust face id to be indexed 1 to 4
    !also check to make sure that all boundary condition faces are defined as zero
    DO i=1,tot_tets
      DO j=1,4
        IF(tet(i)%adj_id(j) .NE. 0)THEN
          tet(i)%adj_face(j)=tet(i)%adj_face(j)+1
        ELSE
          IF(tet(i)%adj_face(j) .NE. 0)STOP 'boundary conditions should have faces defined as zero'
        ENDIF
      ENDDO
    ENDDO

    CLOSE(20)
    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,'(A)')'*'
  ENDSUBROUTINE read_thrm_file
END MODULE read_thrm
