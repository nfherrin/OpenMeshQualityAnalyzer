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
    INTEGER :: i
    !temporary variables
    INTEGER :: temp_int

    OPEN(unit=20,FILE=mesh_infile,ACTION='READ',STATUS='OLD')

    !get number of tets and verts and get past placeholder values
    READ(20,*)num_verts
    READ(20,*)num_tets
    READ(20,*)
    READ(20,*)

    !read in the vertices
    ALLOCATE(vertex(num_verts,3))
    DO i=1,num_verts
      READ(20,*)temp_int,vertex(i,:)
      IF(i .NE. temp_int)STOP 'thrm vertices disordered'
    ENDDO

    !read in the element tags
    ALLOCATE(el_tag(num_tets))
    DO i=1,num_tets
      READ(20,*)temp_int,el_tag(i)
      IF(i .NE. temp_int)STOP 'thrm element tags disordered'
    ENDDO

    !read in the tet compisition (data/vertices)
    ALLOCATE(element(num_tets,4))
    DO i=1,num_tets
      READ(20,*)temp_int,element(i,:)
      IF(i .NE. temp_int)STOP 'thrm element compositions disordered'
    ENDDO

    !read in the boundary face data
    READ(20,*)num_bcf
    ALLOCATE(bc_data(num_bcf,2))
    DO i=1,num_bcf
      READ(20,*)bc_data(i,:)
    ENDDO

    !read in the adjancency data
    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    prog=0
    WRITE(*,'(A)',ADVANCE='NO')'**'
    READ(20,*)
    ALLOCATE(adj_list(num_tets*4,4))
    DO i=1,num_tets*4
      READ(20,*)adj_list(i,:)
      IF(MOD(i,CEILING(num_tets*4*1.0/(max_prog-1.0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO

    CLOSE(20)
    WRITE(*,'(A)')'*'
  ENDSUBROUTINE read_thrm_file
END MODULE read_thrm
