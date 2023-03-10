#-----------------------------------------------------------------------------------
#Compile Options
#-----------------------------------------------------------------------------------
FC := gfortran
EXE := OpenMeshQualityAnalyzer
EXT := .exe

all:			TYPE :=
debug:    TYPE := _Debug
catch:    TYPE := _Catch
coverage: TYPE := _Coverage

all:      OPT := -O3 -Wall -Werror
debug:    OPT := -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow -fcheck=all -Wall -Werror
catch:    OPT := -O0 -g -Wall -W -Wsurprising -Wextra -Werror
coverage: OPT := -O0 -g --coverage -Wall -Werror

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not compiled using the generic rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------
SRC := \
				main.f90 \
				globals.f90 \
				mesh_analyze_3d.f90 \
				mesh_analyze_2d.f90 \
				tet_analyze.f90 \
				tri_analyze.f90 \
				boundary_conditions.f90 \
				out_stats.f90 \
				mesh_types.f90 \
				mesh_2d_construct.f90 \
				read_thrm.f90 \
				read_gmsh.f90

OBJ := $(SRC:.f90=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \

OBJ_FILTER := $(FILTER:.f90=.o)
MOD_FILTER := $(FILTER:.f90=.mod)

#-----------------------------------------------------------------------------------
#Complex external dependencies. Each of these requires a custom build rule.
#Intended for dependencies with their own makefiles or other build system
#-----------------------------------------------------------------------------------
COMP_DEP := \
#Phony targets for cleaning and building

#-----------------------------------------------------------------------------------
#Phony targets for cleaning and building
#-----------------------------------------------------------------------------------
.PHONY: all clean reset

all: $(EXE)

MPparallel: $(EXE)

debug: $(EXE)

catch: $(EXE)

coverage: $(EXE)

#Intended to clean up compilation artifacts but leave executable & coverage
clean:
	rm -f $(OBJ)
	rm -f $(MOD)
	rm -f $(COMP_DEP)

#Intended to reset directory to fresh state with no exe or artifacts
reset: clean
	rm -f *.gcno *.gcda
	rm -f -r $(EXE)*.dSYM
	rm -f ../$(EXE)*

#-----------------------------------------------------------------------------------
#Generics for source files
#-----------------------------------------------------------------------------------
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f90
	$(FC) -c $(OPT) $<

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f90
	$(FC) -c $(OPT) $<

$(EXE): $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ)
	mv ./$(EXE)$(TYPE)$(EXT) ../

#-----------------------------------------------------------------------------------
#Rules for entries in COMP_DEP. Don't forget to add them to make clean / reset
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#Dependency List
#Use [gfortran -M -cpp *.f90] repeatedly until clean compile to update rules below
#-----------------------------------------------------------------------------------
boundary_conditions.o boundary_conditions.mod: boundary_conditions.f90 \
 globals.mod
globals.o globals.mod: globals.f90 mesh_types.mod
main.o: main.f90 globals.mod read_gmsh.mod read_thrm.mod mesh_analyze_3d.mod \
 boundary_conditions.mod out_stats.mod tet_analyze.mod mesh_2d_construct.mod \
 mesh_analyze_2d.mod tri_analyze.mod
mesh_analyze_3d.o mesh_analyze_3d.mod: mesh_analyze_3d.f90 globals.mod
mesh_analyze_2d.o mesh_analyze_2d.mod: mesh_analyze_2d.f90 globals.mod
tet_analyze.o tet_analyze.mod: tet_analyze.f90 globals.mod
tri_analyze.o tri_analyze.mod: tri_analyze.f90 globals.mod
read_gmsh.o read_gmsh.mod: read_gmsh.f90 globals.mod
read_thrm.o read_thrm.mod: read_thrm.f90 globals.mod
mesh_2d_construct.o mesh_2d_construct.mod: mesh_2d_construct.f90 globals.mod
out_stats.o out_stats.mod: out_stats.f90 globals.mod
mesh_types.o mesh_types.mod: mesh_types.f90
