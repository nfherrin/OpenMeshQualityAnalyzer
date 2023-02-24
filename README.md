# OpenMeshQualityAnalyzer
A utility to analyze the quality of a given 3D unstructured mesh according to skewness, smoothness, and aspect ratio.
Currently supports THOR mesh and Gmsh mesh formats.

Usage is done by calling the executable followed by the mesh filename, i.e.
```
./OpenMeshQualityAnalyzer.exe 3D_beavrs_2_out.thrm
```

Output is to a csv file with the same name as the input file but with `_stats.csv` appended.
Output gives the 3D mesh data which is total region volume, number of tets, tet volume average and standard deviation, tet skew (one minus the volume of the regular tet in the circumsphere of the original tet divided by the original tet volume) and standard deviation, tet aspect ratio (longest line in tet divided by shortest) and standard deviation, and finally tet smoothness (ratio of larger divided by smaller tet for a given tet with respect to all surrounding tets) and standard deviation.
Similarly, on each flat boundary face along a Cartesian direction those same values are reported for each triangular element.