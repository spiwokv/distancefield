# Distancefield

## Compile:
`gcc -o distancefield.o distancefield.c -lm`

## Usage
`./distancefield.o protein.pdb 35 35 35 2 15 21 14 > distances.txt`

Protein must be centered in the periodic box with coordinates x = 0 &mdash; 70, y = 0 &mdash; 70 and z = 0 &mdash; 70
(35x2 in Angstroms). Grid will be separated by 2 Angstroms (total number of grid points must not exceed 100,000!).
Ligands must be removed. Proteins containing non-proteinous elements (e.g. P) are not suported.
Element symbol (N, C, O or S) must be in the column number 14. Only atoms strating by ATOM will be 
considered. The binding site is located in x_i = 15, y_i = 21 and z_i = 14 (i.e. [31.0,43.0,29.0] Angstroms).

The code will separate the box by a grid with 35x35x35 (for the above presented seting). Points coliding with the protein will
be ignored. Next it will calculate a minimal distance of each point from the binding site using Dijkstra's algorithm 
(will take couple of minutes) and print them to the output.

