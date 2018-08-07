# Distancefield

## Compile:
`gcc -o distancefield.o distancefield.c -lm`

## Usage
`./distancefield.o protein.pdb 70 70 70 2 10 10 10 > distances.txt`

Protein must be centered in the periodic box with coordinates x = 0 -- 70, y = 0 -- 70 and z = 0 -- 70 (in Angstroms).
Grid will be separated by 2 Angstroms (low value may cause huge number of grid points leading to lack of memory!).
The binding site is located in x_i = 10, y_i = 10 and z_i = 10 (i.e. [21,21,21] Angstroms).

The code will separate the box by a grid with 35x35x35 (for the above presented seting). Points coliding with the protein will
be ignored. Next it will calculate a minimal distance using Dijkstra's algorithm and print them to the output.

