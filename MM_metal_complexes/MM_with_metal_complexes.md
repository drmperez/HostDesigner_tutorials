**MM with metal complexes**

Included are SDF Molfiles for six metal complexes taken from the Cambridge Structural Database. 
In each case the atomic coordinates are the original x-ray coordinates.

![lewis structures of complexes taken from CSDB](https://github.com/drmperez/HostDesigner_tutorials/blob/014bc6ec17a4bae5c12f66c644449ce1bd8c05cb/MM_metal_complexes/HW%236.tif)

1. Using the information in the HD manual, pages 92-96, create pcmodel files for each of these complexes using points-on-a-sphere metal ions and optimize them with MM3 using mengine. So that there is no confusion regarding the identity of the ligands, I have included a single graphic showing a chemdraw sketches for them. When you run these, mengine will estimate parameters and produce the pcmod.par file. Look at this file to see what parameters were estimated. The input pcmodel files you create should contain original x-ray coordinates, do not subsequently overwrite them with optimized coordinates. This will allow you to perform a superposition of optimized structure versus crystal structure to see how well the model is working (see 2 below).

1. PCModel10 can be used to perform a superposition of two structures.

1. Open the original input file
2. Select Compare from the Display drop down menu
3. Use default setting of non-volatile atoms (he means non-hydrogen atoms) which will orient the next molecule to get the best superposition of heavy atoms
4. Hit the Next Str button to open the optimized structure
5. Hit the Calculate button
6. The superimposed structures are displayed. If you are using the default stick view, then will be easier to see the molecules by using the Draw Tubes option under the View menu. The Average Difference and the Root Mean Squared Differences for the pairwise displacements of all heavy atoms are shown in the right side pane
7. Done

Perform superpositions for all complexes. For comparison, below are the results that I obtained. If you are getting something else, then one of us has made an error. Check the HD User Manual and verify that all is correct with the input:
```
Complex   RMSD

BOHMOB    0.1087
JAVHET    0.1194
DUCQAT    0.2712
NEVTAK    0.1575
QELQAA    0.1785
YADTOM    0.2600
```

1. The MMFF94 contains default parameters for two of the metal ions found in these complexes. Using the electrostatic approach (only charge-charge and vdw interactions with the metal), make MMFF94 input files for these two complexes, optimize them, and perform superpositions. Which model gives the best result?
