**RESULTS: MM with metal complexes**
```
Complex   RMSD (Hay)	RMSD (MPG)

BOHMOB    0.1087	0.1085
JAVHET    0.1194	0.1195  Note: on PCModel 10, atom type 318 was not allowed. Must change in input file.
DUCQAT    0.2712	0.2838	Note: Gd labeled as W in PCModel 10 but has correct atom type for Gd of 355.
NEVTAK    0.1575	0.1466 Note: square planar structures need to add dummy "atoms" to get the correct structure. This example can't make tetrahedral.
QELQAA    0.1785	0.1785
YADTOM    0.2600	0.1321

Hay = Ben Hay
MPG = Marilu Perez Garcia
```

*The MMFF94 contains default parameters for two of the metal ions found in these complexes. Using the electrostatic approach (only charge-charge and vdw interactions with the metal), make MMFF94 input files for these two complexes, optimize them, and perform superpositions. Which model gives the best result?*

**Results**

JAVHET: No explicit N-Metal bonds in input file.

RMSD with PCM10 0.3534

BOHMOB:

RMSD with PCM10 0.1972

DUCQAT: optimizing in pcm10 after H-A/D (H-add/delete) made program hang.

RMSD with PCM10 N/A 

NEVTAK: 

RMSD with PCM10 0.2536

QELQAA:

RMSD with PCM10 0.2851


 
