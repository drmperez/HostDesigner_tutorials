## Making mengine input files

The mengine code is a command line executable that does molecular mechanics calculations.
It is essentially the computational part of the old PCModel software, version 9.
Because of this, mengine reads PCModel formatted input files. 
HostDesigner uses mengine to do molecular mechanics post-processing using either the MM3 or MMFF94 force fields.
This topic is covered in detail in the HostDesigner User’s Manual, Section 6.

The procedure is as follows:  

(a) HostDesigner writes a PCModel input file for the structure and writes a conpcm file telling the mengine what to do, 
(b) HostDesigner calls mengine, 
(c) mengine does the calculation and writes an output file containing structure and energy information, 
(d) HostDesigner reads the mengine output file.  In order for this post-processing procedure to work correctly, the user must assign force field atom types in the HostDesigner input fragment files.  

PCModel input files are provided for a series of structures using either MM3 atom types and, where possible MMFF94 atom types, in order to assist with learning how to assign atom types (supporting information). These structures are shown in the graphics below.  The learning process involves (a) using PCModel to create the structure with a correct Lewis structure, (b) write a PCModel formatted file with either MM3 or MMFF94 atom types, and (c) edit this file to check and/or correct the atom type assignments to make sure they are correct.   The learner can check their version against the provided reference PCModel files.  The user should refer to the HostDesigner User Manual as a guide for how to assign types.  In addition to atom typing, note that the MM3 model requires the correct marking of pi systems whereas the MMFF94 model requires the correct specification of bond orders.  If you are working with metal complexes, then you will probably be using the MM3 model.  Note that one should never assume PCModel atom type assignments are correct.  They should always be checked and this is especially true for metal complexes.

