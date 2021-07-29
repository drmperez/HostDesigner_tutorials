# Conformational Searching

GLOBAL-MMX (GMMX) is a steric energy minimization program which uses the currently selected force 
field to search for the global energy minimum and other low energy local minima.  Processing
by GMMX is done in two stages.  The first stages randomly searches over the selected rings and
rotatable bonds and keeps all conformers within 3.5 Kcal of the lowest energy conformer found
(Eminim) during the minimization.  The second cycle reminimizes the structures found in the
first cycle and keeps only those which are within 3.0 Kcal of the lowest energy conformer.  The
default energy windows for the two cycles can be changed.  In addition to an output file containing
the coordinates of the final structures, a textual summary file called '<OUTPUT_FILENAME>.pkm'
is produced which lists the energy  and Boltzmann distribution of each low energy conformation.
This file may also include a listing of query operations (PMR coupling constants, distances,
angles, and dihedrals) for each structure as well as a Boltzmann averaged summary.

Although the conformational searching techniques that follow are unique in their approach, the methods
described by the Still group (see M. Saunders, K.N. Houk, Y-D Wu, W.C. Still, M. Lipton, G.
Chang, and W. Guida J. Am. Chem. Soc. 112, 1419 (1990) and references to previous work cited
there in.) were the inspirational source for this work.  The search techniques in GMMX are based
on the methods used in BAKMDL, developed by Professor Kosta Steliou of Boston University, and
ported to the MMX force field by Mark Midland of UC Riverside and Joe Gajewski and Kevin Gilbert
of Indiana University.  We wish to thank Professors Still and Steliou for sharing with us unpublished
work and code that greatly enhanced our routines.

Conformational space can be searched by GMMX in three ways.  These are the Mixed method, the Bonds method and
the Cartesian method.  The Bonds method randomly selects a subset of the bonds designated by
the user for rotation.  Bond rotation can cause large changes in the shape of a molecule.  The
Cartesian method randomly moves a subset of all heavy and nonvolatile atoms in 3d space causing
small changes in the shape of the molecule.  The Mixed method alternates between the Bonds 
and Cartesian methods and is the preferred method.  This procedure is especially efficient
with cyclic structures having side chains.  You may select ring bonds for rotation as well
as the side chain bonds but restricting the bond selection to the side chain bonds usually
is all that is necessary if only 3 to 7 member rings are present. The coordinate (Cartesians)
movements will apply to the entire structure.
