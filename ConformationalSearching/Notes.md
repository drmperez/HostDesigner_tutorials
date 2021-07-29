## Compute times

1. Longest molecule, undecane, ran for 8 min.
2. Second longest, decane, ran for about 3 min.
3. Nonane ran for about 1 min.
4. Octane and shorter ran for less than a min.

## Cyclohexane

1. Templates -> C6
2. Compute -> conformer search
3. Add ring -> accept
4. Start
5. Runs for a while but no "conformational search done" message pops up.
6. Compute -> conformer search
7. Add ring -> accept
8. Start
9. Runs ok.

10. Repeat 1-9 and found exact same behavior.

## Hexane_2

1. Used Build feature
2. added carbons to both sides of ethanyl group.
3. Compute -> conformational
4. Set up Bonds
5. Default values are 0 1, 0 4, and 1 2.
6. Conformer search results in bent structure - see hexane_2.png.

## Pentane

1. build
2. add c's
3. Compute -> conformations
4. set up bonds
5. add default
6. change one bond from 2 3 to 3 4
7. Run and got fully extended pentane.


# Errors/Bugs/Suggestions

1. Help menu has inconsistent fonts.
2. Can text in help menu be increased?
3. The text refers to a "save job" button but only "Write job file" button exists.
4. Can more bonds be added? Other than the default ones which are always 2 bonds less
than the number of bonds in the chain.
5. Editing bonds is unclear and I'm not sure if I truly edited them. Looking at the job
file doesn't clarify.
6. When adding atoms to both sides of the ethanyl group using "Build" feature, the
numbering in the "Setup Bonds" window are odd. See Notes on Hexane_2 above.
7. Many structures complete with a geometry that is not (or should not be) the lowest
energy conformation.
8. The lowest energy conformation seems to be found when modifying the default
bond values originally populated in the "Bond Setup" menu.