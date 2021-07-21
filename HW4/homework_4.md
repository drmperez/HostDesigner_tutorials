# Homework 4

Study pages 78 thru 167 in the MM3 manual.  Also, study pages 91 – 105 and 112 – 115 of the 
HostDesigner manual.  

This assignment uses dimethyloxime as an example.  If you recall, this structure exhibited 
some squirrely behavior upon minimizing in the MM3 model. 

Try it again.  Open this dimethyloxime.pcm file in PCModel10, optimize the structure with 
MM3, and what you see is that the C=N-O angle goes linear and the C=N-O-H torsion goes to 
90 deg.  When you see behavior like this, either the parameters are incorrectly assigned 
or else they are missing entirely.

What do you do?

### 1.	Investigate to see what parameters are being used for this molecule. In PCModel10:

	a. Open the dimethoxyoxime.pcm file, 

	b. select MM3 under the ForceField menu, 

	c. under the Compute menu, click on Full Energy Printout, 

	d. under the Compute menu, click on Minimize, 

	e. under pcm10.app, click on Quit pcm10.app, 

	f. click on Save Log File in the popup window and you will be prompted to give \
\	a file name and location for the log file, 

	g. open and read the log file with a text editor.

Look at this log file.  Given the MM3 user manual, can you make sense out of what is being 
displayed on each interaction line?  Normally you will not be messing with van der Waals, 
charge, dipole, and cross-term (stretch-bend, bend-bend, stretch-torsion) parameters.  
Thus, you want to focus on the stretching, bending, and torsion interactions.  
These are listed under Bond Terms, Angle Terms, and Torsion Terms. 

* What is the code using for the C=N-O angle?
* What is the code using for the C=N-O-H torsion angle?

### 2.	Next, try to do this with MENGINE.  You can get the analogous full printout by adding 
the following line to the conpcm file:

	iprint 1

After the calculation runs, the pcmod.out file will contain a full printout of all 
interactions. Compare the contents of the pcmod.out file with the log file from PCModel10.  

* Are there any differences?
* Does optimized geometry look different from what PCModel10 produced?
* Why is it different?

The behavior in MENGINE is different than in PCModel10.  The difference is that MENGINE 
will attempt to estimate any missing parameters and, if it estimates any, then they will 
be used in the calculation and written to a file named pcmod.par.  If a file named 
pcmod.par is produced by MENGINE, then this tells you that at least one parameter was 
missing for the structure and if you view pcmod.par with a text editor you can see 
what was missing and what was assigned.


### 3.	Now that you see what was missing, you could read in an added parameter file that
contains the missing parameters prior to optimization.  To do this with MENGINE, rename the
pcmod.par file to add.prm and include the following line in the conpcm file:

	addpar add.prm

Now run MENGINE again.  This time pcmod.out should show the same parameters used as before, 
but now you should not produce a pcmod.par file because we read in the missing parameters and
nothing was missing.

Note that you can name the added parameter file anything you like.  Here the name add.prm is
arbitrary and was used so it would be different from pcmod.prm.  This was so you can see
that the pcmod.par is not produced when the add.prm file is read in.

Note:
You can create your own added parameter file by copying any parameter line from a 
default parameter file, either mm3.prm or mmff94.prm, and pasting it into the file.
In other words, the format for each type of interaction in this added parameter 
file is the same as the format seen in the default parameter file.

Try this out: Make an ethane molecule and optimize it with MM3.  Create an added parameter 
file in which you change the default C-C distance from 1.5247 angstroms to 1.600 angstroms.  
Optimize the structure again using the added parameters and verify that the C-C 
distance gets longer.

### 4.	You can also read in added parameters when using PCModel10.  Try it.  After getting set to minimize the structure, but before you hit Compute/Minimize, go to the Options menu and select Use Added Parameters.  A file open window will appear asking you to select the added parameter files.  Select the add.prm file and optimize the structure.  Is the resulting structure different?  Examine the log file and verify that the added parameters are now being applied.

### 5.	We want to verify that the parameters are giving a structural features that are consistent with experiment.  One way to do this is compare bond lengths, bond angles, and dihedral angles with average values observed in the Cambridge Database.  Here is what the average experimental data is from 1118 dialkylated oxime structures says:

- C=N distance 1.28 +/- 0.02

- N-O distance 1.40 +/- 0.02

- C-N-O angle 113 +/- 2 deg

- C=N-O-H dihedral angle is 180 +/- 7 deg

Use the Query feature in PCModel10 to compare these values to what the MM3 model is currently giving when using add.prm.  What you should see is that the C=N distance is pretty good, the N-O distance is computed too short, the C=N-O angle is computed too large, and the C=N-O-H torsion angle is correct.

(6)	 You could decide to adjust the parameters to get a better agreement with experiment.  The PCModel10 interface provides a fast way to do this.  The approach is to open the input file, choose the use added parameters, optimize the structure, and query the N-O and C=N-O angles so that you can see the values on the screen.  

Open a terminal window and open the add.prm file in a text editor.  Modify one or more parameters by changing their values.  Save the add.prm file.  Tell PCModel10 to use this modified version of the added parameters.  You have to first toggle the Use Added Parameters off, pick it again, and then choose the updated add.prm file.  Now minimize the structure.  The query values should change showing the effect of the changed parameters.  Iteratively modify the stretch and bend parameters to achieve a better agreement with the experimental data.

You could also do the parameter adjustment iterations using MENGINE and PCModel10.  The process would be (a) edit added parameter file, (b) run MENGINE, (c) open output file in PCModel10, (d) query the structural features, (e) repeat.  Note that you should only need to query the structural features once because if you leave PCModel10 with the queried structure showing on the screen, then if you open the structure again, the queries will still be applied and should show the updated values.  Although the above process sounds like it might be tedious, actually it can actually go quickly.  You can keep one terminal window open with the add.prm file displayed.  After you modify and save changes, keep this file open to see the current parameters and to be ready for the next changes.  Have another terminal window open to allow you to run MENGINE.  Finally, have your environment set so that double clicking the output pcm file produced by MENGINE will open it in PCModel10.  

What parameter values do you finally come up with?

