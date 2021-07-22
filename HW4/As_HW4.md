### 1.ii.a) C=N-O Angle Terms in PCModel 10 output:


##### PCModel output `full_energy_printout.out`, Line 52:
```
	     At1      At2        At3       Angle     Thet0   Tconst    Ebend
Angle 1-3:  C(  1) -  N(  2) -  O(  3) :  111.276 180.000   0.200 =   60.952
```

### 1.ii.b) C=N-O-H Torsion Terms in PCModel 10 output:
	
##### PCModel output `full_energy_printout.out`, Line 87:	
```
	At1      At2      At3         At4             Types          Angle     V1      V2        V3         Etor
Tor:  C(  1) -  N(  2) -  O(  3) -  H(  6) :     2  108  145   21 :  179.999   0.000   0.000   0.000 =    0.000
```

### 2.i) Differences in the mengine output compared with PCModel 10 output.


First off, mengine prints out all of the parameters that were estimated in order to run the input molecule.

##### mengine output `pcmod.out`, Lines 2 - 10:

```
Following Parameters Were Estimated

bond         2  108     9.0000     1.2710   // equiv bond   2 37
dipole     108  145     0.2516  // estimated
bond       108  145     6.2457     1.3442   // estimated
angle        2  108  145     0.6950   110.0000     0.0000     0.0000   // estimated
angle       21  145  108     0.5430   110.0000     0.0000     0.0000   // estimated
torsion      2  108  145   21    0.000 +1   -2.400 -2    0.000 +3  // estimated values
torsion      1    2  108  145    0.000 +1    2.470 -2   -0.600 +3  // estimated values
```

Secondly, the estimated values are different from those used in PCModel 10. 

##### mengine output `pcmod.out`, Line 31:

```
	     At1      At2      At3     Angle     Thet0   Tconst     Ebend
Angle 1-3:  C(1  )-  N(2  )-  O(3  )  111.276   110.000  0.6950   = 0.0244
```

Equilibrium (minimum energy) bond angle (Theta_0 in degrees) between this combination of atoms making up a bond angle:

*Theta_0 is 110.000 instead of 180.000*

Bending constant (Tconst in md * Angstroms/rad^2) for this combination of atoms making up a bond angle:

*Tconst is 0.6950 instead of 0.200*

In-plane bending energy (Ebend in kcal/mol) for this combination of atoms making up a bond angle:

*Ebend is 0.0244 instead of 60.952*


##### mengine output `pcmod.out`, Line 66:
```
	At1      At2      At3     At4          Types           Angle     V1      V2        V3         Etor
Tor:   C(1  )-  N(2  )-  O(3  )-  H(6  )    2 108 145 21       179.98    0.000    -2.400   0.000    = -0.0000
```

Given the definition of V2 below (and the form of the torsion angle energy equation) the negative V2 value inverts the phase of the second order term, thus the minimum is at 90 degrees and teh maxima are at 0 and 180 degrees. The absolute value of the constant then determines the degree to which the potential energy curve is resembles the second order term. Indeed, when changing the V2 term to -8.40, the angle goes from 101.24 to 93.96 degrees.

*V2 is -2.400 instead of 0.000. V1 and V3 are 0.000 for both mengine and PCModel 10 estimated values.*

The torsional energy (Etor in kcal/mol) is described by an equation with three parameters, V1, V2, and V3.

V1 is the first order torsional constant. When positive, the energy maximum is at 0 degrees and a minimum is at 180 degrees - otherwise known as staggered/anti conformation in Newman projection terminology.

V2 is the second order torsional constant. When positive, the energy minima are at 0 and 180 degrees, and a maximum at 90 degrees. 

V3 is a third order torsional constant. When positive, the energy minima are at 60 and 180 degrees, while maxima are at 0 and 120 degrees.

### 2.ii) Difference in structures from PCModel 10 and mengine.

Atoms labeled with indices. 

Angles labeled in pink, lines point to the 1st and 3rd or 1st and 4th atom making up the bond angle or torsion angle, respectively. 

Colors: gray = carbon, white = hydrogen, blue = nitrogen, red = oxygen.

![image of pcmodel10](https://github.com/drmperez/HostDesigner_tutorials/blob/main/HW4/images/PCModel10_defaultMM3.png)

Figure above: Angle and torsion for structure optimized with PCModel 10 and default MM3 parameters of `dimeoxime.pcm.` 

![image of mengine](https://github.com/drmperez/HostDesigner_tutorials/blob/main/HW4/images/mengine_defaultMM3.png)

Figure above: Angle and torsion for structure optimized with mengine and default MM3 parameters of `dimeoxime.pcm.`

### 2.iii) Reason for difference between PCModel 10 and mengine results.

The default behavior of PCModel 10 estimates different values than mengine.

### 5. Compare values to average crystal structure values.

The values are as stated in problem 5 except for the torsion angle. The crystal structure average is 180 degrees while the computed value (in PCModel 10 with add.prm from mengine) is 101.43 degrees.

### 6. Adjust parameters.

Adjusting the parameters for the bond length, bond angle, and torsion angle resulted in a geometry with values close to the average crystal structure (see figure below for resulting values). 

User modified parameters:
```
bond       108  145     6.0000     1.4000   // edit - marilu
angle        2  108  145     8.0000   113.0000     0.0000     0.0000   // edit - marilu
torsion      2  108  145   21    0.000 +1    2.400 -2    0.000 +3  // edit - marilu
```


![image of pcmodel10 structure with optimized parameters](https://github.com/drmperez/HostDesigner_tutorials/blob/main/HW4/images/PCModel10_user_opt_MM3parameters.png)

