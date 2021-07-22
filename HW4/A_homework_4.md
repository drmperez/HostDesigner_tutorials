### 1.ii.a) C=N-O Angle Terms:


##### full_energy_printout.out, Line 52:
```
	     At1      At2        At3       Angle     Thet0   Tconst    Ebend
Angle 1-3:  C(  1) -  N(  2) -  O(  3) :  111.276 180.000   0.200 =   60.952
```

### 1.ii.b) C=N-O-H Torsion Terms:
	
##### full_energy_printout.out, Line 87:	
```
	At1      At2      At3         At4             Types          Angle     V1      V2        V3         Etor
Tor:  C(  1) -  N(  2) -  O(  3) -  H(  6) :     2  108  145   21 :  179.999   0.000   0.000   0.000 =    0.000
```

### 2.i) Differences with PCModel 10


##### pcmod.out, Lines 2 - 10:

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

##### pcmod.out, Line 31:

```
	     At1      At2        At3       Angle     Thet0   Tconst    Ebend
Angle 1-3:  C(1  )-  N(2  )-  O(3  )  111.276   110.000  0.6950   = 0.0244
```

Thet9 is 110 instead of 180.000

Tconst is 0.6950 instead of 0.200

Ebend = 0.0244


##### pcmod.out, Line 66:
```
	At1      At2      At3     At4          Types           Angle     V1      V2        V3         Etor
Tor:   C(1  )-  N(2  )-  O(3  )-  H(6  )    2 108 145 21       179.98    0.000    -2.400   0.000    = -0.0000
```

V2 is -2.400 instead of 0.000. 

### 2.ii) Difference in structures from PCModel 10 and mengine.

Atoms labeled with indices.
Colors: gray = carbon, white = hydrogen, blue = nitrogen, red = oxygen.

![image of pcmodel10](https://github.com/drmperez/HostDesigner_tutorials/blob/main/HW4/images/PCModel10_defaultMM3.png)

Figure: Angle and torsion for structure optimized with PCModel 10 and default MM3 parameters of `dimeoxime.pcm.` 

![image of mengine](https://github.com/drmperez/HostDesigner_tutorials/blob/main/HW4/images/mengine_defaultMM3.png)

Figure: Angle and torsion for structure optimized with mengine and default MM3 parameters of `dimeoxime.pcm.`
