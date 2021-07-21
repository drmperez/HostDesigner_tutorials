1.ii.a) C=N-O Angle Terms:

```
full_energy_printout.out, Line 52:
	     At1      At2        At3       Angle     Thet0   Tconst    Ebend
Angle 1-3:  C(  1) -  N(  2) -  O(  3) :  111.276 180.000   0.200 =   60.952
```

1.ii.b) C=N-O-H Torsion Terms:
	
```
full_energy_printout.out, Line 87:	
	At1      At2      At3         At4             Types          Angle     V1      V2        V3         Etor
Tor:  C(  1) -  N(  2) -  O(  3) -  H(  6) :     2  108  145   21 :  179.999   0.000   0.000   0.000 =    0.000
```

2.i) Differences with PCModel 10

```
pcmod.out, Line 31:
	     At1      At2        At3       Angle     Thet0   Tconst    Ebend
Angle 1-3:  C(1  )-  N(2  )-  O(3  )  111.276   110.000  0.6950   = 0.0244
```

Thet9 is 110 instead of 180.000
Tconst is 0.6950 instead of 0.200
Ebend = 0.0244

```
pcmod.out, Line 66:
	At1      At2      At3         At4             Types          Angle     V1      V2        V3         Etor
Tor:   C(1  )-  N(2  )-  O(3  )-  H(6  )    2 108 145 21       179.98    0.000    -2.400   0.000    = -0.0000
2.ii
```

V2 is -2.400 instead of 0.000. 
