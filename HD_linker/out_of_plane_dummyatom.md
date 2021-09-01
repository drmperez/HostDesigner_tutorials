## Adding dummy atom to in-plane fragment

One way to add a dummy atom at 90 degrees from the 
plane of a planar structure is to use z-matrix
coordinate system. 

In MacMolPlt, the instructions are:

1. Open fragment with metal ion using MacMolPlt.
2. From the Subwindow menu, select Coordinates.
3. In the Coord. Type: field, select Z-Matrix
from the dropdown menu.
4. Click "Add" button to add another atom to the 
end of the list.
5. The atm 1, atm 2, and atm 3 indicate the atom
forming a bond, angle, and torsion (or Dihedral in 
MacMolPlt. For atm 1, choose the atom directly 
bonded to the guest (typically a nitrogen, or 
oxygen). For atm 2, select an atom bonded to 
atm 1. For atm 3, select an atom bonded to atm 2. 
6. Enter a Length, set Angle to 90.00 deg, 
and Dihedral to 90.00 deg.
7. Ensure that the bond orders are correct. To 
edit bond order, right click on a bond and select
the "Bond Order", choose the correct bond order
from the drop down menu.
8. Finally, from the File menu, select Export...
and save the file with MDL MolFile with "mol" 
extension.
