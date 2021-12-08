---

title: "Using MM calculations to evaluate ligand structures"

---
## Using MM calculations to evaluate ligand structures

In order to maximize the binding interaction between a host and
guest, the host structure must be organized to bind the
guest. In an ideal situation the binding sites in the
host will match the binding sites on the guest (host
is complementary) and the lowest energy conformation of the host
is (1) the one that binds the guest and (2)
is rigid (host is preorganized). Molecular mechanics calculations can be
used to analyze the degree of complementarity and preorganization.

The degree of complementarity can be measured as the difference
in energy between the bound form of the ligand and
the binding form of the ligand. This energy difference is
computed by first optimizing the host-guest complex, second removing the
guest and computing a single point energy (E<sub>bound</sub>), 
third optimizing the bound form of the ligand to get the energy
of the binding form (E<sub>bind</sub>), and fourth 

∆E<sub>1</sub> = E<sub>bound</sub> – E<sub>bind</sub>. 

The lower the ∆E<sub>1</sub>, the more complementary the
ligand structure.

One measure of the degree of preorganization is the difference
in energy between the lowest energy form of the host
and the binding form of the host.  This energy
difference is computed by first conformer searching the host to
get the energy of the lowest energy form and second

∆E<sub>2</sub> = E<sub>bind</sub> – E<sub>low</sub>.  

The lower the ∆E<sub>2</sub>, the more preorganized the host structure.

The overall degree of host organization is the sum of
the two energy differences, 

∆E<sub>tot</sub> = ∆E<sub>1</sub> + ∆E<sub>2</sub>. 

Note that one can compute this value without consideration of
the binding conformer, 

∆E<sub>tot</sub> = E<sub>bound</sub> – E<sub>low</sub>.

By computing ∆E<sub>tot</sub> over a series of metal ions, for
example, over the lanthanide series, by plotting E<sub>tot</sub> as a
function of metal ionic radius, it is possible to determine
whether the ligand exhibits a steric preference for a specific
size metal and how strong this preference is.  Although
the calculation of ∆E<sub>tot</sub> is straight forward, 
these calculations become tedious when doing an entire series
of metal ions. For this reason, the utility code named scanme 
was created (see **insert link to scanme** for how to
make and install this code).

When scanme is run, it prompts the user for input
twice. The first prompt asks for the name of the
PCModel input file containing the metal-ligand structure. The second prompt
asks whether the user wants to conformer search the ligand.
The reason for this second question is that although some
ligands are not searchable due to excessive freedom, it is
still possible to evaluate the complementarity of the host binding
conformer. If the user chooses not to search the free
ligand, searchme will report only ∆E<sub>1</sub> for each metal and
if the user chooses to search the free ligand, searchme
will report ∆E<sub>tot</sub> for each metal. Note that if the
binding form is the same over the entire metal size
range, then these two methods will give the same shaped
energy versus radius plots, where the ∆E<sub>1</sub> plot is offset
from the ∆E<sub>tot</sub> plot by a constant ∆E<sub>2</sub> value.

To provide an example of how to use scanme, six
example PCModel input files are provided as Supporting Information (in
directory size_scan). The first training assignment is to determine and
compare intrinsic metal size preferences in bis-amine and bis-ether chelates
when the connecting link between the donor groups is varied
over 1,2-ethane, cis-1,2-cyclohexane, and trans-1,2-cyclohexane
links. To do this run scanme for each of the input structures
doing conformer searches for each case and graph the ∆E<sub>tot</sub>
values vs. metal ion radii for all six ligands on the same plot.
