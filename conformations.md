---

title: "Conformational analysis"

---


## Conformational analysis

In addition to providing very rapid geometry optimizations, molecular mechanics models
can also be used to perform conformer searches.  Such calculations can be performed using
either PCModel or mengine. In the latter case creating the conpcm file for a conformer
search is more complicated than for a geometry optimization. For this reason, the utility
code named *searchme* was created (see page 2 for how to make and install this code).  To
provide examples of what *searchme* does, some example PCModel input files are provided as
Supporting Information (in directory conformer_search). To run each example, enter the
command `searchme`, and give the name of the PCModel file when prompted.

The first example is n-butane. Although there are three C-C bonds in this molecule, only
one of these bonds is deemed searchable. This is because rotation of a terminal methyl
group, which is 3-fold symmetric, does not normally produce new conformers and is not
considered by the search algorithm.

The second example is a rigid, dimethylated naphthalene. The *searchme* code detects this
and stops with a message stating that searching this molecule is not necessary.

The third example is the prototype Nikki acyclic ligand. This molecule has 21 searchable
bonds. The *searchme* code detects this and stops with a message that this structure is
too flexible to search. This raises an important point.  When the conformational space
becomes too large, it is not possible to evaluate all of it in a doable time frame.
In this example, if it is assumed that that each rotatable bond has three possible
rotamers, then there are 321 = 3,486,784,401 possible conformers. The *searchme* code
will warn the user for any structure when there are ≥ 9 searchable bonds, ≥ 39 = 19,683
possible conformers, write the conpcm file and stop.  If the user wishes to proceed, then
they can start the search by running *mengine*. If there are more than 20 searchable bonds,
no *conpcm* file is written.
