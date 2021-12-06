===23/02/2012=====================================
GLASCOW 
==================================================

These directories contains sample data in order to test GLASCOW

==================================================
EX1
==================================================

These data consist of 81 genotyped individuals are considered, the
phenotype is a binary trait.  Correlation matrix is computed over the
genome. The parameter file is restricted to the mandatories keywords.

==================================================
EX2
==================================================

These data consist of 308 genotyped individuals. Phenotype is a binary
trait, only 1469 markers are provided in the genotype data, but
correction for stratification is done through the file matrix.txt
(formerly computed with all genome wide SNP data). SCORE and REMLS
options are used for illustrative purposes.

==================================================
EX3
==================================================

These data consist of 762 genotyped individuals. Phenotype is a normal
traits.  Two fixed effects are included in the model. Only 981 markers
data are provided in genotype file, but correction for population
stratification is done through the file matrix.txt (formerly computed
using pedigree). REMLI, REMLT and REMLS options are
used for illustrative purposes.

==================================================
EX4
==================================================

These data are the same as EX1 except that genotypes are splitted in 5
files. This allow to use the parallel i/o of GLASCOW.
