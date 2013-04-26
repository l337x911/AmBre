============
AmBre-design usage
============

Description
========

AmBre-design designs an assay to detect structural variations. Given input
regions where one expects breakpoints for a particular structural variation
(SV), AmBre-design finds primers uniformly tiled across the regions to reliably
amplify DNA harboring the structural variation with PCR.
 
Each structural variation is a rearrangement of a reference DNA segments to
produce a novel adjacency in the donor sample. However, it is possible to know
that a SV occurs, without knowing the exact boundaries of the rearranged
reference segments. In PCR, a forward primer (DNA sequence) designed to the
left and a reverse primer (reverse complement DNA sequence) designed to the
right of the donor sample adjacency would reproduce DNA harboring the SV. PCR
is only capable of reproducing a limited DNA length, so primers can only be
designed a certain distance away from the donor adjacency. Using multiple
forward primers and reverse primers across a breakpoint region ensures that
some forward and reverse primer appear within an amplifiable distance around
the donor sample adjacency.

See Patel et al. (submitted) for details.

Running
========

Designing requires a reference sequence (fasta format) and a set of forward and reverse intervals
to select primers from. The <regions.txt> is a tab-delimited file
having intervals in rows and columns as contig name, interval start, interval end, type (True for forward, False if reverse). See "examples/regions.test" for an example regions file.
 If specified, resulting primer designs are given in
<temp_tag>.out. 

	python ambre_design.py <reference.fasta> <regions.txt> [<temptag>]

Primer designing involves numerous parameters to ensure primers are viable
for PCR. AmBre requires primers to be compatible with one another and
evenly tile supposed regions. Specifying too large or many regions requires
more primers and inefficent PCR. See below for fine-tuning parameters.

Parameters are defined in the the config file "ambre.conf". A
user-defined config file can be used with the "--config" argument.

Temptag is used as an identifier for intermediate files. Recommended
to specify a directory or directory along with a prefix expected for 
AmBre-design intermediate files.

For example, to reprint primer solutions and generate
a summary figure for primer locations on input regions::

	python ambre_design.py -c <reference.fasta> <regions.txt> <temptag>



Parameters
========

