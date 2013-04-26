============
AmBre-design Usage
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

To Run
========

Designing requires a reference sequence (fasta format) and a set of forward and
reverse intervals to select primers from. The <regions.txt> is a tab-delimited
file having intervals in rows and columns as contig name, interval start,
interval end, type (True for forward, False if reverse). See
"examples/regions.test" for an example regions file.  If specified, resulting
primer designs are given in <temp_tag>.out. 

	python ambre_design.py <reference.fasta> <regions.txt> [<temptag>]

Parameters are defined in the the config file "ambre.conf". A user-defined
config file can be used with the "--config" argument.
Primer designing involves numerous parameters to ensure primers are viable for
PCR. AmBre requires primers to be compatible with one another and evenly tile
supposed regions. Specifying too large or too many regions requires more primers
and leads to unreliable PCR. See below for fine-tuning parameters.

Temptag is used as an identifier for intermediate files. 
It is recommended to specify temptag as
a directory or directory along with a prefix expected for AmBre-design
intermediate files.

For example, to reprint primer solutions and generate a summary figure for
primer locations on input regions using temptag.

	python ambre_design.py -c <reference.fasta> <regions.txt> <temptag>


Parameters
========

In the config file, modifying the following parameters
changes the primer designing task.

The approximate distance between primers.

	design_primer_spacing_p=6500

The density of primers in the candidate primer selection task,
that is (Nd/L-1) where **N** is the number of *desired* primers
**d** is the spacing between primers and **L** is the total 
spacing. See Patel et al. (submitted) for details. 

	design_primer_density_p=0.2

Candidate primer generation is performed with Primer3 (Rozen et al. 2000)
and the primer3 parameters file for long range PCR 
used in Patel et al. (submitted) is listed in the ambre/data directory.
Another primer3 parameter file can be used by adding the following
line to the config file.

	primer3_long_p=/<file>/<path>/<primer3>/<parameters.txt>


Fine-tune design parameters
---------

The number of candidate primers to request from primer3
per kilobase

	design_primer3_primers_per_bp_p=50

Filtering for primer3 penalty criteria

	design_max_primer3_penalty_p=1.2

Filtering primers with reference alignments within,
	design_max_amp_dist=20000


