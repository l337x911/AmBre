============
AmBre-analyze Usage
============

Description
========

Analyzes Pacifc Biosciences amplicon resequencing data where amplicons contain
structural variations (SV). AmBre-analyze clusters reads that support the same
SV and then calls breakpoint and consensus amplicon sequences, despite reads
having high insertion and deletion error rates. 

Sequence analysis involves a) local alignment with BLASR, b) alignment segmentation,
 c) geometric SV clustering
approach similar to Sindi et al. (2009). d) breakpoint and consensus sequence 
refinement. 
See Patel et al. (submitted) for details.

To Run 
======== 

AmBre-analyze requires a local alignment file (hard-clipped SAM
format) with the corresponding reference file (fasta format) and contig name.
Current implementation only supports intra-contig SVs. Estimated breakpoints
are reported in "<temptag>.bp" and estimated amplicon sequences are reported
in "<temptag>.fasta".

Given a "aligned_reads.bas.h5" file from PacBio.

	/blasr/alignment/bin/blasr aligned_reads.bas.h5 <reference.fasta> -clipping hard -sam -out <aligned_reads_h5.sam>

	ambre_analyze <reference.fasta> <contig> <aligned_reads.sam> [<temptag>]
	
or from source::

	python -m ambre.analyze.workflow <reference.fasta> <contig> <aligned_reads.sam> [<temptag>]

The output "<temptag>.bp" has each breakpoint as an entry with the following form ::

	#<name>\t\t\t<mode left bp>\t<mode right bp>
	<read_idx>\t<bp on frag idx>\t<left bp>\t<right bp>\t<d>
	...
	...
	=

The next step is to perform amplicon refinement using PacBio Amplicon
Resequencing Protocol.  Using SMRT-Analysis 4.0, call consensus sequencing
using the estimated amplicon "<temptag>.fasta" as the reference and the entire
read set.

Parameters are defined in the the config file "ambre.conf". A user-defined
config file can be used with the "--config" argument. See below for fine-tuning
parameters.

Temptag is used as an identifier for intermediate files.  It is recommended to
specify temptag as a directory or directory along with a prefix expected for
AmBre-analyze (can be the same prefix as AmBre-design) intermediate files. 


Dealing with amplicons with multiple breaks
--------

If independent breakpoints found in the previous section in AmBre-analyze
belong to the same amplicon, then annotate the breakpoint headers in the
"<temptag>.bp" file where instead of "#\t" is "#A549_01\t", "#A549_02", and
"..." to represent ordering of breaks along an amplicon. 

	python ambre_analyze.py --multi-break <reference.fasta> <contig> <aligned_reads.sam> [<temptag>]

The new "<temptag>.fasta" file will have updated amplicon sequences containing
multiple breaks.

Parameters
========

In the config file, modifying the following parameters
changes the sequence analysis task. See Patel et al. (submitted) for details.

Scoring for local alignment segmetnation
and filtering. Mismatch fraction is because CIGAR format only retains aligned bases (M)
 not misaligned.

	analyze_mismatch_fraction_p=0.02
	analyze_match_score_p=1
	analyze_mismatch_penalty_p=-3
	analyze_gapopen_penalty_p=-1
	analyze_gapext_penalty_p=-0.2
	analyze_breakpoint_penalty_p=-50

Filtering for minimum cluster size

	analyze_min_cluster_size_p=25

More filtering criteria for breakpoints

	analyze_min_sv_dist_p=10
	analyze_max_d_p=100
	analyze_default_segment_p=6500

