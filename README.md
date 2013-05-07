---
layout: default
title: AmBre ReadMe
---

===========
AmBre: Primer designing and sequence analysis for amplifying breakpoints.
===========

Developed by Anand D. Patel in Vineet Bafna's lab at University of California, San Diego
patel [dot] anandd [at] gmail [dot] com

AmBre is composed of two parts. A) Primer Designing
B) Long read sequence analysis.

    ambre_design <reference.fasta> <regions.txt> [<temptag>]
    ambre_analyze <reference.fasta> <contig> <mapped_reads.sam> [<temptag>]


Installing
=========

Software Requirements
-------------

AmBre requires Python 2.7.x and numpy

For primer designing, the following software needs to be installed

* Primer3 Version 2.3.0 http://sourceforge.net/projects/primer3/

* Blat Version 34x12 http://hgdownload.cse.ucsc.edu/admin/exe/

* Multiplx Version 1.2 http://bioinfo.ut.ee/?page_id=167

Sequence analysis requires PacBio's BLASR and SMRT pipe.
* BLASR https://github.com/PacificBiosciences/blasr
* PacBio SMRT-Analysis https://github.com/PacificBiosciences/SMRT-Analysis

Installing AmBre
-------------

AmBre can be run from source or installed.

Download the tar or zip ball and install using the setup.py script.
	
	python setup.py install

To configure the above binary dependencies, edit "ambre.conf" file. The config file
is located with

	ambre_test -c

Edit the file paths for primer3_dir, aligner_bin, and multiplx to locate the dependencies.
For primer3_dir set the file path to the directory containing primer3_core.
 In the config file specify number of threads and a temp directory to store intermediate files.
 Note, user-specific config files can be used by running any module with the --config argument. 
 
To verify the installation/configuration is correct and
 run a small test case for AmBre-design and AmBre-analyze (takes ~5min). ::

	ambre_test
	
If all is well, continue to usage.

Alternatively, to run from source directory use ::

	python -m ambre.test.install_unit -c
	python -m ambre.test.install_unit

	
Usage
=========

AmBre-design and AmBre-analyze are fairly independent.

For AmBre-design see [design usage](USEAGE-DESIGN.md) for details.
For AmBre-analyze see [analyze usage](USEAGE-ANALYZE.md) for details.

