===========
AmBre: Primer designing and sequence analysis for amplifying breakpoints.
===========

Developed by Anand D. Patel in Vineet Bafna's lab at University of California, San Diego
patel [dot] anandd [at] gmail [dot] com

AmBre is composed of two computational parts. A) Primer Designing
B) Long read sequence analysis.::

    python ambre_design.py <reference.fasta> <regions.txt> [<temptag>]
    python ambre_analyze.py <reference.fasta> <contig> <mapped_reads.sam> [<temptag>]


Installing
=========

Software Requirements
-------------

AmBre requires Python 2.7.x and numpy

For primer designing, the following software needs to be installed

* Primer3 Version 2.3.0 http://sourceforge.net/projects/primer3/

* Blat Version 34x12 http://hgdownload.cse.ucsc.edu/admin/exe/

* Multiplx Version 1.2 http://bioinfo.ut.ee/?page_id=167

Sequence analysis using PacBio requires BLASR and SMRT pipe.
* BLASR https://github.com/PacificBiosciences/blasr
* PacBio SMRT-Analysis https://github.com/PacificBiosciences/SMRT-Analysis

Installing AmBre
-------------

Download the tar or zip ball and install using the setup.py script. ::
	
	python setup.py install

To configure the above binary dependencies, edit "ambre.conf" file. The config file
is located with ::

	python ambre_test.py -c

Edit the file paths primer3_dir, aligner_bin, and multiplx to point to the correct
dependency folder (primer3) or binary (blat and multiplx). In the config
file specify number of threads and a temp directory to store intermediate files.
 Note, user-specific  config files can be used by running any module with the --config argument. 
 
To verify the installation/configuration is correct and
 run a small test case for AmBre-design and AmBre-analyze (takes ~5min). ::

	python ambre_test.py
	
If all is well, continue to usage.

Usage
=========

AmBre-design and AmBre-analyze are fairly independent.

For AmBre-design see [design usage](USEAGE-DESIGN.md) for details.
For AmBre-analyze see [analyze usage](USEAGE-ANALYZE.md) for details.

