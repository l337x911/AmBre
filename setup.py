from setuptools import setup, find_packages


version = '1.0'

setup(name='AmBre',
      version=version,
     
      description="""Designing primers for and analyzing amplified breakpoints.
Designs primers using a workflow 1) primer generation with Primer3, 2) Primer filtering with Blat, 3) Simulated annealing for primer selection. Analyzes PacBio sequenced amplicons with 1) BLAST local alignment, 2) Dynamic programming to filter and trim local alignments, 3) Sweep-line algorithm for clustering fragments, 4) Summarize breakpoints into amplicon templates""",
      long_description=open('README.md').read(),
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='AmBre PAMP primer PCR SV cancer amplification breakpoints',
      author='Anand D. Patel',
      author_email='patel.anandd@gmail.com',
      url='',
      license=open('license.txt', 'rb').read(),
      packages=find_packages(exclude=['ez_setup',]),
      package_data={'ambre': ['*.txt', '*.md', 'ambre.conf', 
                              'examples/reference.fasta', 'examples/reference.fasta.fai'
                              'examples/aligned_blasr_h5.sam',
                              'examples/regions_ex/*',
                              'examples/regions.test', 'data/*.txt', 'examples/gold/*']},
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
          "numpy>=1.6"
      ],
      extras_require = {
        'MPL' : ["matplotlib>=1.1.1rc"]                 
      },
      entry_points = {
        'console_scripts': [ 
        'ambre_design = ambre.design.workflow:main',
        'ambre_analyze = ambre.analyze.workflow:main',
        'ambre_test = ambre.test.install_unit:install_check',
        ]
      }
      )
