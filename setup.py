from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='AmBre',
      version=version,
      description="Designing primers for and analyzing amplified breakpoints",
      long_description="""\
Designs primers using a workflow 1) primer generation with Primer3, 2) Primer filtering with Blat, 3) Simulated annealing for primer selection. Analyzes PacBio sequenced amplicons with 1) BLAST local alignment, 2) Dynamic programming to filter and trim local alignments, 3) Sweep-line algorithm for clustering fragments, 4) Summarize breakpoints into amplicon templates""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='AmBre PAMP primer PCR SV cancer amplification breakpoints',
      author='Anand D. Patel',
      author_email='patel.anandd@gmail.com',
      url='',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=True,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points = {
        'console_scripts': [ 
        'ambre_design = ambre.design:main',
        'ambre_analyze = ambre.analyze:main',
        ]
      }
      )
