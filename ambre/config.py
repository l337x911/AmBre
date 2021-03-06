'''
#  ambre.config.py
#
#  Copyright March 2013 by Anand D. Patel
#
#  This program is free software; you may redistribute it and/or modify its
#  under the terms of the GNU General Public License as published by the Free
#  Software Foundation; either version 2 of the License or
#  any later version.
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  For any other inquiries send an email to: Anand D. Patel (adp002@ucsd.edu)
'''

from pkg_resources import resource_filename
import os, sys

class Config(object):  
  def __init__(self):
    self.config_fpath = resource_filename(__name__, 'ambre.conf')
    self.bin = dict()
    self.dir = dict()
    self.param = dict()

    # Defaults
    #self.param['primer3_short'] = os.path.abspath(resource_filename(__name__, os.path.join('data','primer3_v1_1_4_default_settings.txt')))
    #self.param['primer3_long'] = os.path.abspath(resource_filename(__name__, os.path.join('data','primer3_lr_A04.txt')))
    self.param['primer3_long'] = os.path.abspath(resource_filename(__name__, os.path.join('data','primer3_lr_A06.txt')))
    #self.param['primer3_long'] = os.path.abspath(resource_filename(__name__, os.path.join('data','primer3_rna_A01.txt')))
    
    self.dir['examples'] = None
  def add_bin(self, bin_name, bin_path):
    self.bin[bin_name]= bin_path
  def add_param(self, param_name, param_path):
    self.param[param_name] = param_path
  def load(self):
    with open(self.config_fpath) as f:
      for line in f:
        if line.startswith('#'): continue
        tokens = line.strip().split('=')
        try:
          key_name, key_path = tokens
        except:
          continue
        if key_name.endswith('_bin'):
          self.bin[key_name[:-4]] = key_path
        elif key_name.endswith('_p'):
          self.param[key_name[:-2]] = key_path
        elif key_name.endswith('_dir'):
          self.dir[key_name[:-4]] = key_path
  def update(self, config_fpath):
    try:
      assert os.path.isfile(config_fpath)
      self.config_fpath = config_fpath
      self.load()
    except AssertionError:
      print "Error: Config path %s does not exist."%config_fpath
      sys.exit(1)

CONFIG = Config()
CONFIG.load()

  
