from pkg_resources import resource_filename, resource_stream
import os

class Config(object):  
  def __init__(self):
    self.config_name = 'ambre.conf'
    self.bin = dict()
    self.dir = dict()
    self.param = dict()

    # Defaults
    self.param['primer3_short'] = resource_filename(__name__, os.path.join('data','primer3_v1_1_4_default_settings.txt'))
    self.param['primer3_long'] = resource_filename(__name__, os.path.join('data','primer3_lr_A04.txt'))

  def add_bin(self, bin_name, bin_path):
    self.bin[bin_name]= bin_path
  def add_param(self, param_name, param_path):
    self.param[param_name] = param_path
  def load(self):
    for line in resource_stream(__name__, self.config_name):
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
  def valid(self):
    pass
  def replace(self):
    pass

CONFIG = Config()
CONFIG.load()

  