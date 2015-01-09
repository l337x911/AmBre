import sys

import numpy as np
import pandas as pd
from itertools import izip


''' Thermodynamics parameter file assumptions
Conditions given at T=37 degC in 1M NaCl

Enthalpy is given in kcal * mol^-1
Entropy is given in cal * degK^-1 * mol^-1
Gibbs given in kcal * mol^-1 

1 thermochemical cal = 4.184 J
1 IST1956 cal = 4.1868 J

0 degC = 273.15 degK

gas constant R
constR = 1.9872041 cal * K^-1 * mol^-1
'''

def conv_Cal_to_Joule(x): return x*4.1868
def conv_degC_to_degK(x): return x+273.15
constR = 1.9872041

PRECISION = 0.000001

'''
Estimates thermodynamic parameters for a 5'-3' oriented primer

Ignores Symmetry (oligonucleotide binds to self)

'''
class EstimateThermo(object):
	def __init__(self, nearest_neighbor_model=None):
		if nearest_neighbor_model is None:
			nearest_neighbor_model = 'santalucia1998.txt'
		
		self.nn_df = pd.read_csv(nearest_neighbor_model, sep='\t', header=0).set_index('PropSeq') 
		self.nn_df[['Enthalpy', 'Entropy', 'Gibbs']] = self.nn_df[['Enthalpy', 'Entropy', 'Gibbs']].astype(float)

		# append 5'->3' nucleotide pair abbrev as indicies
		abbrev = [k for k in self.nn_df.index if '/' in k and len(k)==5]
		#print self.nn_df.ix[abbrev,:].rename(index=lambda x:x[:2])
		abbrev_df = self.nn_df.ix[abbrev,:].rename(index=lambda x:x[:2])
		antiparallel_abbrev = [k for k in abbrev if k[:2]!=k[:-3:-1]]
		
		antiparallel_abbrev_df = self.nn_df.ix[antiparallel_abbrev,:].rename(index=lambda x:x[:-3:-1])

		self.nn_df = pd.concat([self.nn_df, abbrev_df, antiparallel_abbrev_df], axis=0)
			
	def _number_terminalAT(self, primer):
		return sum([primer[0]=='A' or primer[0]=='T', primer[-1]=='A' or primer[-1]=='T'])		

	def _compute_from_table(self, param_type, primer):
		# assumes param is either Gibbs, Entropy, Enthalpy
		param = self.nn_df[param_type]
		tot = param['Initiation']
		for nn_pair in izip(primer[:-1], primer[1:]):
			tot += param[''.join(nn_pair)]
		tot += self._number_terminalAT(primer)*param['TerminalAT']
		return tot

	def compute_enthalpy(self, primer):
		return self._compute_from_table('Enthalpy', primer)
	def compute_entropy(self, primer):
		return self._compute_from_table('Entropy', primer)

	def _derive_gibbs(self, primer, temperature):
		deltaH = self._compute_from_table('Enthalpy', primer)
		deltaS = self._compute_from_table('Entropy', primer)
		degK = conv_degC_to_degK(temperature)
		return deltaH - ( degK * deltaS / 1000.)

	def compute_gibbs(self, primer, temperature=37):
		# returns deltaG in kcal * mol ^-1
		if abs(temperature-37)>PRECISION:
			return	self._derive_gibbs(primer, temperature)
		return self._compute_from_table('Gibbs', primer)
	def compute_dissociation_constant(self, primer, temperature):
		deltaG = self.compute_gibbs(primer, temperature)
		norm_disK = -1000.0 * deltaG / conv_degC_to_degK(temperature)
		#print " derived from deltaG: Rln disK={0:0.2f} kcal * mol^-1 * degK^-1".format(norm_disK)
		return np.exp(norm_disK/constR)

	def compute_all(self, primer, temperature):
		deltaH = self._compute_from_table('Enthalpy', primer)
		deltaS = self._compute_from_table('Entropy', primer)
		deltaG = self.compute_gibbs(primer, temperature)

		print "primer: {0}\n deltaH: {1:0.3f} kcal * mol^-1\n deltaS: {2:0.3f} cal * degK^-1 * mol^-1\n deltaG: {3:0.3f} kcal * mol^-1".format(primer, deltaH, deltaS, deltaG)
		# deltaG = -RT ln disK
		norm_disK = -1.0 * deltaG / conv_degC_to_degK(temperature)
		print " derived from deltaG: Rln disK={0:0.2f} kcal * mol^-1 * degK^-1".format(norm_disK)

		disK = np.exp(1000*norm_disK/constR)
		print " derived from deltaG: disK={0:0.4e}".format(disK)
		print " @{0:0.2f} degC".format(temperature)

	def _fraction_target_bound(self, oligo_c, disK):
		return oligo_c / (oligo_c + disK)
	def plot_fraction_target(self, primer, temperature):
		from matplotlib import pyplot as plt
		disK = self.compute_dissociation_constant(primer, temperature)
		
		oligo_c = np.logspace(-3,35,num=30, base=10)
		f = self._fraction_target_bound(oligo_c, disK)

		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(oligo_c, f)
		ax.set_xscale('log')
		ax.set_xlabel('Oligo []')
		ax.set_ylabel('fraction of target bound')
		ax.set_ylim(0,1)
		plt.show()


def test():
	w = EstimateThermo()
	print 'From Table, Gibbs:', w.compute_gibbs('CGTTGA')
	print 'From Deriv, Gibbs:', w._derive_gibbs('CGTTGA', 37.0)
	w.compute_all('CGTTGA', 37.0)
	w.plot_fraction_target('CGTTGA', 37.0)

	w.compute_all('TTTCTCTCTTAGATTGGAATAATTGGTGGAAC', 25.0)
	w.compute_all('TTTCTCTCTTAGATTGGAATAATTGGTGGAAC', 37.0)
	w.compute_all('TTTCTCTCTTAGATTGGAATAATTGGTGGAAC', 64.0)
	w.compute_all('TTTCTCTCTTAGATTGGAATAATTGGTGGAAC', 94.0)
	#w.plot_fraction_target('TTTCTCTCTTAGATTGGAATAATTGGTGGAAC', 64.0)


if __name__ =='__main__':
	#primer = sys.argv[1]
	test()
