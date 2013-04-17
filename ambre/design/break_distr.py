'''
#  ambre.design.break_distr.py
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

from matplotlib import pyplot as plt
import numpy as na
from itertools import product

plt_xmin, plt_xmax = 10000,15000

def get_length_distribution(fs,rs):
	max_l = max(fs)+max(rs)
	increasing_triangle = na.arange(max_l, dtype=na.int64)
	length_distr = na.zeros(max_l+1, dtype=na.int64)
	
	for f,r in product(fs,rs):
		min_d, max_d = min(f,r), max(f,r)
		diff_d = abs(f-r)
	
		length_distr[:min_d+1] += increasing_triangle[0:min_d+1]
		if diff_d>0:
		  length_distr[min_d+1:min_d+1+diff_d] += min_d
		length_distr[min_d+1+diff_d:2*min_d+diff_d] += min_d-increasing_triangle[1:min_d]
		
		assert f*r== na.sum(increasing_triangle[0:min_d+1])+ na.sum(min_d-increasing_triangle[1:min_d]) + min_d*diff_d
		
	return length_distr

def process_primers(regions_fpath, primers_fpath):
	primers_f, primers_r = [],[]
	with open(primers_fpath, 'rb') as f:
		# quench header
		f.next()
		for line in f:
			tokens = line.strip().split('\t')
			pos = int(tokens[1])
			if tokens[2]=='True':
				primers_f.append(pos)
			else:
				primers_r.append(pos)
	
	delimiters_f, delimiters_r = None,None
	with open(regions_fpath, 'rb') as f:
		for line in f:
			tokens = line.strip().split('\t')
			if tokens[-1] == 'forward':
				delimiters_f = int(tokens[1])+int(tokens[2])
			elif tokens[-1] =='reverse':
				delimiters_r = int(tokens[1])
				
	primers_f.append(delimiters_f)
	primers_r.append(delimiters_r)
	primers_f.sort()
	primers_r.sort()
	primers_f, primers_r = na.array(primers_f, dtype=na.int64), na.array(primers_r, dtype=na.int64)
	primers_f_d = primers_f[1:]-primers_f[:-1]
	primers_r_d = primers_r[1:]-primers_r[:-1]
	return (primers_f_d,primers_r_d)

def plot_lengths(distr, label):
	x = na.arange(len(distr))
	#y = na.array(distr,dtype=na.float)/na.sum(distr)
	y = na.array(na.cumsum(distr),dtype=na.float)/na.sum(distr)
	
	plt.plot(x[plt_xmin:plt_xmax],y[plt_xmin:plt_xmax], label=label)

spacings = process_primers("pamp_regions_A5.out","a6_regions_rho0.2_d06500_1.best")
observed_distr = get_length_distribution(*spacings) 

forward_spacing, reverse_spacing = na.sum(spacings[0])/len(spacings[0]), na.sum(spacings[1])/len(spacings[1])

ideal_spacings = ([forward_spacing]*len(spacings[0]), [reverse_spacing]*len(spacings[1]))
expected_distr = get_length_distribution(*ideal_spacings)

print "Design vs Ideal (f_d = %05d, r_d = %05d)"%(forward_spacing, reverse_spacing)

print "Forward spacing: %06d %06d"%(na.sum(spacings[0]), na.sum(ideal_spacings[0]))
print "Reverse spacing: %06d %06d"%(na.sum(spacings[1]), na.sum(ideal_spacings[1]))
print "Design BPs, Ideal BPs, BPs diff"

print "%.5e, %.5e, %.2e, %.2e"%(na.sum(observed_distr)-(na.sum(spacings[0])*na.sum(spacings[1])), 
												na.sum(expected_distr)-(na.sum(ideal_spacings[0])*na.sum(ideal_spacings[1])),
												abs(na.sum(observed_distr)-na.sum(expected_distr)),
												(na.sum(spacings[0])*na.sum(spacings[1]))-(na.sum(ideal_spacings[0])*na.sum(ideal_spacings[1])))

plot_lengths(observed_distr, 'observed cdf')
plot_lengths(expected_distr, 'ideal cdf')

plt.ylabel('Density')
plt.xlabel('Amplicon length induced by breakpoint')
plt.xlim(plt_xmin,plt_xmax)

plt.legend(loc='lower right')
#plt.legend(loc='upper right')
plt.show()

