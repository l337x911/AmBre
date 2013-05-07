'''
#  ambre.design.sa_cost.py
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
from ambre.config import CONFIG

import numpy as na
import random
import time
import os
import sys
import cStringIO as StringIO
from multiprocessing import Pool 

def multiplx_filter_func(values, stringency=(-4.0, -8.0, -8.0)):
  s1, s2, s3 = float(values[0]), float(values[1]), float(values[2])
  return s1 < stringency[0] or s2 < stringency[1] or s3 < stringency[2]

class TemperatureScheduleLinear(object):
  def __init__(self, m= -8 / 10.0, b=80000):
    self.m, self.b = m, b
  def get_temperature(self, iteration):
    return max(self.m * iteration + self.b, 1)

class PrimerGraph(object):
  '''
  classdocs
  '''
  def __init__(self, edges_fpath, regions_fpath,
        d=6500, rho=0.1, rho_err=0.1):
    '''
    Constructor: Note that primer positions are selected based on
    concatenation of sequences for primer3 selection and position output.
    '''
    # Expects a PrimerPair Scoring File
    # <primer1_loc>\t<primer2_loc>\t<score>\t<score>\t<score>\n 
    self.edges_fpath = edges_fpath
    
    # Expects a tab delimited 
    # <region left most position>\t<region length>\t<direction (forward|reverse)>\n
    # or can be a string
    self.regions_input = regions_fpath
    
    #self.primer_graph = defaultdict(dict)
    self.map_primers_to_idx = None
    self.primer_graph_mat = None
    self.combined_primers = None
    # prob of not sampling a primer at this position 
    # 
    self.primers_pos_density = None
    self.primers = set()
  
    self.regions = []
    self.region_lengths = []
    self.total_length = 0
    self.region_is_forward = []
    self.number_of_regions = 0
    self.map_primer_idx_to_region = None
    self.primers_by_region = []
    
    self.d = d
    self.rho, self.rho_err = rho, rho_err
    
  def add_artificial_primers(self, primers, dimerization_flag=False):
    primers.sort()
    
    if len(primers) > 1:
      for j in xrange(1, len(primers)):
        for i in xrange(j):
          self.primer_graph[primers[i]][primers[j]] = dimerization_flag
          self.primer_graph[primers[j]][primers[i]] = dimerization_flag
      
    
    for primer in primers:
      for i in self.primers:
        self.primer_graph[i][primer] = dimerization_flag
        self.primer_graph[primer][i] = dimerization_flag
    
    for primer in primers:
      self.primers.add(primer)
  
  def load(self, filter_func):
    
    if os.path.isfile(self.regions_input):
      f = open(self.regions_input, 'rb')
    else:
      f = StringIO.StringIO(self.regions_input)
      
    for l in f:
      tokens = l.strip().split('\t')
      a, b, c = int(tokens[0]), int(tokens[1]), tokens[2]
      self.regions.append((a, a + b))
      self.region_lengths.append(b)
      assert c == "forward" or c == "reverse"
      self.region_is_forward.append(c == "forward")
    f.close()
    
    self.total_length = na.sum(self.region_lengths)
    self.number_of_regions = len(self.regions)

    
    with open(self.edges_fpath, 'rb') as f:
      #Quench Header
      f.next()
      for r in f:
        r = r.strip().split('\t')
        p1, p2 = int(r[0]), int(r[1])
        
        self.primers.add(p1)
        self.primers.add(p2)
        
    self.load_primer_maps()
    
    # generate adj matrix
    self.primer_graph_mat = na.zeros((len(self.primers), len(self.primers)), dtype=na.bool)
    
    with open(self.edges_fpath, 'rb') as f:
      #Quench Header
      f.next()
      for r in f:
        r = r.strip().split('\t')
        i_s, j_s = int(r[0]), int(r[1])
        if i_s > self.combined_primers[-1] or j_s > self.combined_primers[-1] or self.map_primers_to_idx[i_s] < 0 or self.map_primers_to_idx[j_s] < 0:
          continue
        i, j = self.map_primers_to_idx[i_s], self.map_primers_to_idx[j_s]
        a = filter_func(r[2:])
        
        self.primer_graph_mat[i, j] = a
        self.primer_graph_mat[j, i] = a
  
  def load_primer_maps(self):
    self.combined_primers = na.array(list(self.primers))
    for i, (a, b) in enumerate(self.regions):
      mask = na.logical_and(self.combined_primers < b, self.combined_primers >= a)  
      self.primers_by_region.append(self.combined_primers[mask])
      if na.sum(mask) == 0:
        print >> sys.stderr, "Error Region %d-%d has no primers to choose from" % (a, b)
        raise
    
    # Discard all primers that do not fall into the regions of interest
    primers_in_regions = set().union(*[set(i) for i in self.primers_by_region])
    self.primers.intersection_update(primers_in_regions)
    self.combined_primers = na.sort(list(self.primers))
    
    self.map_primers_to_idx = na.ones(self.combined_primers[-1] + 1, dtype=na.uint16) * -1
    self.map_primers_to_idx[self.combined_primers] = na.arange(len(self.primers), dtype=na.uint16)
    self.map_primer_idx_to_region = na.empty(len(self.combined_primers), dtype=na.uint16)
    
    # load density
    for i, (a, b) in enumerate(self.regions):
      mask = na.logical_and(self.combined_primers < b, self.combined_primers >= a)  
      self.map_primer_idx_to_region[mask] = i
    
  def remove_primers(self, hit_set):
    self.primers.intersection_update(hit_set)
          
    x, y = na.meshgrid(na.arange(len(self.primers)), na.arange(len(self.primers)))
    self.primer_graph_mat = self.primer_graph_mat[x, y]

  def get_dimerization_density(self):
    #assert len(self.primers)==self.primer_graph_mat.shape[0]
    #return float(na.sum(self.primer_graph_mat))
    return float(na.sum(self.primer_graph_mat))/ (0.5*len(self.primers)*(1+len(self.primers)))
    

  def check_primer_dimers(self, primer_list):    
    primer_idx = self.map_primers_to_idx[primer_list]
    return na.any(self.primer_graph_mat[na.meshgrid(primer_idx, primer_idx)])

  def adj_cost_func(self, sorted_primers_in_regions, dimer_check_flag=False):
    # Every set contains the "two inner primers_cost"
    # //C(P) = \sum_{(p_i,p_j)\in E} (w_p) + \sum_j wc C(p_j,p_j+1)+ \sum_j w\rho T(p_j,p_j+1)+ (|P|*d/rho - L) 
    # C(P) = \sum_{(p_i,p_j)\in E} (w_p) + \sum_j wc C(p_j,p_j+1)+ \sum_j w\rho T(p_j,p_j+1)
    # Primer Dimers are -Inf
    # Formulation where C(p_j,p_j+1) = max{0, |l_u-l_v|-d}
    # and T(p_j,p_j+1} = max{0,-min{0, |l_u-l_v|-d}-((\rho-1)*d)}.

    
    #l_length = len(sorted_l_primer_list)
    #r_length = len(sorted_r_primer_list)
    #cost = self.l_L+self.r_L-((l_length+r_length)*self.d)/(1+self.rho)
  
    combined_primers = na.concatenate(sorted_primers_in_regions)
    
    if dimer_check_flag and self.check_primer_dimers(combined_primers): 
      return na.inf
  
    cost = 0
    for is_forward, p_list, (a, b) in zip(self.region_is_forward, sorted_primers_in_regions, self.regions):
      primers_cost = na.concatenate(([a], p_list, [b]))
      
      primers_cost = (primers_cost[1:] - primers_cost[:-1]) - self.d
      
      if is_forward:
        primers_cost[0] += self.d
      else:
        primers_cost[-1] += self.d
      
      primers_idx = primers_cost > 0
      
      primers_idx_dense = na.logical_not(primers_idx)
      primers_dense = -primers_cost[primers_idx_dense] - self.d * self.rho
        
      cost += na.sum(primers_cost[primers_idx])

      cost += na.sum(primers_dense[primers_dense > 0])
    
    return cost

  def bashir_cost_func(self, sorted_primers_in_regions, dimer_check_flag=False):
    # Every set contains the "two inner primers"
    # C(P) = \sum_{(p_i,p_j)\in E} (w_p) + \sum_j wc C(p_j,p_j+1) + wrho*rho
    # Formulation where C(p_j,p_j+1) = max{0, |l_u-l_v|-d}
    # Primer Dimers are -Inf
    # wrho is -Inf if out of range of rho.
    
    density = na.sum([self.d * float(len(p_list)) / l for l, p_list in zip(self.region_length, sorted_primers_in_regions)])
    
    if abs(density - self.rho - 1) > self.rho_err:
      return na.inf
    
    combined_primers = na.concatenate(sorted_primers_in_regions)

    if dimer_check_flag and self.check_primer_dimers(combined_primers): 
      return na.inf
    
    cost = 0
    for is_forward, p_list, (a, b) in zip(self.region_is_forward, sorted_primers_in_regions, self.regions):
      primers_cost = na.concatenate(([a], p_list, [b]))
      
      primers_cost = primers_cost[1:] - primers_cost[:-1] - self.d
      if is_forward:
        primers_cost[0] += self.d
      else:
        primers_cost[-1] += self.d
        
      cost += na.sum(primers_cost[primers_cost > 0])

    return cost

  def computeRandom(self):
    # Hardcoded maximum initial size has at most 6 primers per region
    
    independent_set_flag = False
    
    num_primers = []
    for i, l in zip(self.primers_by_region, self.region_lengths):
      num_primers.append(max(min(int(na.floor(((self.rho + 1) * l) / self.d)), min(len(i),6)), 1))
    
    primer_indices_by_region = [range(len(r_init)) for r_init in self.primers_by_region]
    
    tries = 0
    while not independent_set_flag:
      primer_init_idx_by_region = [na.sort(random.sample(indices, primer_count)) for primer_count, indices in zip(num_primers, primer_indices_by_region)]
      
      sorted_primers_init_by_region = [na.sort(primers[init_idx]) for primers, init_idx in zip(self.primers_by_region, primer_init_idx_by_region)]
      
      combined_primers = na.concatenate(sorted_primers_init_by_region)
      independent_set_flag = not self.check_primer_dimers(combined_primers)
      tries += 1
      if tries>10000000:
        print >>sys.stderr, "ERROR: Unexpected number of attempts for finding a random initial solution."
        raise 
    curr_cost = self.adj_cost_func(sorted_primers_init_by_region)
    
    return curr_cost, sorted_primers_init_by_region
    
  def computeSimulatedAnnealing(self, T_schedule, max_iterations=1000000000, max_age=300000, output_fpath=None):

    # Randomly select l and r primers that fulfill the primer density criteria.
    # and form an independent set.
    
    curr_cost, regions_init = self.computeRandom()
    curr_primer_sets = tuple([set(r_init) for r_init in regions_init])
    
    min_primer_sets = tuple([set(r_init) for r_init in regions_init])
    min_cost = curr_cost
    
    neighbor_primer_sets = tuple([set() for i in range(self.number_of_regions)])
    sorted_neighbor_sln_lists = [[]] * self.number_of_regions
    iter_count, min_cost_age = 0, 0
    iter_clock = time.time()
    start_clock = time.time()
    
    if output_fpath is None:
      output_info = sys.stderr
    else:
      output_info = open(output_fpath, 'wb')
    
    
    while iter_count < max_iterations:
      
      # get a neighbor.
      new_neighbor_flag = False
      
      while not new_neighbor_flag:
        u = random.choice(self.combined_primers)
        
        # checks if u belongs to one of the primer sets.
        if sum([u in curr_primer_set for curr_primer_set in curr_primer_sets]) == 0:
          new_neighbor_flag = True

        # gets the indices of neighbors of u
        vs_idx = na.flatnonzero(self.primer_graph_mat[self.map_primers_to_idx[u], :])
        
        # get the neighbors of u
        vs = self.combined_primers[vs_idx]
        
        # randomly select a v if u does not generate a new neighbor
        if len(vs) == 0:
          vs_subset = set()
        else:
          vs_subset = set(vs).intersection(curr_primer_sets[0].union(*curr_primer_sets))

        # compute the costs of neighbor solutions
        for neighbor, curr_primer in zip(neighbor_primer_sets, curr_primer_sets):
          neighbor.clear()
          neighbor.update(curr_primer)
        
        if new_neighbor_flag:
          
          neighbor_primer_sets[self.map_primer_idx_to_region[self.map_primers_to_idx[u]]].add(u)
          for v in vs_subset:
            neighbor_primer_sets[self.map_primer_idx_to_region[self.map_primers_to_idx[v]]].discard(v) 
        else:
          neighbor_primer_sets[self.map_primer_idx_to_region[self.map_primers_to_idx[u]]].discard(u)
        
      sorted_neighbor_sln_lists = [na.sort(list(neighbor_primer_set)) for neighbor_primer_set in neighbor_primer_sets]
      
      neighbor_cost = self.adj_cost_func(sorted_neighbor_sln_lists, False)
      
      # compute probability and swap solutions
      
      T = T_schedule.get_temperature(iter_count)
      
      if neighbor_cost < curr_cost or random.random() < na.exp(float(curr_cost - neighbor_cost) / T):
        curr_cost = neighbor_cost
        
        for curr_primer_set, neighbor_primer_set in zip(curr_primer_sets, neighbor_primer_sets):
          curr_primer_set.clear()
          curr_primer_set.update(neighbor_primer_set)
      
      for neighbor_primer_set in neighbor_primer_sets:
        neighbor_primer_set.clear()
      
      if(curr_cost < min_cost):
        min_cost = curr_cost
        min_cost_age = 0
        
        for min_primer_set, curr_primer_set in zip(min_primer_sets, curr_primer_sets):
          min_primer_set.clear()
          min_primer_set.update(curr_primer_set)

      # Guesses that we are in a limit cycle and need to escape!
      # restart at the best solution
      if min_cost_age > max_age:
        min_cost_age = 0
        for min_primer_set, curr_primer_set in zip(min_primer_sets, curr_primer_sets):
          curr_primer_set.clear()
          curr_primer_set.update(min_primer_set)
      
      iter_count += 1
      min_cost_age += 1
      if(iter_count % 100000 == 0):
        print >> output_info, "Iterations\t%010d\t%.4f\t%d\t%d" % (iter_count, time.time() - iter_clock, min_cost, curr_cost)
        iter_clock = time.time()

    print >> output_info, "TerminationClock:%.2f" % (time.time() - start_clock)
    if not output_fpath is None:
      output_info.close()
    else:
      output_info.flush()
    return min_cost, min_primer_sets, iter_count, time.time() - start_clock   

def run_simulated_annealing(graph, schedule, max_iteration, max_age, output=None):
  random.seed(os.getpid()*time.time())
  return graph.computeSimulatedAnnealing(schedule, max_iteration, max_age, output)

def get_primer_graph(edges_fpath,
                regions_input, d=6500,
                rho=0.1):
  graph = PrimerGraph(edges_fpath, regions_input, d=d, rho=rho)
  graph.load(multiplx_filter_func)
  return graph

def multi_vary_temperature_schedule(graph=None,
                 t_ms=(-1, -(10 ** -1), -(10 ** -2), -(10 ** -3), -(10 ** -4)), t_bs=(10 ** 4, 10 ** 5, 10 ** 6),
                 repeats=5,
                 max_iterations=5000000, max_age=300000,
                 output=None):
                
  p = Pool(processes=int(CONFIG.param['threads']))
  if graph is None:
    graph = get_primer_graph()

  results = []
  for t_m in t_ms:
    for t_b in t_bs:
      for i in xrange(repeats):
          schedule = TemperatureScheduleLinear(m=t_m, b=t_b)
          max_iter = min(max_iterations, 2*t_b+abs(int(round(t_b/t_m))))
          max_iter = max(max_iter, 1000000)
          time.sleep(2)
          if output is None:
            r = p.apply_async(run_simulated_annealing, (graph, schedule, max_iter, max_age))
          else:
            r = p.apply_async(run_simulated_annealing, (graph, schedule, max_iter, max_age, "%s.%03d.%02d.%02d" % (output.name, int(na.log10(-1 * schedule.m)), int(na.log10(schedule.b)), i)))
          results.append((schedule, i, r))
  
  if output is None:
    output = sys.stdout
  
  for schedule, i, result in results:
    cost, primer_sets, iter_count, runtime = result.get()
    print >> output, "%03d,%02d,%02d\t%08d\t%.2f\t%s" % (int(na.log10(-1 * schedule.m)), int(na.log10(schedule.b)), i, cost, runtime,
                    ":".join([",".join(['%d' % s for s in na.sort(list(primer_sets[r_idx]))]) for r_idx in range(len(primer_sets))]))    
  if output is None:
    output.flush()
    
if __name__ == '__main__':
  multi_vary_temperature_schedule(t_ms=(-(10 ** -2), -(10 ** -3)), t_bs=(10 ** 4,), repeats=3)
  multi_vary_temperature_schedule()
