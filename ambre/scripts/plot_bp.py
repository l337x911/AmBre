from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib import ticker as tkr
from collections import defaultdict
import numpy as na


def comma_func(x,pos):
  return '{:20,d}'.format(int(x)/1000)

rounding = 10000
comma_format = tkr.FuncFormatter(comma_func)

QUICK_COUNT = {'CEM':18821, 'T98G':3021, 'Detroit562':1662, 'MCF7':1087, 'A549_1':186+201, 'A549_2':217+175}

fix = [('A549_2', 21907127, 22123317), ('A549_1',21832461,21907452)]

def resize(s):
  return s/8.

def map_amplicon_size_to_float(abcd):
  return (abcd-2000)/(13000.0-2000.0)

def roundtoNice(x):
  return round(x,-3) if x>1000 else round(x,-2)

class Breakpoints(object):
  def __init__(self):
    pass
  def parse_breakpoints(self, fpath):
    self.samples_for = defaultdict(list)
    self.samples_rev = defaultdict(list)
    with open(fpath, 'rb') as f:
      for s in f:
        if s[0] == '#':
          continue
        
        tokens = s.strip().split('\t')
        sample_name = tokens[0]
        pos = map(int, tokens[2:])
        if pos[4] > 0:
          self.samples_for[sample_name].append(pos[4])
        else:
          self.samples_for[sample_name].append(pos[0])
          self.samples_for[sample_name].append(pos[1])
        if pos[5] > 0:
          self.samples_rev[sample_name].append(pos[5])
        else:
          self.samples_rev[sample_name].append(pos[2])
          self.samples_rev[sample_name].append(pos[3])
  def _plot_interval(self, p):
    points = sorted(p)
    plot_start = points[0] - rounding
    plot_end = points[-1] + rounding
    
    plot_start = na.floor(plot_start / rounding) * rounding
    plot_end = na.ceil(plot_end / rounding) * rounding
    return plot_start, plot_end

  def plot(self):
    fig = plt.figure()
    gs = gridspec.GridSpec(1,2,width_ratios=[4,1])
    ax = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    sample_names = set(self.samples_rev.keys()).union(set(self.samples_for.keys()))
    sample_names.discard('A549')
    

    for s,a,b in fix:
      sample_names.add(s)
      self.samples_for[s].append(a)
      self.samples_rev[s].append(b)

    sample_names = sorted(sample_names)

    l = [self.samples_for[s][0] for s in sample_names]
    r = [self.samples_rev[s][0] for s in sample_names]
    
    counts = map(resize, [QUICK_COUNT[s] for s in sample_names])
    
    ax.scatter(r,l, s=counts, c=map(map_amplicon_size_to_float, counts) , linewidths=0)

    ax.xaxis.set_major_formatter(comma_format)
    ax.yaxis.set_major_formatter(comma_format)

    xlabels = ax.get_xticklabels() 
    for label in xlabels: 
      label.set_rotation(30) 
    ylabels = ax.get_yticklabels() 
    for label in ylabels: 
      label.set_rotation(30) 

   
    plot_start, plot_end = self._plot_interval(l)
    ax.set_yticks(na.arange(plot_start, plot_end + 2, 20000))
    ax.set_ylabel('Left breakpoint (Kbp)')
    ax.set_ylim(plot_start, plot_end)

    plot_start, plot_end = self._plot_interval(r)
    ax.set_xticks(na.arange(plot_start, plot_end + 2, 20000))
    ax.set_xlabel('Right breakpoint (Kbp)')
    ax.set_xlim(plot_start, plot_end)
    text_orient = 'right'

    #number_of_scales = 6
    number_of_scales = 100

    x = na.linspace(100,na.ceil(max(QUICK_COUNT.values()) / rounding+1) * rounding , number_of_scales)
    #y =  na.power(1.2,na.arange(number_of_scales))

    y =  na.power(1.001,na.arange(number_of_scales))
    print x,y
    ax2.scatter([0]*len(x),y, s=map(resize,x), linewidths=0)
    for text_depth,text_pos in zip(x,y):
      ax2.text(-1,text_pos, "%5d"%roundtoNice(text_depth), horizontalalignment=text_orient, verticalalignment='center')

    
    x = na.linspace(2001,13001, number_of_scales)
    
    
    blocks = na.zeros((y.size+1,))
    blocks[1:-1]= (y[1:]-y[:-1])/2
    blocks[0] = y[0]-blocks[1]
    blocks[-1] = y[-1]+blocks[-2]
    blocks[1:-1] += y[:-1]
    cmap = plt.get_cmap('jet')
    print blocks
    
   
    for bs, bm, be, v in zip(blocks[:-1], y, blocks[1:], x):
    
      ax2.plot((7,7),(bs,be),c=cmap(map_amplicon_size_to_float(v)), lw=16)
    for bs, bm, be, v in zip(blocks[:-1], y, blocks[1:], x)[::10]:
      ax2.text(6,bm, "%5d"%roundtoNice(v), horizontalalignment=text_orient, verticalalignment='center') 

    ax2.set_xlim(-3,8)

    ax2.axis('off')
 

    plt.show()
    

bp = Breakpoints()
bp.parse_breakpoints('/home/anand//Projects/PAMP/A2_breakpoints_hg19_empirical.txt')

bp.plot()

