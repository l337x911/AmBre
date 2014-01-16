import sys
from ambre.analyze.frag_clustering import AssemblyTAlign

talign_fpath = sys.argv[1]
bp_fpath = sys.argv[2]


if __name__ == '__main__':
  clustering = AssemblyTAlign(talign_fpath)
  bps = clustering.condense_fragments_clustering()
  clustering.output_breakpoints(bps,
                                breakpoint_fpath=bp_fpath,
                                min_cluster_size=1)

