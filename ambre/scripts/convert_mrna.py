from ambre.utils import reference
import sys
import pandas as pd
import numpy as na
import bisect

refseq_fpath = '/media/T01/data/hg/hg19/mrna-refseqali.txt'
mrna_fasta_fpath = '/media/T01/data/hg/hg19/human2.rna.fna'

def has_qnames(df):
  # EGFR, PDGFRA, SEPT14, PTPRZ1, MET, FGFR3, TACC3, LANCL2, RP11-745C15.2
  qname = pd.DataFrame({'qName':["NM_005228", "NM_006206", "NM_207366",
"NM_002851", "NM_001127500", "NM_000142", "NM_006342", "NM_018697", "NR_110040"]})
  valid_length = len(qname)
  valid = pd.merge(qname, df, how='inner')
  join_length = len(valid)
  #print valid.values
  return join_length==valid_length, valid_length-join_length

mrna_ref = reference.Reference(mrna_fasta_fpath)
refids,contig_names,csizes = zip(*[(c.split('|')[3].split('.')[0],c,mrna_ref.get_contig_length(c)) for c in mrna_ref.get_contig_names()])
refgene = pd.DataFrame({'qName':refids,'qSize':csizes, 'contigs':contig_names})

def test_lengths():
  for c in contig_names:
    try:
      assert mrna_ref.get_contig_length(c)==len(mrna_ref.get_contig(c))
    except AssertionError:
      print mrna_ref.get_contig_length(c)-len(mrna_ref.get_contig(c))

refseq = pd.read_csv(refseq_fpath, header=0, sep='\t')

print "RefSeq Valid?", has_qnames(refseq)
print "RefGene Valid?", has_qnames(refgene)

def get_contig_name(qname):
  return refgene[refgene.qName==qname].contigs.values[0]

#print refseq['qName'].describe()
#print refseq['qName']
#print refids['qName']
#match_ref_to_fq = pd.merge(refgene, refseq, on='qName')
#eq_size_rna = match_ref_to_fq['qSize.x']-match_ref_to_fq['qSize.y']
#print eq_size_rna[eq_size_rna!=0]
#test_lengths()

def comma_split_int(row):
  return na.array(map(int, row.split(',')[:-1]))

def get_exons(refid, print_flag=False):
  g = refseq[refseq.qName==refid]
  exons = g.ix[:,('blockSizes','qStarts','tStarts',)].applymap(comma_split_int)

  qstarts = exons['qStarts'].values[0]
  tstarts = exons['tStarts'].values[0]
  blocks = exons['blockSizes'].values[0]

  if g['strand']=='+':
    tends = tstarts+blocks
  else:
    # if the strand is minus
    # then the alignment should follow 5-3 along the rna
    # tstarts should be tends reversed

    qstarts = na.concatenate(([0,],abs(qstarts[::-1]-g.ix[:,'qSize'].values)[:-1]))
    tstarts = (tstarts+blocks)[::-1]
    blocks = blocks[::-1]
    tends = tstarts-blocks

  qends = qstarts+blocks
  
  if print_flag:
    print qstarts
    print qends
    print tstarts
    print tends
    print blocks
    print blocks[:-1]+blocks[1:]

  return qstarts, qends, tstarts, tends, blocks

#print pd.merge(refids, refseq, on=['qName', 'qSize'])
def get_closest_exon(refid,pos, start5=False):
  g = refseq[refseq.qName==refid]
  exons = g.ix[:,('blockSizes','qStarts','tStarts',)].applymap(comma_split_int)
  #print exons['qStarts'].values

  qstarts = exons['qStarts'].values[0]
  tstarts = exons['tStarts'].values[0]
  blocks = exons['blockSizes'].values[0]

  tends = tstarts+blocks
  if start5:
    idx = na.argmin(na.abs(tends-pos))
    v = tends[idx]
    qv = (qstarts+blocks)[idx]
  else:
    idx = na.argmin(na.abs(tstarts-pos))
    v = tstarts[idx]
    qv = (qstarts)[idx]
 
  return idx+1, abs(v-pos), v, qv


def ambredel_input(refid, end5, start3, w): 
  g = refseq[refseq.qName==refid]
  qstarts, qends, tstarts, tends, blocks = get_exons(refid)

  end5_idx = end5-2
  start3_idx = start3
  #print end5_idx, start3_idx, tends[end5_idx], tstarts[start3_idx]
  print "#", g['tName'].values, g['qSize'].values, qends[end5_idx], qstarts[start3_idx], tends[end5_idx], tstarts[start3_idx], qstarts[start3_idx]-qends[end5_idx]
  qs1, qw1 = max(qends[end5_idx]-w,0), min(w,qends[end5_idx])
  print "%s\t%d\t%d\tforward"%(get_contig_name(refid), qs1,qw1)
  print "%s\t%d\t%d\treverse"%(get_contig_name(refid), qstarts[start3_idx],min(w, g['qSize']-qstarts[start3_idx]))

def ambrefus_input(refid,exon_c, w, start5=False):
  g = refseq[refseq.qName==refid]
  qstarts, qends, tstarts, tends, blocks = get_exons(refid)
  #exon_c, diff, v, qv  = get_closest_exon(refid,pos,start5)
  idx = exon_c-1

  if start5:
    orient = "forward"
    qv = qends[idx]
    pos = max(qv-w,0)
    qw = min(w,qv)
  else:
    orient = "reverse"
    qv = qstarts[idx]
    pos = qv
    qw = min(w,g['qSize']-qv)
    
  print "#", g['qName'].values, g['qSize'].values, qstarts[idx], qends[idx], tstarts[idx], tends[idx]

  print "%s\t%d\t%d\t%s"%(get_contig_name(refid), pos,qw, orient)

W = 400
# EGFR, PDGFRA, SEPT14, PTPRZ1, MET, FGFR3, TACC3, LANCL2, RP11-745C15.2

#print get_closest_exon("NM_005228", "chr7", 55268106, start5=True)
#print get_closest_exon("NM_006206", "chr7", 55268106, start5=True)
#EGFR vIII
#print_exons('NM_005228', (2,7))
ambredel_input('NM_005228', 2,7,W)
#EGFR vII
#print_exons('NM_005228', (14,15))
ambredel_input('NM_005228', 14,15,W)
#PDGFR
#print_exons('NM_006206', (8,9))
ambredel_input('NM_006206', 8,9,W)
#SEPT14 3 exon 3, 7, 10
ambrefus_input("NM_207366", 3, W, start5=False)
ambrefus_input("NM_207366", 7, W, start5=False)
ambrefus_input("NM_207366", 10, W, start5=False)
#print_exons("NM_207366")
# PTPRZ1 exon 1, 2
#print_exons("NM_002851")
ambrefus_input("NM_002851", 1, W, start5=True)
ambrefus_input("NM_002851", 2, W, start5=True)
# MET exon 2 
#print_exons("NM_001127500")
ambrefus_input("NM_001127500", 2, W, start5=False)
# FGFR3 5 exon 17
ambrefus_input("NM_000142", 17, W, start5=True)
#print_exons("NM_000142")
# TACC3 3 exon 8, 10, 11
#print_exons("NM_006342")
ambrefus_input("NM_006342", 8, W, start5=False)
ambrefus_input("NM_006342", 10, W, start5=False)
ambrefus_input("NM_006342", 11, W, start5=False)
# LANCL2 5 exon 1, 6 
#print_exons("NM_018697")
ambrefus_input("NM_018697", 1, W, start5=True)
ambrefus_input("NM_018697", 6, W, start5=True)
# RP11-745C15.2
#print_exons("NR_110040")

