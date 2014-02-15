from ambre.utils import reference

r = reference.Reference('/media/T01/data/hg/hg19/human2.rna.fna')
input_fpath = "/home/anand/Projects/ambre/gbm/gbm-rna-screen_cDNA04.txt"

rc_dict = {'forward':False, 'reverse':True}

with open(input_fpath, 'rb') as f:
  for line in f:
    if line.startswith('#'): continue
    tokens = line.strip().split('\t')
    start = int(tokens[1])
    end = start+int(tokens[2])
    rev_complement = tokens[3]=='reverse'
    #print rev_complement
    print ">%s-%d-%d\n%s"%(tokens[0].split('|')[3], start, end, r.get_sequence(tokens[0], start, end, rev_complement=rev_complement))
    #assert r.get_contig(tokens[0])[start:end] ==  r.get_sequence(tokens[0], start, end)
    #print r.get_contig(tokens[0])
