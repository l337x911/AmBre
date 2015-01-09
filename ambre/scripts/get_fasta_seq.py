from ambre.utils import reference
import sys

def aml_eto_fusion():
  ref = reference.Reference('/home/adp002/data/hg/hg19_full.fa')

  chr8 = ref.get_contig('chr8')
  chr21 = ref.get_contig('chr21')

  k = int(sys.argv[1])

  try:
    seq = sys.argv[2]
    seq_rc = seq.translate(reference.BASE_COMPLEMENT)[::-1]
  except:
    seq = ''

  input_coord = [(('chr8', 93048551, 93048569), ('chr21', 36230187, 36230207)),
  (('chr21',36229647,36229667), ('chr8', 93066669, 93066689)), (('chr8', 93040017, 93040039), ('chr21', 36208889, 36208909))]

  output_str = []

  for x,y in input_coord:
    x_c, x_a, x_b = x
    y_c, y_a, y_b = y
    output_str.append((x_c,y_c, eval(y_c)[y_b-1+k:y_a-1:-1].translate(reference.BASE_COMPLEMENT), eval(x_c)[x_b-1:x_a-k-1:-1].translate(reference.BASE_COMPLEMENT)))

  for x,y, x_s, y_s in output_str:
    print x,y,len(y_s), len(x_s), (x_s+ y_s).translate(reference.BASE_COMPLEMENT)[::-1]
    if not seq == '':
      print x_s, y_s
      print x_s.upper().find(seq), x_s.upper().find(seq_rc), y_s.upper().find(seq), y_s.upper().find(seq_rc)

def gbm_fusion():
  ref = reference.Reference('/media/ET01/data/hg/hg19/full.fa')

  chr7 = ref.get_contig('chr7')
  w = 4000

  breakpoints = [(116325515,121602420),(116325490,121613840)]
  # expect sequence [..., b] [a, ...]
  for a,b in breakpoints:
    print ">chr7:{0:d}-{1:d};{2:d}".format(a,b,w)
    print "{0}{1}".format(chr7[b-1-w:b-1], chr7[a-1:a-1+w])

if __name__ == '__main__':
  #aml_eto_fusion()
  gbm_fusion()

