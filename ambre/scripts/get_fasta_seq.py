from ambre.utils import reference
import sys

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


