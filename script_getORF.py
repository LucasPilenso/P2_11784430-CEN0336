import sys
import re               #importando as bibliotecas necessárias
from Bio import SeqIO

#Primeiramente vamos definir uma função para encontrar a maior orf
def longest_orf(seqs):
  orfs = []
  if 'ATG' in seqs:
    for start in re.finditer('ATG' , seqs):   #Aqui restringimos somente
      match = seqs[start.start():]          #a região codificadora
      for stop in re.finditer('TAA|TGA|TAG' , match):
        cds = match[:stop.end()]
        if len(cds) % 3 == 0:         #e garantimos que seja lida em trincas
          orfs.append(cds)
          break
        else:
          orfs.append('no orf')
  else:
    orfs.append('no orf')
   #ordenar a lista por tamanho e retornar o maior(ultimo)
  orfs.sort(key=len)
  return orfs[-1]
#As funções abaixo identificam o quadro de leitura pelo resto da divisão da posição de start por 3
def frame(match,code):
  res = re.search(r'{0}'.format(match),code)
  global quadro, position
  if res:
    position = res.span()
    if res.start() % 3 == 0:
      quadro = 1
      return 'frame', quadro, position
    elif res.start() % 3 == 1:
      quadro = 2
      return 'frame', quadro, position
    elif res.start() % 3 == 2:
      quadro = 3
      return 'frame', quadro, position

  else:
    return ''
def rev_frame(rmatch,rcode):    #No caso das orfs reversas:
  resr = re.search(r'{0}'.format(rmatch),rcode)
  global rquadro, rposition
  if resr:
    rposition = resr.span()
    if resr.start() % 3 == 0:
      rquadro = 4
      return 'frame', rquadro, rposition
    elif resr.start() % 3 == 1:
      rquadro = 5
      return 'frame', rquadro, rposition
    elif resr.start() % 3 == 2:
      rquadro = 6
      return 'frame', rquadro, rposition

  else:
    return ''

with open(sys.argv[1], 'r') as fasta_file,  open('ORF.fna', 'w') as orf_file,  open('ORF.faa', 'w') as ptn_file:
  id = []
  forw_rev= []
  final_orf=[]
  for info in SeqIO.parse(fasta_file, "fasta"):
    id.append(info.id)
    forw_rev.append(longest_orf(str(info.seq)))
    forw_rev.append(longest_orf(str(info.seq.reverse_complement())))
  for i in range(0, len(forw_rev), 2):
    sublist = forw_rev[i:i+2]
    index, value = max(enumerate(sublist), key=lambda x: x[1])
  # adiciona o maior elemento da sublista na lista de maiores elementos
    final_orf.append(value)
  for i, orf in enumerate(final_orf):

    if str(orf) in str(info.seq):
      orf_file.write('{} {} \n {}'.format(id[i],frame(longest_orf(str(orf)),str(orf)) , longest_orf(str(orf))))
      ptn_file.write('{} {} \n {}'.format(id[i],frame(longest_orf(str(orf)),str(orf)) , longest_orf(str(orf.translate())))
    else:
      orf_file.write('{} {} \n {}'.format(id[i],rev_frame(longest_orf(str(orf)),str(orf)) , longest_orf(str(orf)))
      ptn_file.write('{} {} \n {}'.format(id[i],rev_frame(longest_orf(str(orf)),str(orf)) , longest_orf(str(orf.translate()))
