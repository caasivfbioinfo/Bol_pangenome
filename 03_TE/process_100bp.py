import sys
from Bio import SeqIO
over_100=open(sys.argv[2],"w")
low_100=open(sys.argv[3],"w")

for rec in SeqIO.parse(sys.argv[1],"fasta"):
    ids = rec.id
    seqs = str(rec.seq)
    lenth = len(seqs)
    if lenth >= 100:
        over_100.write(">"+ids+"\n"+seqs+"\n")
    else:
        low_100.write(">"+ids+"\n"+seqs+"\n")

over_100.close()
low_100.close()