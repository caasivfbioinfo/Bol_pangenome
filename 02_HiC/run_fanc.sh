#prepare genome
fanc fragments -c C01,C02,C03,C04,C05,C06,C07,C08,C09 T4.chr.fasta HindIII T4.chr.fasta.frag.bed
bwa index T4.chr.fasta

#prepare Hi-C data
ln -s /data/lix/data/Cabbage/reseq_data/hic/T4_merge_good_1.fq.gz .
ln -s /data/lix/data/Cabbage/reseq_data/hic/T4_merge_good_2.fq.gz .

#run fanc
fanc auto T4_merge_good_1.fq.gz T4_merge_good_2.fq.gz ./output -i index/T4.chr.fasta -g T4.chr.fasta.frag.bed -n T4_hic -t 60
fanc insulation output/hic/binned/T4_hic_100kb.hic T4.chr.hic.100kb.insulation

#show TAD
fancplot --width 6 -o T4.chr.hic.100kb.insulation.C09.20mb-50mb.line.png C09:20mb-50mb -p triangular output/hic/binned/T4_hic_100kb.hic -m 15000000 -vmin 0 -vmax 0.02 -p line T4.chr.hic.100kb.insulation_500kb.bed T4.chr.hic.100kb.insulation_1mb.bed -l "500kb" "1mb"
