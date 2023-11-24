vg construct -r refgenome -v merge.vcf.gz >../merge.vg
vg index -x merge.xg -g merge.gcsa -k 16 -t 160 -b /tmp -Z 10000 merge.vg

for i in $(cat list.txt)
do
vg map  -d merge -f ${i}_1.fastq.gz -f ${i}_2.fastq.gz -t 20 >gam/${i}.aln.gam
vg pack -x merge.xg -g gam/${i}.aln.gam -Q 5 -t 20 -o gam/${i}.aln.pack
vg call merge.xg -k gam/${i}.aln.pack -t 20 -a -s ${i} >vcf/${i}.vcf
bgzip vcf/${i}.vcf
tabix vcf/${i}.vcf.gz
done
