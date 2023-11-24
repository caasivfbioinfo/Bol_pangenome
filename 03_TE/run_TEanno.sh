#Usage: bash run_TEanno.sh genome.fa cds.fa outdir threads 
#All input parameters should be full path; input fasta format: *.fasta; genome seq ID length should be shoter than 12 characters(LTR_harvest required);
reference=$1
cds=$2
outdir=$3
threads=$4

prefix=$(basename ${reference} .fasta)
~/RepeatModeler-2.0.1/BuildDatabase -name ${prefix} -engine ncbi ${reference}
$RepeatModeler -database ${prefix} -engine ncbi -pa 10
#extract unknown sequence
perl process_pickUnknown.pl ${prefix}-families.fa
#align cds
makeblastdb -in ${cds} -dbtype nucl 
blastn -query ${prefix}-families.fa.unknown -db ${cds} -out ukTEs2gene.out -outfmt 6
perl process_filterOverlap.pl ${prefix}-families.fa.unknown ${cds} ukTEs2gene.out ukTE-gene.fa
#self align
makeblastdb -in ukTE-gene.fa -dbtype nucl
blastn -query ukTE-gene.fa -db ukTE-gene.fa -out ukTE-gene2self.out -outfmt 6
perl process_filterOverlap.pl ukTE-gene.fa ukTE-gene.fa ukTEs2gene.out ukTE-gene.self.fa
#combine
cat  ${prefix}-families.fa.known ukTE-gene.self.fa > Repeatmodel_TE.lib
#LTR_finder parallel
ltr_finder_pa -seq ${reference} -finder ltr_finder -t ${threads}

#LTR_harvest
cd $outdir/harvest
time ltr_harvest/gt suffixerator \
	-db ${reference} \
	-indexname ${prefix} \
	-tis -suf -lcp -des -ssp -sds -dna
time ltr_harvest/gt ltrharvest \
	-index ${prefix} \
	-similar 30 -seed 20 -minlenltr 100 -maxlenltr 3500 -mintsd 4 -maxtsd 20 \
	-motif TGCA -vic 60  > ${prefix}.harvest.scn

#LTR_retriver
time ltr_retriver -genome ${prefix}.fasta -inharvest ${prefix}.harvest.scn -infinder ${prefix}.fasta.finder.combine.scn -threads ${threads}


#combine RepeatModeler and LTR_retriver
cat ${prefix}.fasta.LTRlib.fa Repeatmodel_TE.lib > all_LTR.lib
#filter
time python3 process_100bp.py all_LTR.lib over_100.lib low_100.lib
time $cdhit/cd-hit-est -i over_100.lib -o over_100.filter.lib -c 0.8 -aS 0.8 -M 0 -T ${threads} #(80-80-80 rules)
cat low_100.lib over_100.filter.lib > ${prefix}.LTRlib.clust.fa

#RepeatMasker
RepeatMasker -e ncbi -pa 10 -lib ${prefix}.LTRlib.clust.fa -dir ./ -gff ${prefix}.fasta

perl add_class.pl ${prefix}.LTRlib.clust.fa ${prefix}.fasta.out.gff
perl merge_class.pl ${prefix}.fasta.out.gff.class
perl merge_class.pl ${prefix}.fasta.out.gff.total