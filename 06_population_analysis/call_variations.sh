#call variants by GATK.
#usage: bash run.sh ref_genome data_directory thread_number
#all input should be full path.
#The name of ref_genome should be *.fasta.
#The name of reseq data should be *.R1.fastq.gz *.R2.fastq.gz

genome=$1
data=$2
threads=$3

##referfence build index
mkdir reference && cd reference
ln -s ${genome} ./
bwa index *.fasta
samtools faidx *.fasta

refpre=$(basename *.fasta .fasta)
gatk CreateSequenceDictionary -R ${refpre}.fasta -O ${refpre}.dict
cd ../

##alignment
mkdir bam
mkdir gvcf
mkdir vcf
mkdir Trim && cd Trim
ln -s ${data}/* ./
ls *.R1.fastq.gz | cut -d "." -f1 >data.txt

for i in $(cat ./data.txt):
do
java -jar /data/pub/yuhl/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${threads} -phred33 ${i}.R1.fastq.gz ${i}.R2.fastq.gz ${i}.R1.paired.fastq.gz ${i}.R1.unpaired.fastq.gz ${i}.R2.paired.fastq.gz ${i}.R2.unpaired.fastq.gz ILLUMINACLIP:$adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

bwa mem -R '@RG\tID:foo\tLB:library1\tPU:runbarcode\tPL:Illumina\tDS:resequencing\tCN:MajorBio\tSM:'${i} -t ${threads} -M ../reference/${refpre}.fasta ${i}.R1.paired.fastq.gz  ${i}.R2.paired.fastq.gz  >../bam/${i}.sam

samtools view -bS -h -@ ${threads} ../bam/${i}.sam >../bam/${i}.bam

samtools sort -@ ${threads} ../bam/${i}.bam > ../bam/${i}.sorted.bam

rm ../bam/${i}.sam ../bam/${i}.bam

gatk MarkDuplicates -I ../bam/${i}.sorted.bam -O ../bam/${i}.sorted.markdup.bam -M ../bam/${i}.sorted.markdup_metrics.txt

samtools index -m 15 ../bam/${i}.sorted.markdup.bam

gatk HaplotypeCaller --emit-ref-confidence GVCF   -R ../reference/${refpre}.fasta   -I ../bam/${i}.sorted.markdup.bam -O ../gvcf/${i}.g.vcf
done

gatk CombineGVCFs -R ../reference/${refpre}.fasta $(for i in $(ls ../gvcf/*.vcf);do echo "-V $i";done) -O ../vcf/combine.g.vcf

gatk GenotypeGVCFs -R ../reference/${refpre}.fasta -V ../vcf/combine.g.vcf  -O ../vcf/combine.vcf --allow-old-rms-mapping-quality-annotation-data

bgzip ../vcf/combine.vcf

tabix -p vcf ../vcf/combine.vcf