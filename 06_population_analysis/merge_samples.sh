#
gatk CombineGVCFs -R ../reference/genome.fasta $(for i in $(ls ../gvcf/*.vcf);do echo "-V $i";done) -O ../vcf/combine.g.vcf

gatk GenotypeGVCFs -R ../reference/genome.fasta -V ../vcf/combine.g.vcf  -O ../vcf/combine.vcf --allow-old-rms-mapping-quality-annotation-data

bgzip ../vcf/combine.vcf

tabix -p vcf ../vcf/combine.vcf