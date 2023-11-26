#only keep SNPs
vcftools --gzvcf merge.vcf.gz --remove-indels --max-missing 0.5 --mac 3 --minQ 30 --minDP 3 --maf 0.05 --recode  --out merge.filtered.SNP

#only keep InDels
vcftools --gzvcf merge.vcf.gz --keep-only-indels  --max-missing 0.5 --mac 3 --minQ 30 --minDP 3 --maf 0.05 --recode  --out merge.filtered.InDel
