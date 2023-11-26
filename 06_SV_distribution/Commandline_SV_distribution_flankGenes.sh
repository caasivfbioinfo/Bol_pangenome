perl pickGffPresent_extractPepCds_VG.pl T10.chr.gff T10.chr.fasta
perl process_chromosomeSize_MG.pl T10.chr.fasta
perl split_geneUnit_VG.pl T10.chr.gff.repr T10.chr.fasta.chrSize 5000
perl process_geneFlankGaps.pl T10.chr.gff.repr.5000.geneUnits T10.chr.fasta

perl process_count_elements_Nflank_NN_MG.pl T1.PAV.sort.gff T10.chr.gff.repr.5000.geneUnits T10.chr.fasta genes sv

perl process_stat_TEs_MG.pl T1.PAV.sort.gff.genes.sv T10.chr.gff.repr.5000.geneUnits.gaps T1AllGenes 0 &
perl process_SDoutlier_MG.pl T10.chr.gff.repr.5000.geneUnits.gaps.T1AllGenes0.stat.ALL.sort

perl process_genesElements_fig-6lines_MG.pl T10.chr.gff.repr.5000.geneUnits.gaps.T1AllGenes0.stat.ALL.sort.full T10.chr.gff.repr.5000.geneUnits.gaps.T1AllGenes0.stat.ALL.sort.gexon T10.chr.gff.repr.5000.geneUnits.gaps.T1AllGenes0.stat.ALL.sort.gintron T10.chr.gff.repr.5000.geneUnits.gaps.T1AllGenes0.stat.ALL.sort.full T10.chr.gff.repr.5000.geneUnits.gaps.T1AllGenes0.stat.ALL.sort.gexon T10.chr.gff.repr.5000.geneUnits.gaps.T1AllGenes0.stat.ALL.sort.gintron T1.SV.density.flankGenes
