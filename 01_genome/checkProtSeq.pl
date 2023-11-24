#!/usr/bin/perl
=head2 SYNOPSIS

	Description: 检查起始和终止，并确定其中内部是否含有终止密码子
	protein 序列的提取来自于 cufflinks中的gffread gene_prediction.gff -g genome.fa -y protein.fa
	因为该程序的output protein fa中终止密码子以"." 表示
	
	Usage:
		./this_script.pl protein.fa
		
	Output:
		protein.pass  	通过检查的蛋白
		protein.fail	没有通过的蛋白
=head2 AUTHOR

	Kang Zhang  Email:  demut@foxmail.com

=cut

open FA,"<$ARGV[0]" or die;
my @name=split /\./,$ARGV[0];
my $prefix=$name[0];
open PASS,">$prefix.pass";
open FAIL,">$prefix.fail";
while(<FA>){
  chomp;
  if(/^>(\S+)/){
    $prot=$1;
  }else{
    $length{$prot}.=$_;
 }
}

foreach my $prot (sort keys %length){
    my $seq=$length{$prot};
    if($seq=~/^M.*\.$/){
      if($seq=~/\..+\./){
        print FAIL "$prot\tFail\tInner Stop Codon!\n";
      }else{
        print PASS "$prot\n";
      }
    }else{
      print FAIL "$prot\tFail\tUnaccepted Start/Stop Codon!\n";
    }
 }

