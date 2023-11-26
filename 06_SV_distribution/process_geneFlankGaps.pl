#!/usr/bin/perl -w
#Author: Feng Cheng
use warnings;
use strict;
require("/data/chengf/bin/perl_subs.pl");

my $in1 = $ARGV[0];               # gene units
my $in2 = $ARGV[1];               # scaffold_newset.info.sorted0708.chr.fa
my $out1 = $in1.".gaps";

&process($in1,$in2,$out1);

sub process {
  my ($in1,$in2,$out1) = @_;

  my %sca2fas;
  &readFasta($in2,\%sca2fas);

  &output($in1,$out1,\%sca2fas);

}

sub output {
  my ($in1,$out1,$sca2seq) = @_;

  open(my $FR,$in1);
  open(my $FW,">$out1");
  $_=<$FR>;
  chomp;
  print $FW "$_\tLgap\tRgap\n";
  while($_=<$FR>) {
    chomp;
    my @temp = split(/\t/);
    my $len1 = $temp[5] - $temp[3];
    my $seq1 = substr($sca2seq->{$temp[1]},$temp[3],$len1);
    my $len2 = $temp[8] - $temp[6];
    my $seq2 = substr($sca2seq->{$temp[1]},$temp[6],$len2);
    my ($nlen1,$nlen2) = (0,0);
    if($seq1=~/N/) {
      $nlen1 = $seq1=~s/N//g;
    }
    if($seq2=~/N/) {
      $nlen2 = $seq2=~s/N//g;
    }
    
    print $FW "$_\t$nlen1/$len1\t$nlen2/$len2\n";
  }
  close($FR);
  close($FW);
}

