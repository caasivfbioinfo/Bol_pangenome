#!/usr/bin/perl -w
#Author: Feng Cheng
use warnings;
use strict;
#require "/home/chengf/bin/perl_subs.pl";

my $in1 = $ARGV[0];
my $in2 = $ARGV[1];
my $out1 = $in2.".class";
my $out2 = $in2.".total";

&process($in1,$in2,$out1,$out2);

sub process {
  my ($in1,$in2,$out1,$out2) = @_;

  my %fam2class;
  &readFam($in1,\%fam2class);

  &output($in2,$out1,$out2,\%fam2class);

}

sub output {
  my ($in,$out1,$out2,$fam2class) = @_;
  open(my $FR,$in);
  open(my $FW1,">$out1");
  open(my $FW2,">$out2");
  while($_=<$FR>) {
    if(!/^\#/) {
      chomp;
      my @temp = split(/\t/);
      $temp[8] =~ /\"Motif\:(.+)\"\s(-*\d+)\s(\d+)/;
      my ($fam,$s1,$s2) = ($1,$2,$3);
      s/Target\s\"Motif\:(.+)\"\s(-*\d+)\s(\d+)/Target=\"Motif\:$fam\";$s1-$s2/;
      if(exists($fam2class->{$fam})) {
        s/$fam/$fam\#$fam2class->{$fam}/;
        print $FW1 "$_\n";
      }
      #else {
        print $FW2 "$_\n";
      #}
    }
  }

}

sub readFam {
  my ($in,$fam2class) = @_;
  open(my $FR,$in);
  while($_=<$FR>) {
    if(/^>([^#]+)\#(.+)/) {
      $fam2class->{$1} = $2;
    }
  }
  close($FR);
}


