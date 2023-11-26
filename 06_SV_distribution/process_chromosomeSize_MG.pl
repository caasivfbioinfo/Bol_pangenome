#!/usr/bin/perl -w
#Author: Feng Cheng
use warnings;
use strict;
#require "/home/chengf/bin/perl_subs.pl";

my $in = $ARGV[0];
my $out = $in.".chrSize";
$out =~ s/^.*\///;

&process($in,$out);

sub process {
  my ($in,$out) = @_;

  my %sca2seq;
  &readFasta($in,\%sca2seq);

  &output($out,\%sca2seq);

  system("sort -k2,2n $out > $out.sorted");
  system("mv $out.sorted $out");
}

sub output {
  my ($out,$sca2seq) = @_;
  open(my $FW,">$out");
  foreach my $sca (sort {$a cmp $b} keys %$sca2seq) {
    my $len = length($sca2seq->{$sca});
    print $FW "$sca\t$len\n";
  }
  close($FW);
}

sub readFasta {
  my ($in,$id2seq) = @_;
  open(my $SFR,$in);

  my $id;
  while($_=<$SFR>) {
    if(/^>([^\s^\n]+)\s*\n*/) {
      $id = $1;
      $id2seq->{$id} = "";
    }
    else {
      chomp;
      $id2seq->{$id} .= $_;
    }
  }
  close($SFR);
}

