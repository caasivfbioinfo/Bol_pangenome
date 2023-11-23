#!/usr/bin/perl -w
#Author: Feng Cheng
use warnings;
use strict;
#require "/home/chengf/bin/perl_subs.pl";

my $in = $ARGV[0];
my $out1 = $in.".known";
my $out2 = $in.".unknown";

&process($in,$out1,$out2);

sub process {
  my ($in,$out1,$out2) = @_;

  my %sca2seq;
  &readFasta($in,\%sca2seq);

  &output($out1,$out2,\%sca2seq);

}

sub output {
  my ($out1,$out2,$sca2seq) = @_;

  open(my $FW1,">$out1");
  open(my $FW2,">$out2");
  foreach my $sca (sort {$a cmp $b} keys %$sca2seq) {
    if($sca =~ /Unknown/) {
      print $FW2 ">$sca\n$sca2seq->{$sca}\n";
    }
    else {
      print $FW1 ">$sca\n$sca2seq->{$sca}\n";
    }
  }
  close($FW1);
  close($FW2);
}

sub readFasta {
  my ($in,$id2seq) = @_;
  open(my $SFR,$in);

  my $id;
  while($_=<$SFR>) {
    if(/^>([^\s^\n]+)\s*\n*/) {
      $id = $1;
      if(exists($id2seq->{$id})) {
        print "$id\n";
      }
      $id2seq->{$id} = "";
    }
    else {
      chomp;
      $id2seq->{$id} .= $_;
    }
  }
  close($SFR);
}


