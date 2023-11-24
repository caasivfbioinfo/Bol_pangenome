#!/usr/bin/perl -w
#Author: Feng Cheng
use warnings;
use strict;
#require "/home/chengf/bin/perl_subs.pl";

my $in1 = $ARGV[0];  #-query--
my $in2 = $ARGV[1];  #-subject--
my $blast12 = $ARGV[2];
my $out = $ARGV[3];
my ($cover,$ident,$len) = (0.8,0.8,80);
#gene 0 0.8 80
#RNA  0 0.8 50
#self 0.8 0.8 80

&process($in1,$in2,$blast12,$cover,$ident,$len);

sub process {
  my ($in1,$in2,$blast12,$cover,$ident,$len) = @_;

  my %q2seq;
  &readFasta($in1,\%q2seq);

  my %r2seq;
  &readFasta($in2,\%r2seq);

  my %q2rm;
  &readBlast($blast12,$cover,$ident,$len,\%q2seq,\%q2rm);

  &output($out,\%q2seq,\%q2rm);
}

sub output {
  my ($out,$q2seq,$q2rm) = @_;
  open(my $FW,">$out");
  foreach my $q (keys %$q2seq) {
    if(!exists($q2rm->{$q})) {
      print $FW ">$q\n$q2seq->{$q}\n";
    }
  }
  close($FW);
}

sub readBlast {
  my ($blast12,$cover,$ident,$len,$q2seq,$q2rm) = @_;

  my %rm2rec;
  open(my $FR,$blast12);
  while($_=<$FR>) {
    chomp;
    my @temp = split(/\t/);
    if($temp[0] ne $temp[1] && exists($q2seq->{$temp[0]}) && !exists($rm2rec{$temp[0]}) && !exists($rm2rec{$temp[1]})) {
      my $qlen = length($q2seq->{$temp[0]});
      #my $slen = $temp[7]-$temp[6]+1;
      my $slen = $temp[3];
      my $scover = $slen/$qlen;            #query length
      my $sident = $temp[2];
      if($scover>=$cover && $sident>=$ident) {
      #if(($scover>=$cover || $slen>=$len) && $sident>=$ident) {
        $q2rm->{$temp[0]} = "Y";
        print "$scover\t$slen/$qlen\t$sident\t||\t$_\n";
        $rm2rec{$temp[0]} = "Y";
      }
    }
  }
  close($FR);
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
