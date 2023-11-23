#!/usr/bin/perl -w
#Auther: Feng Cheng
use warnings;
use strict;

my $in = $ARGV[0];
my $out = $in.".merge";

&process($in,$out);

sub process {
  my ($in,$out) = @_;

  &output($in,$out);
}

sub output {
  my ($in,$out) = @_;
  open(my $FR,$in);
  open(my $FW,">$out");
  my ($pchr,$pstart,$pstop,$pfam) = ("NaN",0,0,"NaN");
  my ($c2,$c3,$c6,$c7,$c8,$c9n) = ("RepeatMasker","similarity","","",".","");
  while($_=<$FR>) {
    chomp;
    my @temp = split(/\t/);
    my ($chr,$start,$stop) = ($temp[0],$temp[3],$temp[4]);
    $temp[5] =~ s/\s//g;
    $temp[8] =~ /Target="Motif:(.*)";\d+/;
    my $fam = $1;
    if($chr ne $pchr) {
      if($pchr ne "NaN") {
        print $FW "$pchr\t$c2\t$c3\t$pstart\t$pstop\t$c6\t$c7\t$c8\tTarget=\"Motif:$pfam\";$c9n\n";
      }
      ($pchr,$pstart,$pstop,$pfam) = ($chr,$start,$stop,$fam);
      $temp[8] =~ /;(\d+-\d+)$/;
      ($c6,$c7,$c9n) = ($temp[5],$temp[6],$1);
    }
    elsif($start<=$pstop+10 && $pfam eq $fam) {
      if($pstop < $stop) {
        $pstop = $stop;
        $temp[8] =~ /;(\d+-\d+)$/;
        ($c6,$c7,$c9n) = ($c6."/$temp[5]",$c7."/$temp[6]",$c9n."/$1");
      }
    }
    else {
      print $FW "$pchr\t$c2\t$c3\t$pstart\t$pstop\t$c6\t$c7\t$c8\tTarget=\"Motif:$pfam\";$c9n\n";
      ($pchr,$pstart,$pstop,$pfam) = ($chr,$start,$stop,$fam);
      $temp[8] =~ /;(\d+-\d+)$/;
      ($c6,$c7,$c9n) = ($temp[5],$temp[6],$1);
    }
  }
  close($FR);
  close($FW);
}


