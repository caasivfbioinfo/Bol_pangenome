#!/usr/bin/perl -w
#Author: Feng Cheng
use strict;
use warnings;

my $in1 = $ARGV[0];
my $out1 = $in1.".full";
my $out2 = $in1.".gexon";
my $out3 = $in1.".gintron";

&process($in1,$out1,$out2,$out3);

sub process {
  my ($in1,$out1,$out2,$out3) = @_;

  my (%pos2ave,%pos2std);
  &readData($in1,\%pos2ave,\%pos2std);

  my ($med,$uquar,$dquar) = &sortdata(\%pos2ave,\%pos2std);

  my $uwhisk = "NA";
  if($uquar ne "NA" && $dquar ne "NA" && $uquar ne "NA") {
    $uwhisk = ($uquar-$dquar)*1.5 + $uquar;
    print "$dquar\t$med\t$uquar\t$uwhisk\n";
  }

  &output($out1,$out2,$out3,$uwhisk,\%pos2ave,\%pos2std);

}

sub output {
  my ($out1,$out2,$out3,$uwhisk,$pos2ave,$pos2std) = @_;

  my %pos2rec;
  open(my $FW1,">$out1");
  open(my $FW2,">$out2");
  open(my $FW3,">$out3");
  foreach my $pos (sort {$a<=>$b} keys %$pos2ave) {
    my $tpos = int $pos;
    if(!exists($pos2rec{$tpos})) {
      $pos2rec{$tpos} = "Y";

      if(exists($pos2ave->{$pos})) {
        #if($pos2std->{$pos}<=$uwhisk) {                      #it depends**********************************************
        if($pos2ave->{$pos} ne "NA" && $pos2std->{$pos} ne "NA") {
          my ($left,$right) = ($pos2ave->{$pos}-$pos2std->{$pos},$pos2ave->{$pos}+$pos2std->{$pos});
          print $FW1 "ALL\t$pos\t$pos2ave->{$pos}\t$left\t$right\n";
        }
        #}
        #else {
          #print "$pos\t$pos2ave->{$pos}\t-----\t$pos2std->{$pos}\t$uwhisk\n";
        #}
      }
      $tpos += 0.2;
      if(exists($pos2std->{$tpos})) {
        #if($pos2std->{$tpos}<=$uwhisk) {
        if($pos2ave->{$pos} ne "NA" && $pos2std->{$pos} ne "NA") {
          my ($left,$right) = ($pos2ave->{$tpos}-$pos2std->{$tpos},$pos2ave->{$tpos}+$pos2std->{$tpos});
          print $FW2 "ALL\t$pos\t$pos2ave->{$tpos}\t$left\t$right\n";
        }
        #}
        #else {
          #print $FW2 "\tNA\tNA\tNA\n";
        #}
      }
      $tpos += 0.2;
      if(exists($pos2std->{$tpos})) {
        #if($pos2std->{$tpos}<=$uwhisk) {
        if($pos2ave->{$pos} ne "NA" && $pos2std->{$pos} ne "NA") {
          my ($left,$right) = ($pos2ave->{$tpos}-$pos2std->{$tpos},$pos2ave->{$tpos}+$pos2std->{$tpos});
          print $FW3 "ALL\t$pos\t$pos2ave->{$tpos}\t$left\t$right\n";
        }
        #}
        #else {
          #print $FW2 "\tNA\tNA\tNA";
        #}
      }
    }
  }
  close($FW1);
  close($FW2);
  close($FW3);
}

sub sortdata {
  my ($pos2ave,$pos2std) = @_;

  my @std = values %$pos2std;
  my $num = @std;
  my @stds = sort {$a<=>$b} @std;
  my ($med,$uquar,$dquar);
  my ($min1,$min2,$min3) = (1000,1000,1000);
  for(my $i=0;$i<@stds;$i+=1) {
    #print "$stds[$i]\n";
    if($min1>abs(($i+1)/$num-0.5)) {
      $min1 = abs(($i+1)/$num-0.5);
      $med = $stds[$i];
    }
    if($min2>abs(($i+1)/$num-0.75)) {
      $min2 = abs(($i+1)/$num-0.75);
      $uquar = $stds[$i];
    }
    if($min3>abs(($i+1)/$num-0.25)) {
      $min3 = abs(($i+1)/$num-0.25);
      $dquar = $stds[$i];
    }
  }
  return ($med,$uquar,$dquar);
}

sub readData {
  my ($in,$pos2ave,$pos2std) = @_;
  open(my $FR,$in);
  while($_=<$FR>) {
    chomp;
    my @temp = split(/\t/);
    my $std = "NA";
    #if($temp[0]!=0) {
    if(!exists($pos2ave->{$temp[0]})) {
      if(@temp<=4) {
        if($temp[1]=~/\d+/) {
          if($temp[3] ne "NA" && $temp[1] ne "NA") {
            $std = $temp[3] - $temp[1];
            $pos2ave->{$temp[0]} = $temp[1];
            $pos2std->{$temp[0]} = $std;
          }
        }
      }
      else {
        my $ic = 0;
        for(my $i=1;$i<=7;$i+=3) {
          my $tpos = $temp[0] + $ic*0.2;
          if($temp[$i]=~/\d+/) {
            if($temp[$i+2] ne "NA" && $temp[$i] ne "NA") {
              $std = $temp[$i+2] - $temp[$i];
              $pos2ave->{$tpos} = $temp[$i];
              $pos2std->{$tpos} = $std;
            }
          }
          $ic += 1;
          #print "$tpos\n";
        }
      }
    }
  }
  close($FR);
}


