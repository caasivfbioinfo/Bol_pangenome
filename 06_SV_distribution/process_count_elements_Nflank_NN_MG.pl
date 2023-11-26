#!/usr/bin/perl -w
#Author: Feng Cheng
use warnings;
use strict;
#require("/home/chengf/bin/perl_subs.pl");

my $in1 = $ARGV[0];          #-elements position information--
my $in2 = $ARGV[1];          #-elements list and distances information--
my $in3 = $ARGV[2];          #-genomes.fa--
my $ele = $ARGV[3];
my $afix = $ARGV[4];
my ($upsize,$downsize) = (5000,5000);    #-5000 for genes, 3000 for TEs/SVs--
my ($win,$step) = (100,20);
my ($gcnum,$gwinf) = (20,$win/$step);    #-each element cut to $gcnum bins--
#my $out = $in2.".$ele.".$afix;
my $out = $in1.".$ele.".$afix;
$out =~ s/^.*\///;

&process($in1,$in2,$in3,$out,$upsize,$downsize,$win,$step,$gcnum,$gwinf,$ele);

sub process {
  my ($in1,$in2,$in3,$out,$upsize,$downsize,$win,$step,$gcnum,$gwinf,$ele) = @_;
  
  my (%pos2tags,$tagnum);
  &readElements($in1,\%pos2tags,\$tagnum,$ele);

  my %sca2seq;
  &readFasta($in3,\%sca2seq);
  
  &output($in2,$out,\%pos2tags,\%sca2seq,$upsize,$downsize,$win,$step,$tagnum,$gcnum,$gwinf);
}


sub output {
  my ($in,$out,$pos2tags,$sca2seq,$upsize,$downsize,$win,$step,$tagnum,$gcnum,$gwinf) = @_;
  open(my $FR,$in);
  open(my $FW,">$out");
 
  #--output the position index as the headline---------------
  print $FW "$tagnum";
  my ($start,$stop,$mid);
  my $ss = 0 - $upsize - 1;
  for(my $i=$ss;$i<0;$i+=$step) {
    $start = $i - int($win/2);
    if($start<$ss) {
      $start = $ss;
    }
    $stop = $i + int($win/2);
    if($stop>=0) {
      $stop = -1;
    }
    $mid = int(($start+$stop)/2);
    print $FW "\t$mid";
  }

  my $slen = 2000;          #do not change the numbner, it's for generating intervals--2000 for gene, 1000 for TEs/SVs--------
  my ($sstep,$swing) = &getInterval($slen,$gcnum,$gwinf);
  for(my $i=1;$i-1<=$slen;$i+=$sstep) {
    my ($start,$stop) = &getWinStep($i,$swing,$slen);
    my $mid = int(($start+$stop)/2);
    print $FW "\t*$mid";
  }

  $ss = 1;
  for(my $i=$ss;$i<=$downsize+1;$i+=$step) {
    $start = $i - int($win/2);
    if($start<$ss) {
      $start = $ss;
    }
    $stop = $i + int($win/2);
    if($stop>$downsize+1) {
      $stop = $downsize+1;
    }
    $mid = int(($start+$stop)/2);
    print $FW "\t$mid";
  }
  print $FW "\n";

  #--compute gene one by one-----------------
  my $gcc = 0;
  while($_=<$FR>) {
    if(!/^#/) {
      $gcc += 1;
      #print "calculationg the $gcc genes...\n";
      chomp;
      my @temp = split(/\t/);
      print $FW "$temp[0]";
      my ($left,$right) = ($temp[3],$temp[8]);
      my ($sca,$ori,$gstart,$gstop) = ($temp[1],$temp[2],$temp[5],$temp[6]);
      my $cds = $temp[11];                    #--for multiple exons-------
      #my $cds = $gstart."-".$gstop;          #--for TEs/SVs------------------
      my $pout = "";

      &winScan($gstart,$gstop,$cds,$upsize,$downsize,$win,$step,$sca,$ori,$sca2seq,\$pout,$left,$right,$pos2tags,$gcnum,$gwinf);

      print $FW $pout;
      my @aa = split(/\t/,$pout);
    }
  }
  close($FR);
  close($FW);
}


sub winScan {
  my ($gstart,$gstop,$cds,$upsize,$downsize,$win,$step,$sca,$ori,$sca2seq,$pout,$left,$right,$pos2tags,$gcnum,$gwinf) = @_;
  
  if($ori eq "+") {
    &foreward5($gstart,$upsize,$win,$step,$sca,$sca2seq,$pout,$left,$pos2tags);          # 5'

    &forewardG($cds,$sca,$ori,$sca2seq,$pout,$pos2tags,$gcnum,$gwinf);                   # forward gene

    &foreward3($gstop,$downsize,$win,$step,$sca,$sca2seq,$pout,$right,$pos2tags);        # 3'
  }
  elsif($ori eq "-") {
    &backward5($gstop,$upsize,$win,$step,$sca,$sca2seq,$pout,$right,$pos2tags);          # 5'

    &backwardG($cds,$sca,$ori,$sca2seq,$pout,$pos2tags,$gcnum,$gwinf);                   # back gene

    &backward3($gstart,$downsize,$win,$step,$sca,$sca2seq,$pout,$left,$pos2tags);        # 3'
  }
  else {
    print "$_\n";
  }
  
}

sub backwardG {
  my ($cds,$sca,$ori,$sca2seq,$pout,$pos2tags,$gcnum,$gwinf) = @_;

  my (@exons,@introns,@genes);
  &extractCDs($cds,\@exons,\@introns,\@genes);
  my $elen = @exons;
  my $ilen = @introns;
  my $tlen = $elen + $ilen;
  
  my ($tstep,$twing) = &getInterval($tlen,$gcnum,$gwinf);
  my ($estep,$ewing) = &getInterval($elen,$gcnum,$gwinf);
  my ($istep,$iwing) = &getInterval($ilen,$gcnum,$gwinf);

  my $cc = 0;
  my ($j,$k) = ($elen-1,$ilen-1);
  for(my $i=$tlen-1;$i+2>0;$i-=$tstep) {
    $$pout .= "\t";
    my ($stop,$start) = &getWinStep($i,$twing,$tlen);
    &gbackScan($start,$stop,\@genes,$sca2seq,$pos2tags,$sca,$pout);

    $$pout .= ";";
    $j -= $estep;
    ($stop,$start) = &getWinStep($j,$ewing,$elen);
    &gbackScan($start,$stop,\@exons,$sca2seq,$pos2tags,$sca,$pout);

    $$pout .= ";";
    $k -= $istep;
    ($stop,$start) = &getWinStep($k,$iwing,$ilen);
    &gbackScan($start,$stop,\@introns,$sca2seq,$pos2tags,$sca,$pout);

    $cc += 1;
  }
  #print "--$cc\n";

}

sub gbackScan {
  my ($start,$stop,$genes,$sca2seq,$pos2tags,$sca,$pout) = @_;

  my ($realWin,$ntags) = (0,0);
  for(my $i=$start;$i>=$stop;$i-=1) {
    my $posi = $genes->[int($i)];                            # get integer here---

    my $nucl = substr($sca2seq->{$sca},$posi-1,1);

    &calTags($nucl,\$realWin,$pos2tags,\$ntags,$sca,$posi);
  }
  
  $$pout .= "$ntags/$realWin";
}

sub forewardG {
  my ($cds,$sca,$ori,$sca2seq,$pout,$pos2tags,$gcnum,$gwinf) = @_;

  my (@exons,@introns,@genes);
  &extractCDs($cds,\@exons,\@introns,\@genes);
  my $elen = @exons;
  my $ilen = @introns;
  my $tlen = $elen + $ilen;
  
  my ($tstep,$twing) = &getInterval($tlen,$gcnum,$gwinf);
  my ($estep,$ewing) = &getInterval($elen,$gcnum,$gwinf);
  my ($istep,$iwing) = &getInterval($ilen,$gcnum,$gwinf);

  my $cc = 0;
  my ($j,$k) = (0,0);
  for(my $i=0;$i-1<$tlen;$i+=$tstep) {
    $$pout .= "\t";
    my ($start,$stop) = &getWinStep($i,$twing,$tlen);
    &gforeScan($start,$stop,\@genes,$sca2seq,$pos2tags,$sca,$pout);

    $$pout .= ";";
    $j += $estep;
    ($start,$stop) = &getWinStep($j,$ewing,$elen);
    &gforeScan($start,$stop,\@exons,$sca2seq,$pos2tags,$sca,$pout);

    $$pout .= ";";
    $k += $istep;
    ($start,$stop) = &getWinStep($k,$iwing,$ilen);
    &gforeScan($start,$stop,\@introns,$sca2seq,$pos2tags,$sca,$pout);

    $cc += 1;
  }

}

sub gforeScan {
  my ($start,$stop,$genes,$sca2seq,$pos2tags,$sca,$pout) = @_;

  my $ng = @$genes;
  my ($realWin,$ntags) = (0,0);
  for(my $i=$start;$i<=$stop;$i+=1) {
    my $posi = $genes->[int($i)];                            # get integer here---

    my $nucl = substr($sca2seq->{$sca},$posi-1,1);

    &calTags($nucl,\$realWin,$pos2tags,\$ntags,$sca,$posi);
  }
  
  $$pout .= "$ntags/$realWin";
}

sub getWinStep {
  my ($i,$wing,$len) = @_;

  my ($start,$stop);
  $start = $i - $wing/2;
  if($start<0) {
    $start = 0;
  }
  $stop = $i + $wing/2;
  if($stop>=$len) {
    $stop = $len-1;
  }

  return $start,$stop;
}

sub getInterval {
  my ($len,$gcnum,$gwinf) = @_;

  my $step = $len/$gcnum;
  my $wing = $step*$gwinf;

  return $step,$wing;
}

sub extractCDs {
  my ($cds,$exons,$introns,$genes) = @_;

  my @cdss = split(/;/,$cds);                    #considering sort cds--------------
  my ($ii,$jj,$kk,$is) = (0,0,0,"NA");
  for(my $i=0;$i<@cdss;$i+=1) {
    my @ss = split(/-/,$cdss[$i]);
    if($is ne "NA") {
      for(my $k=$is;$k<$ss[0];$k+=1) {
        $introns->[$jj] = $k;
        $jj += 1;

        $genes->[$kk] = $k;
        $kk += 1;
      }
    }
    my $j;
    for($j=$ss[0];$j<=$ss[1];$j+=1) {
      $exons->[$ii] = $j;
      $ii += 1;

      $genes->[$kk] = $j;
      $kk += 1;
    }
    $is = $j;
  }
}


sub calTags {
  my ($nucl,$realWin,$pos2tags,$ntags,$sca,$i,$pout) = @_;

  if($nucl=~/[ATCGatcg]/) {
    $$realWin += 1;
    my $pos = $sca.":".$i;
    if(exists($pos2tags->{$pos})) {
      $$ntags += 1;                     #--its depend--site counts----------------
      #$$ntags += $pos2tags->{$pos};      #--its depend--site level adding----------------
    }
  }
}


sub compCount {
  my ($seq,$realWin,$pos2tags,$ntags,$sca,$i)= @_;

  $seq =~ s/^(\w)//;
  my $nucl = $1;

  &calTags($nucl,$realWin,$pos2tags,$ntags,$sca,$i);

}

sub backScan {
  my ($pout,$start,$stop,$win,$sca2seq,$pos2tags,$sca) = @_;

  my ($lenW,$seq,$realWin,$ntags);
  $lenW = $win;
  $seq = substr($sca2seq->{$sca},$stop-1,$lenW);
  $seq = reverse($seq);
  $realWin = 0;
  $ntags = 0;
  for(my $i=$start;$i>=$stop;$i-=1) {
    &compCount($seq,\$realWin,$pos2tags,\$ntags,$sca,$i);
  }
  $$pout .= "$ntags/$realWin";
}

sub backward3 {
  my ($gstart,$downsize,$win,$step,$sca,$sca2seq,$pout,$left,$pos2tags) = @_;

  my ($start,$stop);
  my $ss = $gstart - 1;
  for(my $i=$ss;$i>=$gstart-$downsize-1;$i-=$step) {
    $start = $i + int($win/2);
    if($start>$ss) {
      $start = $ss;
    }
    $stop = $i - int($win/2);
    if($stop<$gstart-$downsize-1) {
      $stop = $gstart-$downsize-1;
    }
    $$pout .= "\t";
    if($stop>=$left) {
      &backScan($pout,$start,$stop,$win,$sca2seq,$pos2tags,$sca);
    }
    else {
      $$pout .= "-/-";
    }
  }
  $$pout .= "\n";
}

sub backward5 {
  my ($gstop,$upsize,$win,$step,$sca,$sca2seq,$pout,$right,$pos2tags) = @_;

  my ($start,$stop);
  my $ss = $gstop + $upsize+1;
  for(my $i=$ss;$i>$gstop;$i-=$step) {
    $start = $i + int($win/2);
    if($start>$ss) {
      $start = $ss;
    }
    $stop = $i - int($win/2);
    if($stop<=$gstop) {
      $stop = $gstop+1;
    }
    $$pout .= "\t";
    if($start<=$right) {
      &backScan($pout,$start,$stop,$win,$sca2seq,$pos2tags,$sca);
    }
    else {
      $$pout .= "-/-";
    }
  }
}

sub foreScan {
  my ($pout,$start,$stop,$win,$sca2seq,$pos2tags,$sca) = @_;

  my ($lenW,$seq,$realWin,$ntags);
  $lenW = $win;
  $seq = substr($sca2seq->{$sca},$start-1,$lenW);
  $realWin = 0;
  $ntags = 0;
  for(my $i=$start;$i<=$stop;$i+=1) {
    &compCount($seq,\$realWin,$pos2tags,\$ntags,$sca,$i);
  }
  $$pout .= "$ntags/$realWin";

}

sub foreward3 {
  my ($gstop,$downsize,$win,$step,$sca,$sca2seq,$pout,$right,$pos2tags) = @_;

  my ($start,$stop);
  my $ss = $gstop + 1;
  for(my $i=$ss;$i<=$gstop+$downsize+1;$i+=$step) {
    $start = $i - int($win/2);
    if($start<$ss) {
      $start = $ss;
    }
    $stop = $i + int($win/2);
    if($stop>$gstop+$downsize+1) {
      $stop = $gstop+$downsize+1;
    }
    $$pout .= "\t";
    if($stop<=$right) {
      &foreScan($pout,$start,$stop,$win,$sca2seq,$pos2tags,$sca);
    }
    else {
      $$pout .= "-/-";
    }
  }
  $$pout .= "\n";
}

sub foreward5 {
  my ($gstart,$upsize,$win,$step,$sca,$sca2seq,$pout,$left,$pos2tags) = @_;

  my ($start,$stop);
  my $ss = $gstart - $upsize-1;
  for(my $i=$ss;$i<$gstart;$i+=$step) {
    $start = $i - int($win/2);
    if($start<$ss) {
      $start = $ss;
    }
    $stop = $i + int($win/2);
    if($stop>=$gstart) {
      $stop = $gstart-1;
    }
    $$pout .= "\t";
    if($start>=$left) {
      &foreScan($pout,$start,$stop,$win,$sca2seq,$pos2tags,$sca);
    }
    else {
      $$pout .= "-/-";
    }
  }
}


sub readElements {
  my ($in,$pos2tags,$tagnum,$ele) = @_;

  $$tagnum = 0;
  open(my $FR,$in);
  while($_=<$FR>) {
    chomp;
    my @temp = split(/\t/);
    #my @sdata = split(/:/,$temp[9]);
    #if($ele eq "ALL" || $ele=~/$temp[4]/) {
    if(/DEL/ || /CPL/) {
      #my $i = $temp[1];
      #for(my $i=$temp[1];$i<=$temp[2];$i+=1) {        #--it's depend------------
      #for(my $i=$temp[3];$i<=$temp[4];$i+=1) {        #--it's depend------------
      #for(my $i=$temp[5];$i<=$temp[6];$i+=1) {        #--it's depend------------
      for(my $i=$temp[7];$i<=$temp[8];$i+=1) {        #--it's depend------------
        #my $pos = $temp[0].":".$i;
        #my $pos = $temp[1].":".$i;
        my $pos = $temp[6].":".$i;
        if(!exists($pos2tags->{$pos})) {
          $pos2tags->{$pos} = 1;                      #--it's depend------------
          #$pos2tags->{$pos} = $temp[3]/$temp[2];       #--it's depend--site weighted----------
        }
        #else {
        #  $pos2tags->{$pos} += 1;
        #}
      }
      $$tagnum += 1;
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
      $id2seq->{$id} = "";
    }
    else {
      chomp;
      $id2seq->{$id} .= $_;
    }
  }
  close($SFR);
}



