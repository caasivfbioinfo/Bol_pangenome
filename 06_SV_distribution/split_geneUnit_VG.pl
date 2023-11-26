#!/usr/bin/perl -w
#Author: Feng Cheng
use warnings;
use strict;

my $in1 = $ARGV[0];               # scaffold_gene.gff3.new_chr0708_v1.2 or gene position sort
my $in2 = $ARGV[1];               # scaffold_newset.info.sorted0708.chr.fa.ChrSize;
my $length = $ARGV[2];            # 1500;
my $out1 = $in1.".".$length.".geneUnits";

&process($in1,$in2,$length,$out1);

sub process {
  my ($in1,$in2,$length,$out1) = @_;

  my ($index,%index2gene,%gene2strand,%gene2start,%gene2stop,%gene2chr,%pos2pgene,%pos2agene,%chr2end,%chr2seq);
  my $len5 = $length;
  my $len3 = $length;

  my %mrna2cdsPos;
  &readGff($in1,\$index,\%index2gene,\%gene2chr,\%gene2strand,\%gene2start,\%gene2stop,\%mrna2cdsPos);
  &sortCDS(\%mrna2cdsPos);  
  #&readGpos($in1,\$index,\%index2gene,\%gene2chr,\%gene2strand,\%gene2start,\%gene2stop);

  &readChr($in2,\%chr2end);

  &output($out1,$len5,$len3,\$index,\%index2gene,\%gene2chr,\%gene2strand,\%gene2start,\%gene2stop,\%chr2end,\%chr2seq,\%mrna2cdsPos);

}

sub output {
  my ($out1,$len5,$len3,$index,$index2gene,$gene2chr,$gene2strand,$gene2start,$gene2stop,$chr2end,$chr2seq,$mrna2cdsPos) = @_;

  my ($strand,$gene,$pos,$chr);
  open(my $FW,">$out1");
  print $FW "#Gene\tChr\tStrand\tLeftMost\tL2GenesMid\tGeneStart\tGeneStop\tR2GenesMid\tRightMost\tLoverlap\tRoverlap\tCDs\n";
  for(my $i=0;$i<$$index;$i+=1) {
    $gene = $index2gene->{$i};
    $strand = $gene2strand->{$gene};
    my $start = $gene2start->{$gene};
    my $stop = $gene2stop->{$gene};
    $chr = $gene2chr->{$gene};
    my ($rstart,$rstop) = (0,0);
    my ($rrstart,$rrstop) = (0,0);
    my ($lover,$rover) = (0,0);
    my ($pgene,$agene);
    
    if($i==0) {
      $rstart = $start-1-$len5+1;
      if($rstart<=0) {
        $rstart = 1;
      }
      $rrstart = 1;
    }
    else {
      my $pstop = $start;
      my $ii = $i;
      while($pstop>=$start) {                         # pay attention-----
        $ii -= 1;
        $pgene = $index2gene->{$ii};
        $pstop = $gene2stop->{$pgene};
        if($gene2chr->{$pgene} ne $chr) {
          $pstop = $start - 1;
          $rstart = $start-1-$len5+1;
          if($rstart<=0) {
            $rstart = 1;
          }
          $rrstart = 1;
        }
      }
      if($rstart==0) {
        if($start-1-($pstop+1)+1<2*$len5) {
          my $tlen = $start-1-($pstop+1)+1;
          $rstart = int($start-1 - $tlen*($len5/($len5+$len3)) +1);
        }
        else {
          $rstart = $start-1 - $len5+1;
        }
        $rrstart = $pstop + 1;
      }
      $pgene = $index2gene->{$i-1};
      $pstop = $gene2stop->{$pgene};
      if($pstop>=$start && $gene2chr->{$pgene} eq $chr) {
        $lover = 1;
      }
    }

    if($i==$$index-1) {
      $rstop = $stop+1+$len3-1;
      if($rstop>$chr2end->{$chr}) {
        $rstop = $chr2end->{$chr};
      }
      $rrstop = $chr2end->{$chr};
    }
    else {
      my $astart = $stop;
      my $ii = $i;
      while($astart<=$stop) {
        $ii += 1;
        my $agene = $index2gene->{$ii};
        $astart = $gene2start->{$agene};
        if($gene2chr->{$agene} ne $chr) {
          $astart = $stop + 1;
          $rstop = $stop+1+$len3-1;
          if($rstop>$chr2end->{$chr}) {
            $rstop = $chr2end->{$chr};
          }
          $rrstop = $chr2end->{$chr};
        }
      }
      if($rstop==0) {
        if($astart-1-($stop+1)+1<2*$len3) {
          my $tlen = $astart-1-($stop+1)+1;
          $rstop = int($stop + $tlen*($len3/($len5+$len3)));
        }
        else {
          $rstop = $stop+1 + $len3-1;
        }
        $rrstop = $astart - 1;
      }
      $agene = $index2gene->{$i+1};
      $astart = $gene2start->{$agene};
      if($astart<=$stop && $gene2chr->{$agene} eq $chr) {
        $rover = 1;
      }
    }

    #print $FW "$gene\t$chr\t$strand\t$rstart\t$start\t$stop\t$rstop\n";
    print $FW "$gene\t$chr\t$strand\t$rrstart\t$rstart\t$start\t$stop\t$rstop\t$rrstop\t$lover\t$rover";

    my $temp = $mrna2cdsPos->{$gene};
    print $FW "\t$temp";
    #my @pos = split(/;/,$temp);
    #for(my $i=0;$i<@pos;$i+=1) {
    #  @subpos = split(/,/,$pos[$i]);
    #}
    print $FW "\n";
  }
  close($FW);

  my $out2 = $out1.".sort";

  #system("sort -k2,2 -k5,5n -k6,6n $out1 > $out2");
}

sub readChr {
  my ($in2,$chr2end) = @_;

  open(my $FR2,$in2);
  while($_=<$FR2>) {
    chomp;
    my @temp = split(/\t/);
    $chr2end->{$temp[0]} = $temp[1];
  }
  close($FR2);
}

sub readGpos {
  my ($in1,$index,$index2gene,$gene2chr,$gene2strand,$gene2start,$gene2stop) = @_;

  $$index = 0;
  my ($strand,$gene,$pos,$chr);
  open(my $FR1,$in1);
  while($_=<$FR1>) {
    chomp;
    my @temp = split(/\t/);
    $gene = $temp[0];
    $index2gene->{$$index} = $gene;
    $$index += 1;
    $gene2chr->{$gene} = $temp[1];
    $gene2strand->{$gene} = $temp[4];
    $strand = $temp[4];
    if($temp[2]<$temp[3]) {
      $gene2start->{$gene} = $temp[2];
      $gene2stop->{$gene} = $temp[3];
    }
    else {
      $gene2start->{$gene} = $temp[3];
      $gene2stop->{$gene} = $temp[2];
    }
  }
  close($FR1);
}

sub readGff {
  my ($in1,$index,$index2gene,$gene2chr,$gene2strand,$gene2start,$gene2stop,$mrna2cdsPos) = @_;

  $$index = 0;
  my ($strand,$gene,$mrna,$pos,$chr,%mrna2rec);
  open(my $FR1,$in1);
  while($_=<$FR1>) {
  #if(/\tgene\t/) {
    if(/\tmRNA\t/) {
      chomp;
      my @temp = split(/\t/);
      $temp[8] =~ /ID=([^;\n]+)/;
      $mrna = $1;
      if(/Parent=([^;\n]+)/) {        #depends
        $gene = $1;
      }
      else {
        $gene = $mrna;
      }
      $index2gene->{$$index} = $mrna;
      $$index += 1;
      $gene2chr->{$mrna} = $temp[0];
      $gene2strand->{$mrna} = $temp[6];
      $strand = $temp[6];
      if($temp[3]<$temp[4]) {
        $gene2start->{$mrna} = $temp[3];
        $gene2stop->{$mrna} = $temp[4];
      }
      else {
        $gene2start->{$mrna} = $temp[4];
        $gene2stop->{$mrna} = $temp[3];
      }
      $mrna2rec{$mrna} = $gene;
    }
    if(/\tCDS\t/) {
      if(/Parent=([^;\n,\s]+)/) {   #-it depends----
        $mrna = $1;
        if(exists($mrna2rec{$mrna})) {
          my @temp = split(/\t/);
          #my $pos = $temp[0].",".$temp[3].",".$temp[4];
          my $pos = $temp[3]."-".$temp[4];
          if(!exists($mrna2cdsPos->{$mrna})) {
            $mrna2cdsPos->{$mrna} = $pos;
          }
          else {
            $mrna2cdsPos->{$mrna} .= ";".$pos;
          }
        }
      }
    }
  }
  close($FR1);
}

sub sortCDS {
  my ($gene2cdsPos) = @_;
  foreach my $gene (keys %$gene2cdsPos) {
    my $cdsPos = $gene2cdsPos->{$gene};
    my @pos = split(/;/,$cdsPos);
    my @temp1 = map{[$_,((split(/-/))[0]+(split(/-/))[1])/2]} @pos;
    my @temp2 = sort{$a->[1] <=> $b->[1]} @temp1;
    my @sorted = map{$_->[0]} @temp2;
    $gene2cdsPos->{$gene} = join(";",@sorted);
  }
  print "Sort CDS is OK!\n";
}

