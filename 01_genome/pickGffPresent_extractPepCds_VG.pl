#!/usr/bin/perl -w
#Author: Feng Cheng
use strict;
use warnings;
use Getopt::Long;
#require "/home/chengf/bin/perl_subs.pl";

my $gff = $ARGV[0];  #-gff--
my $genoRef = $ARGV[1];
my $out1 = $gff.".repr";
my $out11 = $gff.".repr.position";
my $out2 = $gff.".repr.cds";
my $out3 = $gff.".repr.pep";
#system("awk '($3=="mRNA" || $3=="CDS")' $gff > $gff.temp");
#system("mv $gff.temp $gff");

my %co2aa = ("TTT"=>"F","TTC"=>"F","TTA"=>"L","TTG"=>"L","CTA"=>"L","CTT"=>"L","CTC"=>"L","CTG"=>"L","ATT"=>"I",
             "ATC"=>"I","ATA"=>"I","ATG"=>"M","GTT"=>"V","GTC"=>"V","GTA"=>"V","GTG"=>"V","TCT"=>"S","TCC"=>"S",
             "TCA"=>"S","TCG"=>"S","AGC"=>"S","AGT"=>"S","CCT"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P","ACT"=>"T",
             "ACC"=>"T","ACA"=>"T","ACG"=>"T","GCT"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A","TAT"=>"Y","TAC"=>"Y",
             "CAT"=>"H","CAC"=>"H","CAA"=>"Q","CAG"=>"Q","AAT"=>"N","AAC"=>"N","AAA"=>"K","AAG"=>"K","GAT"=>"D",
             "GAC"=>"D","GAA"=>"E","GAG"=>"E","TGT"=>"C","TGC"=>"C","TGG"=>"W","CGT"=>"R","CGC"=>"R","CGA"=>"R",
             "CGG"=>"R","AGA"=>"R","AGG"=>"R","GGT"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G",
             "TAA"=>"*","TAG"=>"*","TGA"=>"*"); #corrected


&process($gff,$genoRef,$out1,$out11,$out2,$out3);

sub process {
  my ($gff,$genoRef,$out1,$out11,$out2,$out3) = @_;

  my %sca2seq;
  &readGenoFas($genoRef,\%sca2seq);

  my (%mrna2cdsNum,%mrna2orient,%mrna2cdsPos,%gene2mrna);
  &readGff($gff,\%mrna2cdsNum,\%mrna2orient,\%mrna2cdsPos,\%gene2mrna);

  &sortCDS(\%mrna2cdsPos);

  &getSeq($gff,\%mrna2cdsPos,\%mrna2orient,\%mrna2cdsNum,\%gene2mrna,\%sca2seq,$out1,$out11,$out2,$out3);

}

sub getSeq {
  my ($in,$mrna2cdsPos,$mrna2orient,$mrna2cdsNum,$gene2mrna,$sca2seq,$out1,$out11,$out2,$out3) = @_;

  my %mrna2repr;
  my ($gene,$mrna);
  foreach $gene (keys %$gene2mrna) {
    $mrna = $gene2mrna->{$gene};
    if(!exists($mrna2repr{$mrna})) {
      $mrna2repr{$mrna} = $gene;
    }
    else {
      print "what??? $mrna\t$mrna2repr{$mrna}\n";
    }
  }

  open(my $FR,$in);
  open(my $FW,">$out1");
  open(my $FW0,">$out11");
  while($_=<$FR>) {
    my @temp = split(/\t/);
    if(/\tgene\t/) {
      print $FW $_;
    }
    #if(/\tgene\t/) {
      #/ID=gene:([^;\n,]+)/;           #-depends--
    elsif(/\tmRNA\t/) {           #-depends--
      #/mRNA ([^;\n,\s]+)/;           #-depends--
      /ID=([^;\n,]+)/;           #-depends--
      #/Name=([^;\n,]+)/;           #-depends--
      $mrna = $1;
      if(exists($mrna2repr{$mrna})) {
        print $FW $_;
        print $FW0 "$mrna\t$temp[0]\t$temp[3]\t$temp[4]\t$temp[6]\n";
      }
    }
    elsif(/\tCDS\t/) {
      #/CDS ([^;\n,\s]+)/;  #-depends--
      #/Parent=gene:([^;,\n]+)/;  #-depends--
      #/Parent=mRNA:([^;,\n]+)/;  #-depends--
      /Parent=([^;,\n]+)/;  #-depends--
      $mrna = $1;
      if(exists($mrna2repr{$mrna})) {
        print $FW $_;
      }
    }
  }
  close($FR);
  close($FW);
  close($FW0);
  system("sort -k2,2 -k3,3n -k4,4n $out11 > $out11.sort");

  open(my $FW1,">$out2");
  open(my $FW2,">$out3");
  foreach my $mrna (sort keys %mrna2repr) {
    my $temp = $mrna2cdsPos->{$mrna};
    my @pos = split(/;/,$temp);
    my $ori = $mrna2orient->{$mrna};
    my $num = $mrna2cdsNum->{$mrna};
    my ($sca,$seq);
    my $start = 0;
    my $stop = 0;
    my $tseq = "";
    my @subpos;
    for(my $i=0;$i<@pos;$i+=1) {
      @subpos = split(/,/,$pos[$i]);
      $sca = $subpos[0];
      $start = $subpos[1];
      my $starti = $start - 1;
      $stop = $subpos[2];
      my $len = $stop - $start + 1;
      if(exists($sca2seq->{$sca})) {
        $seq = substr($sca2seq->{$sca},$starti,$len);
        $tseq .= $seq;
      }
      else {
        print "$sca not existed?\n";
      }
    }
    if($ori eq "-") {
      $tseq =~ tr/ATCG/TAGC/;
      $tseq = reverse($tseq);
    }

    my $cds = $tseq;
    my $pseq = "";
    while($tseq=~s/^(\w{3})//) {
      if(exists($co2aa{$1})) {
        $pseq .= $co2aa{$1};
      }
      else {
        #print $1,"\n";
        $pseq .= "X";
      }
    }

    print $FW1 ">$mrna\n$cds\n";
    print $FW2 ">$mrna\n$pseq\n";
  }
  close($FW1);
  close($FW2);
}

sub sortCDS {
  my ($gene2cdsPos) = @_;
  #open(FW,">test");
  foreach my $gene (keys %$gene2cdsPos) {
    #print FW "$gene\t$gene2orient{$gene}\t$gene2cdsNum{$gene}\t$gene2cdsPos{$gene}\n";
    my $cdsPos = $gene2cdsPos->{$gene};
    my @pos = split(/;/,$cdsPos);
    my @temp1 = map{[$_,((split(/,/))[1]+(split(/,/))[2])/2]} @pos;
    my @temp2 = sort{$a->[1] <=> $b->[1]} @temp1;
    my @sorted = map{$_->[0]} @temp2;
    $gene2cdsPos->{$gene} = join(";",@sorted);
    #print FW "$gene\t$gene2orient{$gene}\t$gene2cdsNum{$gene}\t$gene2cdsPos{$gene}\n";
  }
  #close(FW);
  print "Sort CDS is OK!\n";
}

sub readGff {
  my ($in,$mrna2cdsNum,$mrna2orient,$mrna2cdsPos,$gene2mrna) = @_;

  open(my $FR1,$in);
  my ($gene,$mrna);
  my (%mrna2rec,%mrna2len);
  while(<$FR1>) {
    #if(/\tgene\t/) {
      #/ID=gene:([^;\n]+)/;       #depends
    if(/\tmRNA\t/) {            #depends
      #/mRNA ([^;\n,\s]+)/;         #depends
      /ID=([^;\n]+)/;         #depends
      #/Name=([^;\n]+)/;         #depends
      $mrna = $1;
      #/Parent=gene:([^;\n]+)/;   #depends
      if(/Parent=([^;\n]+)/) {        #depends
        $gene = $1;
      }
      else {
        $gene = $mrna;
      }
      $mrna2rec{$mrna} = $gene;
    }
    if(/\tCDS\t/) {
      #if(/Parent=gene:([^;,\n]+)/) {   #-it depends----
      #if(/Parent=transcript:([^;\n]+)/) {   #-it depends----
      #if(/Parent=mRNA:([^;\n]+)/) {   #-it depends----
      #if(/CDS ([^;\n,\s]+)/) {   #-it depends----
      if(/Parent=([^;\n,\s]+)/) {   #-it depends----
        $mrna = $1;
        if(exists($mrna2rec{$mrna})) {
          my @temp = split(/\t/);
          my $pos = $temp[0].",".$temp[3].",".$temp[4];
          if(!exists($mrna2cdsPos->{$mrna})) {
            $mrna2cdsPos->{$mrna} = $pos;
            $mrna2cdsNum->{$mrna} = 1;
            $mrna2orient->{$mrna} = $temp[6];
            $mrna2len{$mrna} = $temp[4] - $temp[3] + 1;
          }
          else {
            $mrna2cdsPos->{$mrna} .= ";".$pos;
            $mrna2cdsNum->{$mrna} += 1;
            $mrna2len{$mrna} += $temp[4] - $temp[3] + 1;
          }
        }
      }
    }
  }
  close($FR1);

  foreach $mrna (keys %mrna2rec) {
    $gene = $mrna2rec{$mrna};
    if(!exists($gene2mrna->{$gene})) {
      $gene2mrna->{$gene} = $mrna;
    }
    else {
      my $tmrna = $gene2mrna->{$gene};
      if($mrna2len{$tmrna}<$mrna2len{$mrna}) {
        $gene2mrna->{$gene} = $mrna;
      }      
    }
  }
}

sub readGenoFas {
  my ($in,$sca2seq) = @_;
  my $sca;
  open(my $FR,$in);
  while(defined($_=<$FR>)) {
    if(/^>([^\s]+)\s*/) {
      $sca = $1;
      $sca2seq->{$sca} = '';
    }
    if(!/^>/) {
      chomp;
      $sca2seq->{$sca} .= uc($_);
      my $lenSca += length($_);
    }
  }
  close($FR);
  print "Genome fas is OK!\n";
}

