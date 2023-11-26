#!/usr/bin/perl -w
#Author: Feng Cheng
use strict;
use warnings;

my $in = $ARGV[0];  #-"BraV2.0PacBioGene.gff.100utr.sort.intergenic.TEs.step";
my $in1 = $ARGV[1]; #-"3genomes_full_tandem_AKBr_v1.2" or gene list;
my $afix = $ARGV[2]; 
my $col = $ARGV[3];
my $out = $in1.".$afix$col.stat";

&process($in,$in1,$out,$col);

sub process {
  my ($in,$in1,$out) = @_;

  my %gene2list;
  &readGene($in1,\%gene2list,$col);

  my (%pos2fdata,%pos2gdata);
  &readData($in,\%pos2fdata,\%pos2gdata,\%gene2list,"ALL");
  &output($out,\%pos2fdata,\%pos2gdata,"ALL");

  #my %gene2sub;
  #&readSub($in1,\%gene2sub);

  #my %pos2data;
  #&readData($in,\%pos2data,\%gene2sub,"ALL");
  #&output($out,\%pos2data,"ALL");
  #
  #%pos2data = ();
  #&readData($in,\%pos2data,\%gene2sub,"LF");
  #&output($out,\%pos2data,"LF");
  #
  #%pos2data = ();
  #&readData($in,\%pos2data,\%gene2sub,"MF1");
  #&output($out,\%pos2data,"MF1");
  #
  #%pos2data = ();
  #&readData($in,\%pos2data,\%gene2sub,"MF2");
  #&output($out,\%pos2data,"MF2");
}

sub output {
  my ($out,$pos2fdata,$pos2gdata,$sub) = @_;

  $out .= ".$sub";
  open(my $FW,">$out");

  my %sub2count;
  $sub2count{$sub} = 0;
  foreach my $pos (sort {$a<=>$b} keys %$pos2fdata) {
    my ($aveData,$ldata,$rdata) = &getStat($pos2fdata,$pos,\%sub2count,$sub);
    print $FW "$pos\t$aveData\t$ldata\t$rdata\n";
  }

  my (%pos2gt,%pos2ge,%pos2gi);
  foreach my $pos (sort {$a<=>$b} keys %$pos2gdata) {
    &getDataGS($pos2gdata,$pos,\%pos2gt,\%pos2ge,\%pos2gi);

    my ($aveData,$ldata,$rdata) = &getStat(\%pos2gt,$pos,\%sub2count,$sub);
    print $FW "$pos\t$aveData\t$ldata\t$rdata";

    ($aveData,$ldata,$rdata) = &getStatG(\%pos2ge,$pos);
    print $FW "\t$aveData\t$ldata\t$rdata";
    ($aveData,$ldata,$rdata) = &getStatG(\%pos2gi,$pos);
    print $FW "\t$aveData\t$ldata\t$rdata\n";
  }

  close($FW);
  system("sort -k1,1n $out > $out.sort");
}

sub getStatG {
  my ($pos2data,$pos) = @_;

  my @data = split(/-/,$pos2data->{$pos});
  my $isize = 0;  #-@data----
  my $sumData = &sumdata(\@data,\$isize,11);
  my ($aveData,$ldata,$rdata) = ("NA","NA","NA");
  if($isize!=0) {
    $aveData = $sumData/$isize;
    my $std = &stddata(\@data);
    ($ldata,$rdata) = ($aveData-$std,$aveData+$std);
  }

  return $aveData,$ldata,$rdata;
}

sub getDataGS {
  my ($pos2data,$pos,$pos2gt,$pos2ge,$pos2gi) = @_;

  my @data = split(/-/,$pos2data->{$pos});
  for(my $i=0;$i<@data;$i+=1) {
    my @dts = split(/;/,$data[$i]);

    #if()                                                      #-considering NA data, which were not summed in following analysis
    if(!exists($pos2gt->{$pos})) {
      $pos2gt->{$pos} = $dts[0];
    }
    else {
      $pos2gt->{$pos} .= "-$dts[0]";
    }
    if(!exists($pos2ge->{$pos})) {
      $pos2ge->{$pos} = $dts[1];
    }
    else {
      $pos2ge->{$pos} .= "-$dts[1]";
    }
    if(!exists($pos2gi->{$pos})) {
      $pos2gi->{$pos} = $dts[2];
    }
    else {
      $pos2gi->{$pos} .= "-$dts[2]";
    }
  }
}

sub getStat {
  my ($pos2data,$pos,$sub2count,$sub) = @_;

  my @data = split(/-/,$pos2data->{$pos});
  my $isize = 0;  #@data;
  my $sumData = &sumdata(\@data,\$isize,22);
  $sub2count->{$sub} = $isize;
  my ($aveData,$ldata,$rdata) = ("NA","NA","NA");
  if($isize!=0) {
    $aveData = $sumData/$isize;
    my $std = &stddata(\@data);
    if($std ne "NA") {
      ($ldata,$rdata) = ($aveData-$std,$aveData+$std);
    }
  }

  return $aveData,$ldata,$rdata;
}

sub readData {
  my ($in,$pos2fdata,$pos2gdata,$gene2rec,$sub) = @_;
  
  open(my $FR,$in);
  $_ = <$FR>;
  chomp;
  my @pos = split(/\t/);
  my ($gstart,$gstop);
  &revisePos(\@pos,\$gstart,\$gstop);

  my %sub2count;
  $sub2count{$sub} = 0;
  while($_=<$FR>) {
    chomp;
    my @data = split(/\t/);
    if(exists($gene2rec->{$data[0]})) {                       #--No tandem--
      if($gene2rec->{$data[0]} eq $sub || $sub eq "ALL") {
        $sub2count{$sub} += 1;
        for(my $i=1;$i<@pos;$i+=1) {
          if($data[$i]=~/;/) {               #>=$gstart <=$gstop
            &getGdata($data[$i],$pos2gdata,$pos[$i]);
          }
          else {
            &getFdata($data[$i],$pos2fdata,$pos[$i]);
          }
        }
      }
    }
  }
  close($FR);
  print "$sub\t$sub2count{$sub}\n";
}


sub getFdata {
  my ($data,$pos2fdata,$posi) = @_;

  if($data =~ /^(\d*\.*\d+)\/(\d+)$/) {
    my ($tag,$ref) = ($1,$2);
    if($ref != 0) {
      my $rtag = $tag/$ref;
      $rtag = sprintf "%.4f",$rtag;
      if(exists($pos2fdata->{$posi})) {
        $pos2fdata->{$posi} .= "-".$rtag;
      }
      else {
        $pos2fdata->{$posi} = $rtag;
      }
    }
  }
  #else {
  #  print "222: $data\n";
  #}
}

sub getGdata {
  my ($data,$pos2gdata,$posi) = @_;

  my @ds = split(/;/,$data);
  my $tdata = "";
  for(my $j=0;$j<@ds;$j+=1) {
    if($ds[$j] =~ /^(\d*\.*\d+)\/(\d+)$/) {
      my ($tag,$ref) = ($1,$2);
      if($ref != 0) {
        my $rtag = $tag/$ref;
        $rtag = sprintf "%.4f",$rtag;
        $tdata .= $rtag.";";
      }
      else {
        $tdata .= "NA;";
        #print "what 2?? $data\n";
      }
    }
    else {
      print "$ds[$j]\n";
      exit;
    }
  }
  $tdata =~ s/;$//;
  if(exists($pos2gdata->{$posi})) {
    $pos2gdata->{$posi} .= "-$tdata";
  }
  else {
    $pos2gdata->{$posi} = $tdata;
  }
}

sub revisePos {
  my ($pos,$gstart,$gstop) = @_;

  my ($gflag,$gflage) = (0,0);
  my ($gst,$gpr,$glen);
  for(my $i=1;$i<@$pos;$i+=1) {
    if($pos->[$i]=~/^\*\d+/ && $gflag==0) {
      $pos->[$i] =~ s/\*//;
      $gflag = 1;
      $gst = $pos->[$i];
      $$gstart = $i;
    }
    elsif($pos->[$i]=~/^\*\d+/ && $gflag==1) {
      $pos->[$i] =~ s/\*//;
      $gpr = $pos->[$i];
    }
    elsif($gflag==1 && $pos->[$i]!~/^\*\d+/) {
      $glen = $gst + $gpr;
      $gflag = 0;
      $$gstop = $i - 1;
      $gflage = 1;
    }
    if($gflage==1) {
      $pos->[$i] += $glen;
    }
  }
}

sub readGene {
  my ($in,$genes,$col) = @_;
  open(my $FR,$in);
  while($_=<$FR>) {
    chomp;
    my @temp = split(/\t/);
    $genes->{$temp[$col]} = "Y";
  }
  close($FR);
}

sub readSub {
  my ($in,$gene2sub) = @_;
  open(my $FR,$in);
  while($_=<$FR>) {
    chomp;
    my @temp = split(/\t/);
    #if(($temp[4]=~/Bra/ && $temp[5]=~/Bra/) || ($temp[4]=~/Bra/ && $temp[6]=~/Bra/) || ($temp[5]=~/Bra/ && $temp[6]=~/Bra/)) {
    if($temp[4]=~/Bra/) {    #--Tandem genes?
      $gene2sub->{$temp[4]} = "LF";
    }
    if($temp[5]=~/Bra/) {
      $gene2sub->{$temp[5]} = "MF1";
    }
    if($temp[6]=~/Bra/) {
      $gene2sub->{$temp[6]} = "MF2";
    }
    #}
  }
  close($FR);
}

sub sumdata {
  my ($data,$cc,$lab) = @_;

  my $sumData = 0;
  for(my $i=0;$i<@$data;$i+=1) {
    if($data->[$i]=~/\d/) {
      if($data->[$i]=~/e/) {
        print "$lab\n@$data\n";
        exit;
      }
      $sumData+=$data->[$i];
      $$cc += 1;
    }
  }
  return $sumData;
}

sub stddata {
  my ($data) = @_;
  my $isize = 0;  #@data;
  my $sumData = &sumdata($data,\$isize,33);
  my $aveData = $sumData/$isize;
  my $stdData = 0;
  my $count = 0;
  for(my $i=0;$i<@$data;$i+=1) {
    if($data->[$i]=~/\d/) {
      $stdData += ($data->[$i]-$aveData)*($data->[$i]-$aveData);
      $count += 1;
    }
  }
  if($count<=1) {
    $stdData = "NA";
  }
  else {
    $stdData /= ($count-1);
    $stdData = sqrt($stdData);
  }
  return $stdData;
}


