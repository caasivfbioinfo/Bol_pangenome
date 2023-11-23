#!/usr/bin/perl -w
#Author: Feng Cheng
use strict;
use warnings;
#require("/home/chengf/bin/perl_subs.pl");

my $in1 = $ARGV[0]; # five column, data in 2 to 5 column
my $out = $in1.".feTest";

&process($in1,$out);

sub process {
  my ($in1,$out) = @_;
  ##CAN ADD OTHER FORMATING HERE IS NEEDED##------------------------

  &output($in1,$out);
}


sub output {
  my ($in,$out,$At2go,$go2info,$BrT2Br,$At2Brs,$Br2At) = @_;
  
  ##CAN ADD OTHER FORMATING HERE IS NEEDED##------------------------

  &fisherR($in,$out);
}


sub fisherR {
  my ($in,$out) = @_;

  my ($in1t,$in2t) = ($in.".data",$in.".info");
  system("cut -f 1-5 $in > $in1t");
  #system("cut -f 6 $in > $in2t");

  my $out1 = $out.".p";
  my $outR = $in.".rsc";
  if(-e $out1) {
    system("rm $out1");
  }
  open(my $FW,">$outR");
  print $FW "read.table(\"$in1t\",as.is=1) -> data;\n";
  print $FW "for (i in 1:length(data[,1])) {\n";
  print $FW "  x = c(data[i,2],data[i,3],data[i,4],data[i,5]);\n";
  print $FW "  alle <- matrix(x,nrow=2);\n";
  print $FW "  #print(alle);\n";
  print $FW "  out<-fisher.test(alle,alternative =\"two.sided\");\n";  #-- it depends
  print $FW "  #write.table(c(data[i,1],data[i,2],data[i,3],data[i,4],data[i,5],out\$p.value),file=\"$out1\",append=TRUE,sep=\"\t\",row.name=FALSE,col.name=FALSE);\n";
  print $FW "  cat(data[i,1],\"\\t\",data[i,2],\"\\t\",data[i,3],\"\\t\",data[i,4],\"\\t\",data[i,5],\"\\t\",out\$p.value,\"\\n\",sep=\"\",file=\"$out1\",append=TRUE);\n";
  print $FW "}\n";
  close($FW);
  system("R --vanilla --slave < $outR");

  #system("paste $out1 $in2t > $out");
  $out = $out1;

  my $outs = $out.".sort";
  open(my $FR,$out);
  open(my $FW1,">$outs");
  my @input1 = map{[$_,(split(/\t/))[5]]} <$FR>;
  my @input2 = sort{$a->[1] <=> $b->[1]} @input1;
  my @sorted = map{$_->[0]} @input2;
  for(@sorted) {
    print $FW1 "$_";
  }
  close($FR);
  close($FW1);

  system("rm $in1t");
  #system("rm $in2t");
  #system("rm $out1");
  system("rm $outR");
  system("rm $out");
}

