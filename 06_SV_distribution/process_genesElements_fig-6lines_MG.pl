#!/usr/bin/perl
#Author: Feng Cheng
use warnings;
use strict;
#require("/home/chengf/bin/perl_subs.pl");
use SVG;
#use GD;
#use GD::Arrow;
#use GD::SVG;

my $in1 = $ARGV[0];
my $in2 = $ARGV[1];
my $in3 = $ARGV[2]; #- the same as $ARGV[0] and $ARGV[1];
my $in4 = $ARGV[3];
my $in5 = $ARGV[4];
my $in6 = $ARGV[5];
my $pic = $ARGV[6].".svg";
my $chr = "ALL";

my ($widtn,$height) = (1000,600);
&process($widtn,$height,$in1,$in2,$in3,$in4,$in5,$in6,$pic,$chr);

sub process {
  my ($widtn,$height,$in1,$in2,$in3,$in4,$in5,$in6,$pic,$chr) = @_;
  my ($xw,$yh,$x0,$y0,$im,$x1,$y1,$x2,$y2,@color,$maxx,$maxy);
  &setCanvas($widtn,$height,\$xw,\$yh,\$x0,\$y0,\$im);
  &initialColor(\@color);
  &drawAxis($xw,$yh,$x0,$y0,$x1,$y1,$x2,$y2,\$im);
  #&drawAxisArrow($xw,$yh,$x0,$y0,$x1,$y1,$x2,$y2,\$im);

  &mainDraw($in1,$in2,$in3,$in4,$in5,$in6,$chr,$maxx,$maxy,$xw,$yh,$x0,$y0,$x1,$y1,$x2,$y2,\@color,\$im);

  open(PNG,">$pic");
  print PNG $im->xmlify;
  close(PNG);
  #system("/home/chengf/bin/distributing_svg_4.74/svg2xxx_release/svg2xxx $pic -t png");
}

sub mainDraw {
  my ($in1,$in2,$in3,$in4,$in5,$in6,$chr,$maxx,$maxy,$xw,$yh,$x0,$y0,$x1,$y1,$x2,$y2,$color,$im) = @_;
  my ($x,$y);

  my ($minx,$miny) = (10000,10000);
  ($maxx,$maxy) = (0,0);
  my %chr2len;
  my (%index2pos1,%index2y1);
  &readMeasure($in1,\%chr2len,\%index2pos1,\%index2y1,\$minx,\$miny,\$maxx,\$maxy,$chr);
  my (%index2pos2,%index2y2);
  &readMeasure($in2,\%chr2len,\%index2pos2,\%index2y2,\$minx,\$miny,\$maxx,\$maxy,$chr);
  my (%index2pos3,%index2y3);
  &readMeasure($in3,\%chr2len,\%index2pos3,\%index2y3,\$minx,\$miny,\$maxx,\$maxy,$chr);

  my (%index2pos4,%index2y4);
  &readMeasure($in4,\%chr2len,\%index2pos4,\%index2y4,\$minx,\$miny,\$maxx,\$maxy,$chr);
  my (%index2pos5,%index2y5);
  &readMeasure($in5,\%chr2len,\%index2pos5,\%index2y5,\$minx,\$miny,\$maxx,\$maxy,$chr);
  my (%index2pos6,%index2y6);
  &readMeasure($in6,\%chr2len,\%index2pos6,\%index2y6,\$minx,\$miny,\$maxx,\$maxy,$chr);
  $maxy = 0.03;  # it depends
  $miny = 0;

  my (%chr2y);
  &drawXYTicks($x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$minx,$miny,$maxx,$maxy,$im,\%chr2y);

  &drawMeasure($x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$minx,$miny,$maxx,$maxy,$im,\%index2pos1,\%index2y1,"red",0);
  &drawMeasure($x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$minx,$miny,$maxx,$maxy,$im,\%index2pos2,\%index2y2,"green",0);
  &drawMeasure($x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$minx,$miny,$maxx,$maxy,$im,\%index2pos3,\%index2y3,"blue",0);

  &drawMeasure($x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$minx,$miny,$maxx,$maxy,$im,\%index2pos4,\%index2y4,"red",5);
  &drawMeasure($x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$minx,$miny,$maxx,$maxy,$im,\%index2pos5,\%index2y5,"green",5);
  &drawMeasure($x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$minx,$miny,$maxx,$maxy,$im,\%index2pos6,\%index2y6,"blue",5);
}

sub readGene {
  my ($in,$gene2chr,$gene2start,$gene2stop) = @_;
  open(my $SFR,$in);
  while($_=<$SFR>) {
    chomp;
    my @temp = split(/\t/);
    $gene2chr->{$temp[0]} = $temp[1];
    $gene2start->{$temp[0]} = $temp[2];
    $gene2stop->{$temp[0]} = $temp[3];
  }
  close($SFR);
}

sub readMeasure {
  my ($in,$chr2len,$index2pos,$index2y,$minx,$miny,$maxx,$maxy,$chr) = @_;

  open(FR,$in);
  my ($index,$pchr) = (1,"A01");
  while($_=<FR>) {
    chomp;
    my @temp = split(/\t/);
    if($temp[0] eq $chr) {
      if($temp[2] eq "NA") {
        $temp[2] = 0;
      }
      if($temp[2] ne "NaN" && $temp[2] ne "NA") {
      $index2pos->{$index} = $temp[1];
        if($$maxx<$index2pos->{$index}) {
          $$maxx = $index2pos->{$index};
        }
        if($$minx>$index2pos->{$index}) {
          $$minx = $index2pos->{$index};
        }
        $index2y->{$index} = $temp[2];
        if($$maxy<$index2y->{$index}) {
          $$maxy = $index2y->{$index};
        }
        if($$miny>$index2y->{$index}) {
          $$miny = $index2y->{$index};
        }
        $index += 1;
        $pchr = $temp[0];
      }
    }
  }
  $chr2len->{$pchr} = $index2pos->{$index-1};
  close(FR);
}

sub drawMeasure {
  my ($x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$minx,$miny,$maxx,$maxy,$im,$index2pos,$index2y,$color,$num) = @_;
  my @tempk = keys %$index2pos;
  my ($maxxl,$maxyl) = ($maxx-$minx,$maxy-$miny);
  my $index = @tempk;
  my ($indext,$ratioL) = (0,0.3);
  my ($lines,$x3,$y3,$x4,$y4,$px1,$py1);
  for(my $i=1;$i<=$index;$i+=1) {
    $indext += 1;
    ($x,$y) = ($index2pos->{$i}-$minx,$index2y->{$i}-$miny);
    &transRatio(\$x,\$y,$xw,$yh,$maxxl,$maxyl);
    ($x1,$y1) = ($x,$y);
    ($x2,$y2) = (0,0);
    &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
    if($indext==1) {
      $lines = "M$x1 $y1";
    }
    else {
      my $tempL = $x1 - $px1;
      ($x3,$y3) = ($tempL*$ratioL + $px1,$py1);
      ($x4,$y4) = ($x1 - $tempL*$ratioL,$y1);
      $lines .= ", C$x3 $y3 $x4 $y4 $x1 $y1";
    }
    ($px1,$py1) = ($x1,$y1);

    #($x,$y) = ($index2pos->{$i},$index2y->{$i});
    #&transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
    #($x1,$y1) = ($x,$y);
    #&transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
    #$$im->circle(cx=>$x1,cy=>$y1,r=>0.3,stroke=>'red','stroke-width'=>0.3,fill=>'red');
  }
  if($num==0) {
    $$im->path(
      d=>"$lines",
      style=>{
        stroke=>$color,
        fill=>"none"
      }
    );
  }
  else {
    $$im->path(
      d=>"$lines",
      style=>{
        "stroke-dasharray"=>"$num,$num",
        stroke=>$color,
        fill=>"none"
      }
    );
  }
}

sub drawXYTicks {
  my ($x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$minx,$miny,$maxx,$maxy,$im,$chr2y) = @_;
  my $flag = 1;
  for(my $i=int($minx/1000);$i*1000<=$maxx;$i+=1) {
    ($x,$y) = ($i*1000-$minx,0-$miny);
    &transRatio(\$x,\$y,$xw,$yh,$maxx-$minx,$maxy-$miny);
    ($x1,$y1,$x2,$y2) = ($x,0,$x,-10);
    &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
    $$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'black','stroke-width'=>1,fill=>'black');
    if($flag==1) {
      ($x1,$y1) = ($x-20,-10-15);
      &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
      $$im->text(
        x=>$x1, y=>$y1,
        style => {
          'font-family' => 'Arial',
          'font-size' => '9pt'
          #'stroke' => 'black',
          #'font-weight' => 'bold'
        }
      )->cdata($i." Kb");
    }
    if($flag==1) {
      $flag = 0;
    }
    else {
      $flag = 1;
    }
  }
  ($x1,$y1,$x2,$y2) = ($xw,0,$xw,$yh);
  &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
  $$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'black','stroke-width'=>1,fill=>'black');

  $flag = 0;
  #for(my $i=int($miny*10)/10;$i<=(int($maxy*10)+1)/10;$i+=((int($maxy*10)+0.2)/10)/10) {
  my $tn = 1;
  if($maxy=~/0\.(0*)[^0]/) {
    $tn = length($1)+1;
  }
  for(my $i=int($miny*(10**$tn))/(10**$tn);$i<=(int($maxy*(10**$tn))+1)/(10**$tn);$i+=((int($maxy*(10**$tn))+0.2)/(10**$tn))/10) {
    ($x,$y) = (0-$minx,$i-$miny);
    &transRatio(\$x,\$y,$xw,$yh,$maxx-$minx,$maxy-$miny);
    ($x1,$y1,$x2,$y2) = (-10,$y,0,$y);
    &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
    $$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'black','stroke-width'=>1,fill=>'black');

    if($flag==1) {
      if($i==0) {
        ($x1,$y1) = (-10-8,$y-4);
      }
      else {
        ($x1,$y1) = (-10-8,$y-14);
      }
      &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
      $$im->text(
        x=>$x1, y=>$y1,
        style => {
        'font-family'  => 'Arial',
        #'stroke' => 'black',
        'font-size' => '9pt'
        #'font-weight' => 'bold'
        },
        transform => "rotate(-90,$x1,$y1)"
      )->cdata($i."");
    }
    if($flag==1) {
      $flag = 0;
    }
    else {
      $flag = 1;
    }
  }

  ($x,$y) = (0,int($maxy/2));
  &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
  ($x1,$y1) = (-10-40,$y+106);
  &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
  $$im->text(
    x=>$x1, y=>$y1,
    style => {
    'font-family'  => 'Arial',
    #'stroke' => 'black',
    'font-size' => '10pt'
    #'font-weight' => 'bold'
    },
   transform => "rotate(-90,$x1,$y1)"
  )->cdata("Average number of xx bp in 100 bp window");

  ($x1,$y1,$x2,$y2) = (0,$yh,$xw,$yh);
  &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
  $$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'black','stroke-width'=>1,fill=>'black');
}

sub drawCentromere {
  my ($x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$maxx,$maxy,$im,$gene2start,$gene2stop,$chr2y) = @_;
  my $cent = "centromere";
  open(CE,$cent);
  while($_=<CE>) {
    if(/^A/) {
      chomp;
      my @temp = split(/\s+/);
      ($x,$y) = ($gene2start->{$temp[3]},$chr2y->{$temp[0]});
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x1,$y1) = ($x,$y);
      ($x,$y) = ($gene2start->{$temp[3]},$chr2y->{$temp[0]}+0.6);
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x2,$y2) = ($x,$y);
      &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
      my ($xp1,$yp1,$xp2,$yp2) = ($x1,$y1,$x2,$y2);
      ($x,$y) = ($gene2stop->{$temp[4]},$chr2y->{$temp[0]});
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x1,$y1) = ($x,$y);
      ($x,$y) = ($gene2stop->{$temp[4]},$chr2y->{$temp[0]}+0.6);
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x2,$y2) = ($x,$y);
      &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
      my ($xp3,$yp3,$xp4,$yp4) = ($x2,$y2,$x1,$y1);

      my $xv = [$xp1,$xp2,$xp3,$xp4];
      my $yv = [$yp1,$yp2,$yp3,$yp4];
      my $points = $$im->get_path(
        x=>$xv, y=>$yv,
        -type=>'polygon'
      );
      $$im->polygon(
        %$points,
        fill => 'grey'
      );
    }
  }
  close(CE);
}

sub drawCandidateRegion {
  my ($in3,$x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$maxx,$maxy,$im,$chr2y) = @_;

  open(my $SFR3,$in3);
  while($_=<$SFR3>) {
    if(/^A/) {
      chomp;
      my @temp = split(/\s+/);
      ($x,$y) = ($temp[1],$chr2y->{$temp[0]});
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x1,$y1) = ($x,$y);
      ($x,$y) = ($temp[1],$chr2y->{$temp[0]}+0.6);
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x2,$y2) = ($x,$y);
      &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
      my ($xp1,$yp1,$xp2,$yp2) = ($x1,$y1,$x2,$y2);
      ($x,$y) = ($temp[2],$chr2y->{$temp[0]});
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x1,$y1) = ($x,$y);
      ($x,$y) = ($temp[2]+5000,$chr2y->{$temp[0]}+0.6);
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x2,$y2) = ($x,$y);
      &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
      my ($xp3,$yp3,$xp4,$yp4) = ($x2,$y2,$x1,$y1);

      my $xv = [$xp1,$xp2,$xp3,$xp4];
      my $yv = [$yp1,$yp2,$yp3,$yp4];
      my $points = $$im->get_path(
        x=>$xv, y=>$yv,
        -type=>'polygon'
      );
      $$im->polygon(
        %$points,
        'stroke-width' => 1,
        stroke => 'orange',
        fill => 'orange'
      );
    }
  }
  close($SFR3);
}

sub drawCandidateGenes {
  my ($ings,$x1,$y1,$x2,$y2,$xw,$yh,$x0,$y0,$x,$y,$maxx,$maxy,$im,$gene2chr,$gene2start,$gene2stop,$chr2y) = @_;
  open(my $SFR,$ings);
  while($_=<$SFR>) {
    chomp;
    my @temp = split(/\t/);
    my $gene = $temp[0];
    #print "$gene\t$gene2chr{$gene}\t$gene2start{$gene}\t$gene2stop{$gene}\n";
    if($gene2chr->{$gene}=~/^A/) {
      my $pos = ($gene2start->{$gene}+$gene2stop->{$gene})/2;
      ($x,$y) = ($pos,$chr2y->{$gene2chr->{$gene}});
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x1,$y1) = ($x,$y-2);
      ($x,$y) = ($pos,$chr2y->{$gene2chr->{$gene}});
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x2,$y2) = ($x-2,$y-8);
      &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
      my ($xp1,$yp1,$xp2,$yp2) = ($x1,$y1,$x2,$y2);
      #$$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'red','stroke-width'=>1,fill=>'red');
      ($x,$y) = ($pos,$chr2y->{$gene2chr->{$gene}});
      &transRatio(\$x,\$y,$xw,$yh,$maxx,$maxy);
      ($x2,$y2) = ($x+2,$y-8);
      &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
      my ($xp3,$yp3) = ($x2,$y2);

      my $xv = [$xp1,$xp2,$xp3];
      my $yv = [$yp1,$yp2,$yp3];
      my $points = $$im->get_path(
        x=>$xv, y=>$yv,
        -type=>'polygon'
      );
      $$im->polygon(
        %$points,
        fill => 'red'
      );
    }
  }
  close($SFR);
}

sub setCanvas {
  my ($widtn,$height,$xw,$yh,$x0,$y0,$im) = @_;
  $$im = SVG->new(width=>$widtn,height=>$height);
  ($$xw,$$yh,$$x0,$$y0) = ($widtn-200,$height-200,100,$height-100);
}

sub drawAxis {
  my ($xw,$yh,$x0,$y0,$x1,$y1,$x2,$y2,$im) = @_;
  my ($x,$y);
  ($x1,$y1,$x2,$y2) = (0,0,$xw,0);
  &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
  $$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'black','stroke-width'=>0.5,fill=>'black');
  ($x1,$y1,$x2,$y2) = (0,0,0,$yh);
  &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
  $$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'black','stroke-width'=>0.5,fill=>'black');
}

sub drawAxisArrow {
  my ($xw,$yh,$x0,$y0,$x1,$y1,$x2,$y2,$im) = @_;
  #arrows----
  my $lenArrow = 8;
  my $widthArrow = $lenArrow/3;
  ($x1,$y1,$x2,$y2) = ($xw-$lenArrow,0+$widthArrow,$xw,0);
  &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
  $$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'black','stroke-width'=>1,fill=>'black');
  ($x1,$y1,$x2,$y2) = ($xw-$lenArrow,0-$widthArrow,$xw,0);
  &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
  $$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'black','stroke-width'=>1,fill=>'black');
  ($x1,$y1,$x2,$y2) = (0+$widthArrow,$yh-$lenArrow,0,$yh);
  &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
  $$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'black','stroke-width'=>1,fill=>'black');
  ($x1,$y1,$x2,$y2) = (0-$widthArrow,$yh-$lenArrow,0,$yh);
  &transPos(\$x1,\$y1,\$x2,\$y2,$x0,$y0);
  $$im->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,stroke=>'black','stroke-width'=>1,fill=>'black');
}

sub transPos {
  my ($x1,$y1,$x2,$y2,$x0,$y0) = @_;
  $$x1 += $x0;
  $$y1 = $y0 - $$y1;
  $$x2 += $x0;
  $$y2 = $y0 - $$y2;
}

sub transRatio {
  my ($x,$y,$xw,$yh,$maxx,$maxy) = @_;
  $$x /= $maxx;
  $$x *= $xw;
  $$y /= $maxy;
  $$y *= $yh;
}

sub initialColor {
  my ($color) = @_;
  $color->[0] = "#".sprintf("%02x",255).sprintf("%02x",99).sprintf("%02x",71);
  $color->[1] = "#".sprintf("%02x",100).sprintf("%02x",149).sprintf("%02x",237);
  $color->[2] = "#".sprintf("%02x",50).sprintf("%02x",205).sprintf("%02x",50);
  $color->[3] = "#".sprintf("%02x",255).sprintf("%02x",255).sprintf("%02x",0);
  $color->[4] = "#".sprintf("%02x",255).sprintf("%02x",105).sprintf("%02x",180);
  $color->[5] = "#".sprintf("%02x",0).sprintf("%02x",191).sprintf("%02x",255);
  $color->[6] = "#".sprintf("%02x",255).sprintf("%02x",228).sprintf("%02x",225);
  $color->[7] = "#".sprintf("%02x",192).sprintf("%02x",192).sprintf("%02x",192);
  $color->[8] = "#".sprintf("%02x",220).sprintf("%02x",20).sprintf("%02x",60);
  $color->[9] = "#".sprintf("%02x",0).sprintf("%02x",0).sprintf("%02x",139);
  $color->[10] = "#".sprintf("%02x",0).sprintf("%02x",128).sprintf("%02x",0);
  $color->[11] = "#".sprintf("%02x",255).sprintf("%02x",215).sprintf("%02x",0);
  $color->[12] = "#".sprintf("%02x",255).sprintf("%02x",0).sprintf("%02x",255);
  $color->[13] = "#".sprintf("%02x",30).sprintf("%02x",144).sprintf("%02x",255);
  $color->[14] = "#".sprintf("%02x",105).sprintf("%02x",105).sprintf("%02x",105);
  $color->[15] = "#".sprintf("%02x",250).sprintf("%02x",128).sprintf("%02x",114);
  $color->[16] = "#".sprintf("%02x",65).sprintf("%02x",105).sprintf("%02x",255);
  $color->[17] = "#".sprintf("%02x",46).sprintf("%02x",139).sprintf("%02x",87);
  $color->[18] = "#".sprintf("%02x",255).sprintf("%02x",165).sprintf("%02x",0);
  $color->[19] = "#".sprintf("%02x",255).sprintf("%02x",20).sprintf("%02x",147);
  $color->[20] = "#".sprintf("%02x",0).sprintf("%02x",0).sprintf("%02x",128);
  $color->[21] = "#".sprintf("%02x",128).sprintf("%02x",128).sprintf("%02x",128);
}


