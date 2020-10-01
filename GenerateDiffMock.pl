#!/usr/bin/perl
use strict;
use List::Util qw(shuffle);

my $DiffOTUPer=0.1;###percentage of differece at OTU level
my $num=12;
my $seedmock="mock.txt";

my %data;
my @otus;
open(SEED,$seedmock)||die"#cannot open $seedmock";
while(<SEED>){
  chomp;
  my @items=split("\t",$_);
  $data{$items[0]}=$items[1];
  push(@otus,$items[0]);
}
close SEED;

for(my $i=1;$i<=$num;$i++){
  @shuffotus=shuffle @otus;
  my $end=$DiffOTUPer*$#otus;
  my @diffotus=shuffle @shuffotus[0..$end];
  my %otumap;
  for(my $i=0;$i<=$end;$i++){
    $otumap{$shuffotus[$i]}=$diffotus[$i];
  }
  open(OUT,">mock$i.txt")||die"#2 cannot write $mock$i.txt";
  foreach my $otu(@otus){
    if(!$otumap{$otu}){
      print OUT "$otu\t$data{$otu}\n";
    }
    else{
      print OUT "$otu\t$data{$otumap{$otu}}\n";
    }
  }
  close OUT;
}
