#!/usr/bin/perl
use strict;
use List::Util qw(sum shuffle);
use Statistics::R;

my %profile;
my %otus;
my %abundance = &GenerateMock( 10000, 100000000, "lnorm", "mock.txt" );
for ( my $i = 1; $i <= 100; $i++ ) {
  my %abundance_rs = &RandomSampling( \%abundance, 30000 );
  %{ $profile{"sample$i"} } = %abundance_rs;
  foreach my $otu ( keys %abundance_rs ) {
    $otus{$otu} = 1;
  }
  print "$i\n";
}

my @samples = sort keys %profile;
open( OUT, ">mock_profile_s300.txt" ) || die "#1\n";
print OUT "OTU\t", join( "\t", @samples ), "\n";
foreach my $otu ( keys %otus ) {
  print OUT "OTU$otu";
  foreach my $sample (@samples) {
    print OUT "\t$profile{$sample}{$otu}" if $profile{$sample}{$otu};
    print OUT "\t0"                       if !$profile{$sample}{$otu};
  }
  print OUT "\n";
}
close OUT;

sub GenerateMock() {
  my ( $speciesnum, $seqnum, $type, $out ) = @_;
  my $r_sen = qq `library("mobsim")
	mock<-sim_sad($speciesnum,$seqnum,sad_type="$type",sad_coef=list(meanlog=5,sdlog=2))
	#sim_sad($speciesnum,$seqnum,sad_type="power",sad_coef=list(s=2))
	sink("$out")
	mock
	sink()
  `;
  open( R, ">r.tmp" ) || die "#1\n";
  print R "$r_sen";
  close R;
  my $R = "R --vanilla --slave <r.tmp >r.out 2> r.out.tmp";
  system("$R");
  my @abundance;
  open( FILE, "$out" ) || die "#1\n";

  while (<FILE>) {
    if (/^\s+\d+/) {
      my @items = split( /\s+/, $_ );
      push( @abundance, @items[ 1 .. $#items ] );
    }
  }
  close(FILE);
  my %abundance;
  for ( my $i = 1; $i <= $#abundance + 1; $i++ ) {
    $abundance{$i} = $abundance[ $i - 1 ];
  }
  return %abundance;
}

sub RandomSampling() {
  my ( $abundance, $rs ) = @_;
  my %abundance = %$abundance;
  my $sum;
  foreach my $id ( keys %abundance ) {
    $sum += $abundance{$id};
  }
  my %abundance_rs;
  my $R = Statistics::R->new();
  $R->set( 'x',  $sum );
  $R->set( 'rs', $rs );
  $R->run(q'sample_for_perl = sample.int(x, rs)');
  my $Rsample = $R->get('sample_for_perl');
  $R->stop();
  my @array = sort { $a <=> $b } @$Rsample;
  my $i = 1;

  foreach my $gene ( sort { $a <=> $b } keys %abundance ) {
    my @tmp = grep { $_ >= $i && $_ <= ( $abundance{$gene} + $i - 1 ) } @array;
    $abundance_rs{$gene} = @tmp if $#tmp >= 0;
    splice( @array, 0, @tmp );
    $i += $abundance{$gene};
  }
  return %abundance_rs;
}

