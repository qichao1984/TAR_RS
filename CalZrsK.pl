#!/usr/bin/perl
use strict;
use List::Util qw(sum shuffle);
use Statistics::R;
use Statistics::LineFit;

my @df_num = ( 0, 500, 1000, 1500, 2000, 2500, 3000 );#increased species number per increased sampling area
my @rs     = ( 20000 );##sequencing depth
#my @rs     = ( 10000, 20000, 30000, 40000, 50000 );
my @radius = ( 1, 10, 50, 100, 200 );
my $samplenum=4;#number of samples per increased sampling area

my %z;
foreach my $rs (@rs) {
  foreach my $df_num (@df_num) {
	print "Working on sequencing depth $rs with $df_num new species per increased sampling area...\n";
    my $workdir = "S$rs";
    my @mock = ( 1 .. @radius );
    mkdir $workdir if !-e $workdir;
	open(OUT1,">$workdir/slope_report_S$rs\_$df_num.txt")||die"#1 can not open $workdir/slope_report.txt\n";
	print OUT1 "Zexp\tZobs\n";
    system("cp mock1.txt $workdir/mock1_df$df_num.txt");
    for ( my $i = 1; $i <= 4; $i++ ) {
      my $j         = $i + 1;
      my $df        = $df_num;
      my %abundance = &ReadTabMock("$workdir/mock$i\_df$df_num.txt");
      my @otus = shuffle keys %abundance;
      my %diffotus;
      foreach my $otu ( @otus[ 0 .. $df ] ) {
        $diffotus{$otu} = 1;
      }
      my $sum;
      open( OUT, ">$workdir/mock$j\_df$df_num.txt" ) || die "#2 can not open $workdir/mock$j\_df$df_num.txt\n";
      foreach my $otu ( @otus[ 0 .. 9998 ] ) {
        print OUT "$otu\t$abundance{$otu}\n"     if !$diffotus{$otu};
        print OUT "$otu\_$j\t$abundance{$otu}\n" if $diffotus{$otu};
        $sum += $abundance{$otu};
      }
      close OUT;
    }
    foreach my $mock (@mock) {
      my %profile;
      my %otus;
      my %abundance = &ReadTabMock("$workdir/mock$mock\_df$df_num.txt");
      for ( my $i = 1; $i <= 30; $i++ ) {
        my %abundance_rs = &RandomSampling( \%abundance, $rs );
        %{ $profile{"sample$i"} } = %abundance_rs;
        foreach my $otu ( keys %abundance_rs ) {
          $otus{$otu} = 1;
        }
        print "$workdir\tdf$df_num\t$mock\t$i\n";
      }
      my @samples = sort keys %profile;
      open( OUT, ">$workdir/mock$mock\_df$df_num\_profile_s30.txt" ) || die "#4 can not open $workdir/mock$mock\_df$df_num\_profile_s30.txt\n";
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
    }
    my @exp;
    my @files;
    foreach my $mock(@mock){
      push(@files,"mock$mock\_df$df_num.txt");
    }
    for ( my $i = 0; $i <= $#mock; $i++ ) {
      my %otus;
      foreach my $file ( @files[ 0 .. $i ] ) {
      	open( FILE, "$workdir/$file" ) || die "#5 can not open $workdir/$file\n";
        while (<FILE>) {
          my @items = split( "\t", $_ );
          $otus{ $items[0] } = 1;
        }
        close FILE;
      }
      push( @exp, scalar( keys %otus ) );
    }
    my @obs;
    my %otus;
    for ( my $num = 1; $num <= @mock; $num++ ) {
      my $file = "$workdir/mock$num\_df$df_num\_profile_s30.txt";
      my %data;
      open( FILE, "$file" ) || die "#6 can not open $file\n";
      my $line = <FILE>;
      chomp $line;
      my @heads = split( "\t", $line );
      while (<FILE>) {
        chomp;
        my @items = split( "\t", $_ );
        for ( my $i = 1; $i <= $#items; $i++ ) {
          $data{ $heads[$i] }{ $items[0] } = 1 if $items[$i] > 0;
        }
      }
      close FILE;
      my @samples = shuffle @heads[ 1 .. $#heads ];
      @samples = shuffle @samples;
      my @subsamples = @samples[ 0 .. ($samplenum-1) ];
      foreach my $sample (@subsamples) {
        my @otus = keys %{ $data{$sample} };
        foreach my $otu (@otus) {
          $otus{$num}{$otu} = 1;
        }
      }
    }
    foreach my $num1 ( sort keys %otus ) {
      my %temp;
      foreach my $num2 ( keys %otus ) {
        if ( $num2 <= $num1 ) {
          my @otus = keys %{ $otus{$num2} };
          foreach my $otu (@otus) {
            $temp{$otu} = 1;
          }
        }
      }
      print "area$num1\tspecies num:", scalar( keys %temp ), "\n";
      push( @obs, scalar( keys %temp ) );
    }
    my ( @logarea, @logexp, @logobs );
    for ( my $i = 0; $i <= $#radius; $i++ ) {
      push( @logarea, log( 3.14 * ( $radius[$i]**2 ) ) );
      push( @logexp,  log( $exp[$i] ) );
      push( @logobs,  log( $obs[$i] ) );
    }
    my $lineFit = Statistics::LineFit->new();
    $lineFit->setData( \@logarea, \@logexp ) or die "Invalid data";
    my ( $intercept, $slope ) = $lineFit->coefficients();
    print OUT1 "df$df_num\t$slope\t";
    $z{$rs}{$df_num}{"exp"}=$slope;
    $lineFit->setData( \@logarea, \@logobs ) or die "Invalid data";
    ( $intercept, $slope ) = $lineFit->coefficients();
    print OUT1 "$slope\n";
    $z{$rs}{$df_num}{"obs"}=$slope;
  }
}
close OUT1;


print "############BEGIN REPORT############\n";
foreach my $rs(@rs){
  my $zrs;
  if($z{$rs}{"0"}){
    $zrs=$z{$rs}{"0"}{"obs"};
  }
  else{
    die "Please include 0 in df_nums to calculate Zrs\n";
  }
  print "Under $rs sequencing depth:\nThe random sampling noise is Zrs=$zrs\n";
  my @k;
  foreach my $df_num(@df_num[1..$#df_num]){
    my $zobs=$z{$rs}{$df_num}{"obs"};
    my $zexp=$z{$rs}{$df_num}{"exp"};
    my $k=($zobs-$zrs)/$zexp;
	push (@k,$k);
    print "At $df_num new species per new area,the estimated k=$k\n";
  }
  my $k_mean=sum(@k)/@k;
  print "The mean k=$k_mean\n";
}
print "####################################\n";

sub ReadTabMock() {
  my @files = @_;
  my %abundance;
  foreach my $file (@files) {
    open( FILE, "$file" ) || die "#1\n";
    while (<FILE>) {
      chomp;
      my @items = split( "\t", $_ );
      $abundance{ $items[0] } = $items[1];
    }
    close FILE;
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
