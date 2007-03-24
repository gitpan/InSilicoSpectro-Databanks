#!/usr/bin/env perl
use strict;
use Carp;
use Pod::Usage;

=head1 NAME

databanks-decoy.pl - decoy input databanks following several moethods

=head1 DESCRIPTION

Reads input fasta file and produce a decoyed databanks with several methods:

=over 4

=item reverse: simply reverse each the sequence

=item shuffle: shuffle AA in each sequence

=item shuffle & avoid known cleaved peptides: shuffe sequence but avoid producing kown trayptic peptides

=item Markov model: learn Markov model chain distribution of a given level), then produces entries corresponding to this distribution

=back

=head1 SYNOPSIS

#reverse sequences for a local (optionaly compressed) file
databanks-decoys.pl --in=/tmp/uniprot_sprot.fasta.gz --method=reverse

#download databanks from the web | uncompress it and shuffle the sequence
wget -silent -O - ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz | zcat |  databatanks-decoy.pl --method=shuffle

#use a .dat file (with splice forms) as an input
uniprotdat2fasta.pl --in=uniprot_sprot_human.dat | databanks-decoy.pl --method=markovmodel


=head1 ARGUMENTS


=head3 --in=infile.fasta

An input fasta file (will be uncompressed if ending with gz)

=head3 -out=outfile.fasta

A .fasta file [default is stdout]

=head3 --method=(reverse|shuffle|markovmodel)

Set the decoying method

=head1 OPTIONS

=head2 --method=shuffle options

=head3 --shuffle-ejectcleavedpeptides

Re-shuffle peptides of size >=6 that where detected as cleaved one in original databank

=head3 --shuffle-cleaveenzyme=regexp

Set a regular expression for the enzyme [default is trysin: '(?<=[KR])(?=[^P])']

=head3 --shuffle-testenzyme

Just digest entries with the set enzyme and produces space separated peptides (to check the enzyme)

=head2 --method=markovmodel options

=head3 --markovmodel-level=int [default 3]

Set length of the model (0 means only AA distrbution will be respected, 3 means chains of length 3  distribution etc.). Setting a length >3 can deal to memory burnout.

=head2 misc

=head3 --noprogressbar

do not display terminal progress bar (if possible)

=head3 --help

=head3 --man

=head3 --verbose

Setting an environment vartiable DO_NOT_DELETE_TEMP=1 will keep the temporay file after the script exit

=head1 EXAMPLE


=head1 COPYRIGHT

Copyright (C) 2004-2006  Geneva Bioinformatics www.genebio.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

=head1 AUTHORS

Alexandre Masselot, www.genebio.com

=cut



use SelectSaver;
use File::Basename;
use List::Util qw/shuffle/;

use Getopt::Long;
my ($method, $noProgressBar);

my ($shuffle_ejectCleavPept, $shuffle_test);
my $shuffle_enzyme='(?<=[KR])(?=[^P])';

my $inFile='-';
my $outFile='-';
my $inFD;

my($help, $man, $verbose);

if (!GetOptions(
		"in=s"=>\$inFile,
		"out=s" => \$outFile,

		"method=s" => \$method,

		"shuffle-ejectcleavedpeptides"=> \$shuffle_ejectCleavPept,
		"shuffle-testenzyme"=> \$shuffle_test,
		"shuffle-cleaveenzyme"=> \$shuffle_enzyme,

		"noprogressbar" => \$noProgressBar,

                "help" => \$help,
                "man" => \$man,
                "verbose" => \$verbose,
               )
    || $help || $man){


  pod2usage(-verbose=>2, -exitval=>2) if(defined $man);
  pod2usage(-verbose=>1, -exitval=>2);
}

die "no --method=methodname argument (see --help)" unless $method;

#init parsing progress bar
my ($pg, $size, $nextpgupdate, $readsize);
$readsize=0;
$nextpgupdate=0;
eval{
  if ((!$noProgressBar) && ($inFile != /gz$/i) && ($inFile ne '-')&& -t STDIN && -t STDOUT){
    require Term::ProgressBar;
    $size=(stat $inFile)[7];
    $pg=Term::ProgressBar->new ({name=> "parsing ".basename($inFile),
				 count=>$size,
				 ETA=>'linear',
				 remove=>1
				});
  }
};
if ($@){
  warn "could not use Term::ProgressBar (consider installing the module for more interactive use";
}

__setInput();

#set output on default
my $saver;
if ($outFile ne '-'){
  if($outFile=~/\.gz$/i){
    open FDOUT, ">:gzip", $outFile or die "cannot open for writing gziped [$outFile]: $!";
    $saver=new SelectSaver(\*FDOUT);
  }else{
    open (FDOUT, ">$outFile") or die "could not open for writing [$outFile]:$!";
    $saver=new SelectSaver(\*FDOUT);
  }
}

if($method eq 'reverse'){
  my $nbentries=0;
  while((my ($head, $seq)=__nextEntry())[0]){
    $nbentries++;
    $head=~s/>/>REV_/ or die "entry header does not start with '>': $head";
    $seq=reverse $seq;
    print  "$head\n";
    print __prettySeq($seq)."\n";
  }
  print STDERR "reverted $nbentries\n" if $verbose;
}elsif($method eq 'shuffle'){
  unless ($shuffle_ejectCleavPept){
    my $nbentries=0;
    while((my ($head, $seq)=__nextEntry())[0]){
      $nbentries++;
      $head=~s/>/>SHFL_/ or die "entry header does not start with '>': $head";
      $seq=join ('', shuffle  split(//, $seq));
      print  "$head\n";
      print __prettySeq($seq)."\n";
    }
    print STDERR "shuffled $nbentries\n" if $verbose;
  }else{
    my $enz=qr/$shuffle_enzyme/;
    my %pepts;
    my $nbentries=0;
    while((my ($head, $seq)=__nextEntry())[0]){
      $nbentries++;
      foreach (split $enz, $seq) {
	$pepts{$_}++ if length($_)>=6;
      }
    }
    __setInput();
    my @histoEjected;
    my $nreshuffled=0;
    while((my ($head, $seq)=__nextEntry())[0]){
      $head=~s/>/>SHFLEJ_/ or die "entry header does not start with '>': $head";
      $seq=join ('', shuffle (split //, $seq));
      my $newseq="";
      my $nreshuffperseq=0;
      foreach (split $enz, $seq) {
	if (length($_)<6) {
	  $newseq.=$_;
	  next;
	}
	while ($pepts{$_}) {
	  $_=join ('', shuffle (split //, $_));
	  $nreshuffled++;
	  $nreshuffperseq++;
	}
	$newseq.=$_;
      }
      print "$head\n".__prettySeq($newseq)."\n";
      $histoEjected[$nreshuffperseq]++;
    }
    if($verbose){
      print STDERR "reshuffled pept/nb sequences: $nreshuffled/$nbentries\n";
      print STDERR "#seq\t#nb reshuffled peptides\n";
      foreach (0..$#histoEjected) {
	next unless $histoEjected[$_];
	print STDERR "$_\t$histoEjected[$_]\n";
      }
    }
  }
}else{
  die "unimplemented method [$method]";
}


sub __nextEntry{
  local $/="\n>";
  my $contents=<$inFD>;
  return undef unless $contents;

  $readsize+=length $contents;
  $nextpgupdate=$pg->update($readsize) if $pg && $readsize>=$nextpgupdate;

  chomp $contents;
  $contents=">$contents" unless $contents=~/^>/;
  my ($head, $seq)=split /\n/,$contents,2;
  $seq=~s/\s+//g;
  return ($head, $seq);
}

sub __prettySeq{
  $_=$_[0];
  s/(\S{60})(?=\S)/$1\n/g;
  return $_;
}

sub __setInput{
  #set input
  if ($inFile eq '-'){
    $inFD=\*STDIN;
  }elsif($inFile=~/\.gz$/i){
    require PerlIO::gzip;
    open $inFD, "<:gzip", "$inFile" or die $!;
  }else{
    open ($inFD, "<$inFile") or die "cannot open for reading [$inFile]: $!";
  }
}

