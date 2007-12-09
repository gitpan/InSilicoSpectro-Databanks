#!/usr/bin/env perl
use strict;
use Carp;
use Pod::Usage;

=head1 NAME

fasta-decoy.pl - decoy input databanks following several moethods

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
  fasta-decoys.pl --in=/tmp/uniprot_sprot.fasta.gz --method=reverse

  #download databanks from the web | uncompress it and shuffle the sequence
  wget -silent -O - ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz | zcat |  databatanks-decoy.pl --method=shuffle

  #use a .dat file (with splice forms) as an input
  uniprotdat2fasta.pl --in=uniprot_sprot_human.dat | fasta-decoy.pl --method=markovmodel

  #reversing each sequence
  fasta-decoy.pl --ac-prefix=DECOY_ --in=mitoch.fasta --method=reverse --out=mitoch-reverse.fasta

  #drawing amino acid following distribution in original fasta (end of sequence is considered as a learned random event)
  fasta-decoy.pl --ac-prefix=DECOY_ --in=mitoch.fasta --method=markovmodel --markovmodel-level=0 --out=mitoch-markovmodel_0.fasta

  #drawing amino acid with a markov model (here of length 3)
  fasta-decoy.pl --ac-prefix=DECOY_ --in=mitoch.fasta --method=markovmodel --markovmodel-level=3 --out=mitoch-markovmodel_3.fasta

  #each sequence is randomly shuffled
  fasta-decoy.pl --ac-prefix=DECOY_ --in=mitoch.fasta --method=shuffle --out=mitoch-shuffle.fasta

  #idem, but no tryptic peptide (of length>=6) from the original bank must be found in the random one; 
  #crc is the number of bit attributed to builiding the index mechanism, to avoid memory overallocation
  #crc=33 means 1GB of ram will be used to store the index, 32 means 512MB, 34 means 2GB etc..)
  fasta-decoy.pl --ac-prefix=DECOY_ --in=mitoch.fasta --method=shuffle --shuffle-reshufflecleavedpeptides-crc=33 --shuffle-reshufflecleavedpeptides --out=mitoch-shuffleplus.fasta


=head1 ARGUMENTS


=head3 --in=infile.fasta

An input fasta file (will be uncompressed if ending with gz)

=head3 -out=outfile.fasta

A .fasta file [default is stdout]

=head3 --method=(reverse|shuffle|markovmodel)

Set the decoying method

=head1 OPTIONS

=head2 --ac-prefix=string

Set a key to be prepended befor the AC in the randomized banks. By default, it will be dependent on the choosen method)

=head2 --method=shuffle options

=head3 --shuffle-reshufflecleavedpeptides

Re-shuffle peptides of size >=6 that where detected as cleaved one in original databank

=head3 --shuffle-reshufflecleavedpeptides-minlength [default 6]

Set the size of the peptide to be reshuffled is they already exist

=head3 --shuffle-reshufflecleavedpeptides-crc=int

Building a hash of known cleaved peptide can be quite demanding for memory (uniprot_trembl => ~4GB). Thereforea solution is to make an but array containing stating if or not a peptide with corresponding crc code was found.

=head3 --shuffle-cleaveenzyme=regexp

Set a regular expression for the enzyme [default is trypsin: '.*?[KR](?=[^P])|.+$']

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
use String::CRC;
use Bit::Vector;
use List::Util qw/shuffle/;

use Getopt::Long;
my ($acPrefix, $method, $noProgressBar);

my ($shuffle_reshuffleCleavPept, $shuffle_test);
my $shuffle_enzyme='.*?[KR](?=[^P])|.+$';
my $shuffle_reshuffleCleavPept_minLength=6;
my $shuffle_reshuffleCleavPept_CRCLen;

my $markovmodel_level=2;

my $inFile='-';
my $outFile='-';
my $inFD;

my($help, $man, $verbose);

if (!GetOptions(
		"in=s"=>\$inFile,
		"out=s" => \$outFile,

		"method=s" => \$method,
		"ac-prefix=s" =>\$acPrefix,

		"shuffle-reshufflecleavedpeptides"=> \$shuffle_reshuffleCleavPept,
		"shuffle-reshufflecleavedpeptides-minlength=i"=> \$shuffle_reshuffleCleavPept_minLength,
		"shuffle-reshufflecleavedpeptides-crc=i"=> \$shuffle_reshuffleCleavPept_CRCLen,

		"shuffle-testenzyme"=> \$shuffle_test,
		"shuffle-cleaveenzyme"=> \$shuffle_enzyme,

		"markovmodel-level=i"=>\$markovmodel_level,

		"noprogressbar" => \$noProgressBar,

                "help" => \$help,
                "man" => \$man,
                "verbose" => \$verbose,
               )
    || $help || $man){


  pod2usage(-verbose=>2, -exitval=>2) if(defined $man);
  pod2usage(-verbose=>1, -exitval=>2);
}

my ($nbVCRC, $modVCRC, $maxBitCRC);
if($shuffle_reshuffleCleavPept_CRCLen>=32){
  $nbVCRC=1<<$shuffle_reshuffleCleavPept_CRCLen-31;
  $maxBitCRC=1<<31;
}else{
  $nbVCRC=0;
  $maxBitCRC=1<<$shuffle_reshuffleCleavPept_CRCLen;
}
print STDERR "nb CRC vector=$nbVCRC\nmax bit 4 crc=$maxBitCRC\n" if $verbose;

die "no --method=methodname argument (see --help)" unless $method;

#init parsing progress bar
my ($pg, $size, $nextpgupdate, $readsize);
$readsize=0;
$nextpgupdate=0;
my $imaxreshuffle=1000;
my $imaxreshuffleFinal=10000;


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
  $acPrefix ||='REV_';
  my $nbentries=0;
  while((my ($head, $seq)=__nextEntry())[0]){
    $nbentries++;
    $head=~s/>/>$acPrefix/ or die "entry header does not start with '>': $head";
    $seq=reverse $seq;
    print  "$head\n";
    print __prettySeq($seq)."\n";
  }
  print STDERR "reverted $nbentries\n" if $verbose;
} elsif ($method eq 'markovmodel') {
  $acPrefix ||='MARKOVMODEL_';
  my $nbentries=0;
  die "--markovmodel-level=int [$markovmodel_level] should be between 0 and 3" unless $markovmodel_level>=0 && $markovmodel_level<=3;

  my %chain;			#count and prob 
  #counting
  while ((my ($head, $seq)=__nextEntry())[0]) {
    $nbentries++;
    $seq.='*';
    my $buf='';
    while ($seq=~/(.)/g) {
      my $p=pos $seq;
      my $q=$p-$markovmodel_level-1;
      $q=0 if $q<0;
      my $prec=substr $seq, $q, $p-$q-1;
      $chain{$prec || '^'}{$1}++;
    }
  }
  #normalize;
  foreach my $h (values %chain) {
    my $tot=0;
    $tot+=$_ foreach values %$h;
    $h->{$_}=$h->{$_}*1.0/$tot foreach keys %$h;
  }
  #generate:
  eval{
    if ( -t STDIN && -t STDOUT){
      require Term::ProgressBar;
      $size=(stat $inFile)[7];
      $pg=Term::ProgressBar->new ({name=> "markovmodel generation",
				   count=>$nbentries,
				   ETA=>'linear',
				   remove=>1
				  });
    }
    $nextpgupdate=0;
    $readsize=0;
  };


  foreach (1..$nbentries) {
    $readsize++;
    $nextpgupdate=$pg->update($readsize) if $pg && $readsize>=$nextpgupdate;

    print ">$acPrefix$_\n";
    my $seq='';
    my $prec='';
    while (1) {
      my $hNext;
      $hNext=$prec?$chain{$prec}:$chain{'^'};
      my $crand=rand;
      my $nextAA='?';
      foreach $_ (keys %$hNext) {
	$crand-=$hNext->{$_};
	if ($crand<=0) {
	  $nextAA=$_;
	  last;
	}
      }
      last if $nextAA eq '*';
      if ($markovmodel_level>0) {
	$prec.=$nextAA;
	$prec=~s/.+(?=.{$markovmodel_level})//;
      }
      die "NO next AA generated" if $nextAA eq '?';
      $seq.=$nextAA;
    }
    print __prettySeq($seq)."\n";
  }
} elsif ($method eq 'shuffle') {
  unless ($shuffle_reshuffleCleavPept){
    $acPrefix ||='SHFL_';
    my $nbentries=0;
    while ((my ($head, $seq)=__nextEntry())[0]) {
      $nbentries++;
      $head=~s/>/>$acPrefix/ or die "entry header does not start with '>': $head";
      $seq=join ('', shuffle  split(//, $seq));
      print  "$head\n";
      print __prettySeq($seq)."\n";
    }
    print STDERR "shuffled $nbentries\n" if $verbose;
  } else {
    $acPrefix ||='SHFLPLUS_';
    my $enz=qr/($shuffle_enzyme)/;
    my %pepts;
    my @vcrc = $shuffle_reshuffleCleavPept_CRCLen && Bit::Vector->new($maxBitCRC, $nbVCRC||1);
    my $nbentries=0;
    while ((my ($head, $seq)=__nextEntry())[0]) {
      $nbentries++;
      while ($seq=~/$enz/g) {
	$_=$1;
	my $l=length($_);
	if ($l>=$shuffle_reshuffleCleavPept_minLength) {
	  if ($shuffle_reshuffleCleavPept_CRCLen) {
	    my ($i, $c)=crc($_, 64);
	    $i%=$nbVCRC;
	    $c%=$maxBitCRC;
	    $vcrc[$i]->Bit_On($c);
	  } else {
	    $pepts{$_}=undef;
	  }
	}
      }
    }
    __setInput();
    my @histoReshuffled;
    my $nreshuffled=0;
    my $nFinalReshuffle=0;
    my $donotReadNext;
    my ($head, $seq);
    my $seqbak;
    my $nreshuffperseq;
    #donotReadNext is useed to resshuffle the whole sequence without reading the next one
    my $iseq=0;
    while ($donotReadNext || (($head, $seq)=__nextEntry())[0]) {
      #print  "(".__LINE__.") [$head]\n[$seq]\n";
      if ($donotReadNext) {
	$seq=$seqbak;
	undef $donotReadNext;
      } else {
	$nreshuffperseq=0 ;
	$seqbak=$seq;
	$iseq++;
      }
      #die"TEMP END"  if $iseq>100;
      if ($nreshuffperseq>$imaxreshuffleFinal) {
	$head=~/>(\S+)/;
	if ($pg) {
	  $pg->message("reshuffling the whole sequence without control [$1]");
	} else {
	  warn "reshuffling the whole sequence without control [$1]\n";
	}
	$head=~s/>/>$acPrefix/ or die "entry header does not start with '>': $head";
	my $newseq=join ('', shuffle (split //, $seq));
	print "$head\n".__prettySeq($newseq)."\n";

	$histoReshuffled[$nreshuffperseq]++;
	$nFinalReshuffle++;
	undef $donotReadNext;
	next;
      }
      $seq=join ('', shuffle (split //, $seq));
      my $newseq="";
      while ($seq && $seq=~/$enz/) {
	if ($nreshuffperseq && (($nreshuffperseq % $imaxreshuffle) == 0)) {
	  $head=~/>(\S+)/;
	  if ($pg) {
	    $pg->message("reshuffling the whole sequence [$1] ($nreshuffperseq/$imaxreshuffleFinal)");
	  } else {
	    warn "reshuffling the whole sequence [$1] ($nreshuffperseq/$imaxreshuffleFinal)\n";
	  }
	  $nreshuffperseq++;
	  $donotReadNext=1;
	  last;
	}
	my $pept=$1;

	my $reshuffle;
	if($newseq){
	  my $charniere=substr($newseq, length($newseq) -1, 1).substr($pept, 0, 1);
	  if(($charniere=~s/$enz//) && !$charniere){
	    $reshuffle=1;
	  }
	}


	unless ($reshuffle){
	  if (length($pept)<$shuffle_reshuffleCleavPept_minLength) {
	    $newseq.=$pept;
	    $seq=substr($seq, length($pept));
	    next;
	  }
	  if ($shuffle_reshuffleCleavPept_CRCLen) {
	    my ($i, $c)=crc($pept, 64);
	    $i%=$nbVCRC;
	    $c%=$maxBitCRC;
	    $reshuffle =  $vcrc[$i]->bit_test($c);
	  } else {
	    $reshuffle = exists $pepts{$_};
	  }
	}
	if ($reshuffle) {
	  $nreshuffled++;
	  $nreshuffperseq++;
	  if ($seq=~/(.{50})(.+)/) {
	    my ($s1, $s2)=($1, $2);
	    $seq=join ('', shuffle (split //, $s1)).$s2;
	  } else {
	    $seq=join ('', shuffle (split //, $seq));
	  }
	  next;
	} else {
	  $newseq.=$pept;
	  $seq=substr($seq, length($pept));
	  $histoReshuffled[$nreshuffperseq]++;
	  $nreshuffperseq=0;
	}
      }
      unless($donotReadNext){
	$head=~s/>/>$acPrefix/ or die "entry header does not start with '>': $head";
	print "$head\n".__prettySeq($newseq)."\n";
	$histoReshuffled[$nreshuffperseq]++;
	#print "(".__LINE__.") [$head]\n[$newseq]\n";
      }
    }
    if ($verbose) {
      print STDERR "reshuffled pept/nb sequences: $nreshuffled/$nbentries\n";
      print STDERR "nb final (no-check)  seq reshuffling: $nFinalReshuffle\n";
      print STDERR "#seq\t#nb reshuffled peptides\n";
      foreach (0..$#histoReshuffled) {
	next unless $histoReshuffled[$_];
	print STDERR "$_\t$histoReshuffled[$_]\n";
      }
    }
  }
} else {
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
  eval{
    if ((!$noProgressBar) && ($inFile !~ /gz$/i) && ($inFile ne '-')&& -t STDIN && -t STDOUT){
      require Term::ProgressBar;
      $size=(stat $inFile)[7];
      $pg=Term::ProgressBar->new ({name=> "parsing ".basename($inFile),
				   count=>$size,
				   ETA=>'linear',
				   remove=>1
				  });
    }
    $nextpgupdate=0;
    $readsize=0;
  };
  if ($@){
    warn "could not use Term::ProgressBar (consider installing the module for more interactive use";
  }
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

