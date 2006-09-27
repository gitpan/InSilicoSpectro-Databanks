#!/usr/bin/env perl
use strict;
use Carp;
use Pod::Usage;

=head1 NAME

uniprotdata2fasta.pl

=head1 DESCRIPTION

Read a uniprot .dat (sprot, trembl...) and convert it into a fasta file, expanding VAR_SPLIC, CHAIN, PEPT annotations.

FT annotation, such as MOD_RES, VARIANT... are kept in the fasta header coherently with the sequence modifictaion

=head1 SYNOPSIS

uniprotdata2fasta.pl --in=uniprot_sprot.dat --out=uniprot_sprot.fasta

=head1 ARGUMENTS


=head3 --in=file

A .dat file [default is stdin]

=head3 -out=file

A .fasta file [default is stdout]


=head1 OPTIONS

=head3 --noderived

does not produce the derived sequences (VAR_SEQ, CHAIN...)

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



use InSilicoSpectro::Databanks::DBEntryUniprot;
use SelectSaver;
use File::Basename;

use Getopt::Long;
my($datFile, $fastaFile, $noDerivedForm, $help, $man, $verbose);
$datFile='-';
$fastaFile='-';

if (!GetOptions(
		"in=s"=>\$datFile,
		"out=s" => \$fastaFile,

		"noderived"=>\$noDerivedForm,

                "help" => \$help,
                "man" => \$man,
                "verbose" => \$verbose,
               )
    || $help || $man){


  pod2usage(-verbose=>2, -exitval=>2) if(defined $man);
  pod2usage(-verbose=>1, -exitval=>2);
}


my ($pg, $size, $nextpgupdate, $readsize);
$readsize=0;
$nextpgupdate=0;

eval{
  if (($datFile ne '-')&& -t STDIN && -t STDOUT){
    require Term::ProgressBar;
    $size=(stat $datFile)[7];
    $pg=Term::ProgressBar->new ({name=> "parsing ".basename($datFile), count=>$size});
  }
};
if ($@){
  warn "could not use Term::ProgressBar (consider installing the module for more interactive use";
}


open (FDIN, "<$datFile") or die "could not open for reading [$datFile]:$!";

my $saver;
if ($fastaFile ne '-'){
  open (FDOUT, ">$fastaFile") or die "could not open for writing [$fastaFile]:$!";
  $saver=new SelectSaver(\*FDOUT);
}


$/="//\n";
while (<FDIN>){
  $readsize+=length $_;
  $nextpgupdate=$pg->update($readsize) if $pg && $readsize>$nextpgupdate;
  my $dbu=InSilicoSpectro::Databanks::DBEntryUniprot->new;
  $dbu->readDat($_);
  unless($noDerivedForm){
    my @tmp=$dbu->generateDerivedForms();
    foreach (@tmp){
      $_->printFasta;
    }
    unless (@tmp){
     $dbu->printFasta;
   }
  }else{
    $dbu->printFasta;
  }
}
$pg->update($size) if $pg;
exit(0);

