#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;

my $nrFile;
my $outFile='-';
my ($nrFile, $taxodmpFile, $sortedData, $noProgressBar, $help, $verbose);


############### read command line options
if (!GetOptions(
		"in=s"=>\$nrFile,
		"taxo=s"=>\$taxodmpFile,
		"out=s"=>\$outFile,
		"sorted"=>\$sortedData,
		
		"noprogressbar"=>\$noProgressBar,
		
		"help"=>\$help,
		"verbose"=>\$verbose,
	       )
    || $help
   ){
  die <<EOT;
Usage: $0
Arguments
  --in=file              input file
  --taxo=file            gi_taxid_prot.dmp location

Options
  --out=file             output file

  --sorted               if gi id and taxnomy are sorted (musch faster and lower in memory)

  --help                 produces this message
  --verbose              to print verbose messages on STDERR

Examples:

EOT
}

die "no nr fasta (--in=file)" unless $nrFile;
die "notaxo.dmp (--taxo=file)" unless $taxodmpFile;

my ($nextpgupdate, $readsize, $pg);


unless($sortedData){
  setupPG($nrFile);

  my %ac2taxid;
  open (FD,"<$nrFile") or die "cannot open for reading [$nrFile]:$!";
  while (<FD>) {
    $readsize+=length $_;
    $nextpgupdate=$pg->update($readsize) if $pg && $readsize>=$nextpgupdate;
    next unless /^>gi\|(\w+)/;
    $ac2taxid{$1}=undef;
  }
  warn "nb of fasta entries=".scalar (keys %ac2taxid)."\n";
  close FD;

  setupPG($taxodmpFile);

  open (FD,"<$taxodmpFile") or die "cannot open for reading [$taxodmpFile]:$!";
  while (<FD>) {
    $readsize+=length $_;
    $nextpgupdate=$pg->update($readsize) if $pg && $readsize>=$nextpgupdate;
    chomp;
    my ($gi, $taxo)=split;
    next unless exists $ac2taxid{$gi};
    $ac2taxid{$gi}=$taxo;
  }
  close FD;

  my $nbnotaxid=0;
  foreach (keys %ac2taxid) {
    unless ($ac2taxid{$_}) {
      $nbnotaxid++;
      $ac2taxid{$_}=1;
    }
  }

  warn "no taxid found for $nbnotaxid/".scalar (keys %ac2taxid)." entries\n";
  warn "producing output file [$outFile]\n";
  setupPG($nrFile);
  open (FD, "<$nrFile") or die "cannot open for reading [$nrFile]:$!";
  open (OUT, ">$outFile") or die "cannot open for writing: [$outFile]: $!";
  while (<FD>) {
    $readsize+=length $_;
    $nextpgupdate=$pg->update($readsize) if $pg && $readsize>=$nextpgupdate;
    my $gi;
    chomp;
    if (s/^>((?:\w+)\|(\w+))\|(\S+)\|\s*(.*)/>$1 \\ID=$3 \\NCBITAXID=$ac2taxid{$2} \\DE=$4/) {
      $gi=$2;
      s/\cA/;/g;
    }
    unless(/^>/){
      s/[^A-Z]//g;
    }

    print OUT "$_\n";
    $ac2taxid{$gi}=undef if defined $gi;
  }
  close FD;
  close OUT;
} else {
  setupPG($nrFile);
  open (FD, "<$nrFile") or die "cannot open for reading [$nrFile]:$!";
  open (FDTAXO,"<$taxodmpFile") or die "cannot open for reading [$taxodmpFile]:$!";

  open (OUT, ">$outFile") or die "cannot open for writing: [$outFile]: $!";
  local $/="\n>";
  while (<FD>) {
    $readsize+=length $_;
    $nextpgupdate=$pg->update($readsize) if $pg && $readsize>=$nextpgupdate;

    s/^>//;
    my ($line, $seq)=split /\n/, $_, 2;

    my $gi;
    if ($line=~s/^((?:\w+)\|(\w+))\|(\S+)\|\s*(.*)/>$1 \\ID=$3 \\NCBITAXID=__TAXID__ \\DE=$4/) {
      $gi=$2;
      $line=~s/\cA/;/g;
      #locate gi
      my $taxid;
      local $/="\n";
      while(<FDTAXO>){
	my ($gitmp, $taxo)=split;
	if(($gitmp+0)>($gi+0)){
	  warn "rewind at line $line\n";
	  close FDTAXO;
	  open (FDTAXO,"<$taxodmpFile") or die "cannot open for reading [$taxodmpFile]:$!";
	  next;
	}
	if($gitmp eq $gi){
	  $taxid=$taxo;
	  last;
	}
      }
      die "no taxid for gi=[$gi] in $taxodmpFile" unless defined $taxid;
      $line=~s/NCBITAXID=__TAXID__/NCBITAXID=$taxid/;
    }
    $seq=~s/[^A-Z]//ig;

    print OUT "$line\n$seq\n";
  }
}
warn "completed\n";

sub setupPG{
  my $file=shift or die "no arg to setupPG";
  die "$file does not exist" unless -f $file;
  if ((!$noProgressBar)  && ($file ne '-')&& -t STDIN && -t STDOUT){
    require Term::ProgressBar;
    my $size=(stat $file)[7];
    $pg=Term::ProgressBar->new ({name=> "parsing ".basename($file),
				 count=>$size,
				 ETA=>'linear',
				 remove=>1
				});
    $nextpgupdate=0;
    $readsize=0;
  }
}
