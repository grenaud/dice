#!/usr/bin/perl


use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;

my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);



sub exfileExists{
  my ($exeFile) = @_;

  if (!( -x $exeFile)) {
    die "Executable file ".$exeFile." does not exist\n";
  }
}

sub fileExists{
  my ($nFile) = @_;

  if (!( -e $nFile)) {
    die "File ".$nFile." does not exist\n";
  }
}

sub dirExists{
  my ($exeFile) = @_;

  if (!( -d $exeFile)) {
    die "Directory ".$exeFile." does not exist\n";
  }
}


my $bam2dice  = $pathdir."/BAM2DICE";
my $dice      = $pathdir."/dice";
my $log2plot  = $pathdir."/log2plots.R";
my $logs2text = $pathdir."/logs2text.R";

my $help;

my $anchorPop     = "YRI";
my $outputprefix  = "diceout";
my $mappability   = "/home/gabriel/projects/dice/mapability/all.1kregions.gz";
my $alleleFreqNuc = "/home/gabriel/projects/dice/alleleFreqNuc/";
my $contpop       = "ASW,BEB,CDX,CEU,CHB,CHS,CLM,FIN,GBR,IBS,JPT,MXL,PJL,PUR,TSI,YRI";

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\n
 This script produces a Makefile that contains the commands to convert
 a BAM file in native dice format and runs it on the different files.


\n\n usage:\t".$0." <options>  [reference FASTA] [BAM file] > dice.make

 The BAM file should be sorted and indexed. The reference fasta
 should be have an index from faidx. This reference fasta has to be
 the same used for alignment.
 You then run make -f dice.make
 tip: use make -j [cores] to use multiple cores.


 Options:\n".
  "\t--anch\t[anchor population]\tThis population will be used as anchor (Default: ".$anchorPop.")\n".
  "\t--out\t[output prefix]\t\tThis is the output prefix (Default: ".$outputprefix.")\n".
  "\t--reg\t[file]\t\t\tSet of regions to consider (Default: ".$mappability.")\n".
  "\t\t\t\t\tthe format is chr:start-end\n".
  "\t\t\t\t\ta set of regions with high mappability is recommended\n".
  "\t--alfr\t[directory]\t\tDirectory containing the allele frequencies (Default: ".$alleleFreqNuc.")\n".
  "\t--cont\t[pop1,pop2]\t\tShorthand for the populations to use as contaminant, comma separated\n".
    "\t\t\t\t\t(Default: ".$contpop.")\n".
    "";
  exit;
}



if ( @ARGV < 1 or ! GetOptions('help|?' => \$help, 'anch=s' => \$anchorPop,'out=s' => \$outputprefix,'reg=s' => \$mappability,'alfr=s' => \$alleleFreqNuc,'cont=s' => \$contpop) or defined $help ){
  usage();
}


my $bamf   = $ARGV[ $#ARGV   ];
my $fastaf = $ARGV[ $#ARGV-1 ];

exfileExists($bam2dice);
exfileExists($dice);
exfileExists($logs2text);
exfileExists($log2plot);


warn "Found program:\t".$bam2dice."\n";
warn "Found program:\t".$dice."\n";
warn "Found program:\t".$log2plot."\n";
warn "Found program:\t".$logs2text."\n";

fileExists($bamf);
fileExists($fastaf);
warn "Found bam:\t".$bamf."\n";
warn "Found fasta:\t".$fastaf."\n";
fileExists($bamf.".bai");
fileExists($fastaf.".fai");
fileExists($mappability);
warn "Found regions:\t".$mappability."\n";

dirExists($alleleFreqNuc);
warn "Found dir.:\t".$alleleFreqNuc."\n";
fileExists($alleleFreqNuc."/".$anchorPop.".mst.gz");
warn "Found anch:\t".$alleleFreqNuc."/".$anchorPop.".mst.gz\n";

my $contp1="";
foreach my $fifr (split(",",$contpop)){
  fileExists($alleleFreqNuc."/".$fifr.".mst.gz");
  warn "Found cont:\t".$alleleFreqNuc."/".$fifr.".mst.gz\n";
  $contp1=$fifr;
}


my @arrayTargets;
my @arrayFiles;

my $b2dc=$bam2dice." --anch  ".$anchorPop." -o ".$outputprefix." -2p ".$fastaf."  ".$bamf." ".$mappability." ";

foreach my $fifr (split(",",$contpop)){
  $b2dc = $b2dc ." ". $alleleFreqNuc."/".$fifr.".mst.gz";
}

my $targetd2b;
if($contp1 eq $anchorPop){
  $targetd2b = $outputprefix."_Cont_Anch_".$contp1.".dice";
}else{
  $targetd2b = $outputprefix."_Cont_".$anchorPop."_Anch_".$contp1.".dice";
}
push(@arrayTargets,$targetd2b);

foreach my $fifr (split(",",$contpop)){
  my $targetd;
  my $targetp;
  my $target;
  if($fifr eq $anchorPop){
    $targetd = $outputprefix."_Cont_Anch_".$fifr.".dice.out.gz";
    $targetp = $outputprefix."_Cont_Anch_".$fifr.".dice.pdf";
    $target  = $outputprefix."_Cont_Anch_".$fifr.".dice";
  }else{
    $targetd = $outputprefix."_Cont_".$anchorPop."_Anch_".$fifr.".dice.out.gz";
    $targetp = $outputprefix."_Cont_".$anchorPop."_Anch_".$fifr.".dice.pdf";
    $target  = $outputprefix."_Cont_".$anchorPop."_Anch_".$fifr.".dice";

  }
  push(@arrayTargets,$targetd);
  push(@arrayTargets,$targetp);

  push(@arrayFiles,$targetd);
  push(@arrayFiles,$targetp);
  push(@arrayFiles,$target);

}

push(@arrayTargets,$outputprefix.".dice.txt");
push(@arrayFiles,  $outputprefix.".dice.txt");

print "SHELL := /bin/bash\n\nall:\t".join(" ",@arrayTargets)."\n\nclean:\n\trm -vf ".join(" ",@arrayFiles)."\n\n";

print $targetd2b.":\n\t".$b2dc."\n\n";

my @arraydiceout;
foreach my $fifr (split(",",$contpop)){
  my $targetd;
  my $targetp;
  my $target;
  if($fifr eq $anchorPop){
    $targetd = $outputprefix."_Cont_Anch_".$fifr.".dice.out.gz";
    $targetp = $outputprefix."_Cont_Anch_".$fifr.".dice.pdf";
    $target  = $outputprefix."_Cont_Anch_".$fifr.".dice";
  }else{
    $targetd = $outputprefix."_Cont_".$anchorPop."_Anch_".$fifr.".dice.out.gz";
    $targetp = $outputprefix."_Cont_".$anchorPop."_Anch_".$fifr.".dice.pdf";
    $target  = $outputprefix."_Cont_".$anchorPop."_Anch_".$fifr.".dice";

  }

  print $targetd.": ".$targetd2b."\n\t".$dice." -2p ".$target." |gzip > ".$targetd."\n\n";
  print $targetp.": ".$targetd.  "\n\t".$log2plot." ".$targetd."        ".$targetp."\n\n";

  push(@arraydiceout,$targetd);
  #."\n\n";
  #  ~/projects/dice/src/dice -2p -o diceout_Cont_Anch_CEU.dice.out diceout_Cont_Anch_CEU.dice  &
  #~/projects/dice/src/dice -2p -o diceout_Cont_CHB_Anch_CEU.dice.out diceout_Cont_CHB_Anch_CEU.dice &
}


print $outputprefix.".dice.txt: ".join(" ",@arraydiceout)."\n\t".$logs2text." ".join(" ",@arraydiceout). " > ".$outputprefix.".dice.txt";

