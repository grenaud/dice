#!/usr/bin/perl


use strict;
use warnings;
use Data::Dumper;
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
my $logs2plot  = $pathdir."/logs2plot.R";
my $logs2text = $pathdir."/logs2text.R";

my $help;

my $anchorPop     = "YRI";
my $outputprefix  = "diceout";
my $mappability   = $pathdir."/../mapability/all.1kregions.gz";
my $alleleFreqNuc = $pathdir."/../alleleFreqNuc/";
my $contpop       = "ASW,BEB,CDX,CEU,CHB,CHS,CLM,FIN,GBR,IBS,JPT,MXL,PJL,PUR,TSI,YRI";
my $tauRange      = "1e-06,1";

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\n
 This script produces a Makefile that contains the commands to convert
 a BAM file in native dice format and runs it on the different files.


\n\n usage:\t".$0." <options>  [reference FASTA] [BAM file 1] [BAM file 2] ... > dice.make

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
  "\t--tau\t[tauLow,tauHigh]\tRange for tauA and tauC (Default: ".$tauRange.")\n".
    "\n\n".
    "";
  exit;
}



if ( @ARGV < 1 or ! GetOptions('help|?' => \$help, 'tau=s' => \$tauRange, 'anch=s' => \$anchorPop,'out=s' => \$outputprefix,'reg=s' => \$mappability,'alfr=s' => \$alleleFreqNuc,'cont=s' => \$contpop) or defined $help ){
  usage();
}

my @arrayBAMfiles;
my $fastaf;
my $ar=$#ARGV;
while($ar!=-1){

  if( ($ARGV[ $ar   ] =~ /fa$/ ) ||
      ($ARGV[ $ar   ] =~ /fasta$/ ) ){
    $fastaf = $ARGV[ $ar   ];
    last;
  }else{
    push(@arrayBAMfiles, $ARGV[ $ar ]);
    my $bamf = $ARGV[ $ar ];
    fileExists($bamf);
    fileExists($bamf.".bai");
    warn "Found bam:\t".$bamf."\n";
  }
  $ar-=1;
}

#my $bamf   = $ARGV[ $#ARGV   ];


exfileExists($bam2dice);
exfileExists($dice);
exfileExists($logs2text);
exfileExists($log2plot);
exfileExists($logs2plot);


warn "Found program:\t".$bam2dice."\n";
warn "Found program:\t".$dice."\n";
warn "Found program:\t".$log2plot."\n";
warn "Found program:\t".$logs2plot."\n";
warn "Found program:\t".$logs2text."\n";


fileExists($fastaf);

warn "Found fasta:\t".$fastaf."\n";

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
my $stringToPrint="";

my %hashbamkey;
foreach my $bamf (@arrayBAMfiles) {
  my $bamkey = substr($bamf,0,-4);
  my @bamkeya = split("/",$bamkey);
  $bamkey = $bamkeya[ $#bamkeya ];
  warn "key ".$bamkey."\n";
  if(exists $hashbamkey{ $bamkey }){
    die "Problem, 2 BAM files have the same prefix/suffixes ".$bamkey."\n";
  }else{
    $hashbamkey{ $bamkey } = 1;
  }

  my $b2dc=$bam2dice." --anch  ".$anchorPop." -o ".$outputprefix."_".$bamkey." -2p ".$fastaf."  ".$bamf." ".$mappability." ";

  foreach my $fifr (split(",",$contpop)) {
    $b2dc = $b2dc ." ". $alleleFreqNuc."/".$fifr.".mst.gz";
  }

  my $targetd2b;
  if ($contp1 eq $anchorPop) {
    $targetd2b = $outputprefix."_".$bamkey."_Cont_Anch_".$contp1.".dice";
  } else {
    $targetd2b = $outputprefix."_".$bamkey."_Cont_".$anchorPop."_Anch_".$contp1.".dice";
  }
  push(@arrayTargets,$targetd2b);

  foreach my $fifr (split(",",$contpop)) {
    my $targetd;
    my $targetp;
    my $target;
    if ($fifr eq $anchorPop) {
      $targetd = $outputprefix."_".$bamkey."_Cont_Anch_".$fifr.".dice.out.gz";
      $targetp = $outputprefix."_".$bamkey."_Cont_Anch_".$fifr.".dice.pdf";
      $target  = $outputprefix."_".$bamkey."_Cont_Anch_".$fifr.".dice";
    } else {
      $targetd = $outputprefix."_".$bamkey."_Cont_".$fifr."_Anch_".$anchorPop.".dice.out.gz";
      $targetp = $outputprefix."_".$bamkey."_Cont_".$fifr."_Anch_".$anchorPop.".dice.pdf";
      $target  = $outputprefix."_".$bamkey."_Cont_".$fifr."_Anch_".$anchorPop.".dice";
    }

    push(@arrayTargets,$targetd);
    push(@arrayTargets,$targetp);

    push(@arrayFiles,$targetd);
    push(@arrayFiles,$targetp);
    push(@arrayFiles,$target);

  }

  push(@arrayTargets,$outputprefix."_".$bamkey.".dice.txt");
  push(@arrayFiles,  $outputprefix."_".$bamkey.".dice.txt");
  push(@arrayTargets,$outputprefix."_".$bamkey."_c.pdf");
  push(@arrayFiles,  $outputprefix."_".$bamkey."_c.pdf");






  my @arraydiceout;
  my @arraydiceToZip;
  $stringToPrint.= $targetd2b.":\n\t".$b2dc."\n\n";
  foreach my $fifr (split(",",$contpop)){
    my $targetd;
    my $targetp;
    my $target;
    if($fifr eq $anchorPop){
      $targetd = $outputprefix."_".$bamkey."_Cont_Anch_".$fifr.".dice.out.gz";
      $targetp = $outputprefix."_".$bamkey."_Cont_Anch_".$fifr.".dice.pdf";
      $target  = $outputprefix."_".$bamkey."_Cont_Anch_".$fifr.".dice";
    }else{
      $targetd = $outputprefix."_".$bamkey."_Cont_".$fifr."_Anch_".$anchorPop.".dice.out.gz";
      $targetp = $outputprefix."_".$bamkey."_Cont_".$fifr."_Anch_".$anchorPop.".dice.pdf";
      $target  = $outputprefix."_".$bamkey."_Cont_".$fifr."_Anch_".$anchorPop.".dice";
    }

    $stringToPrint.= $targetd.": ".$targetd2b."\n\t".$dice."  -tA ".$tauRange." -tC ".$tauRange." -2p ".$target." |gzip > ".$targetd."\n\n";
    $stringToPrint.= $targetp.": ".$targetd.  "\n\t".$log2plot." ".$targetd."        ".$targetp."\n\n";

    push(@arraydiceout,  $targetd);
    push(@arraydiceToZip,$target );

    #."\n\n";
    #  ~/projects/dice/src/dice -2p -o diceout_Cont_Anch_CEU.dice.out diceout_Cont_Anch_CEU.dice  &
    #~/projects/dice/src/dice -2p -o diceout_Cont_CHB_Anch_CEU.dice.out diceout_Cont_CHB_Anch_CEU.dice &
  }

  $stringToPrint.=  $outputprefix."_".$bamkey.".dice.txt: ".join(" ",@arraydiceout)."\n\t".$logs2text." ".join(" ",@arraydiceout). " > ".$outputprefix."_".$bamkey.".dice.txt\n";

$stringToPrint.=  $outputprefix."_".$bamkey."_c.pdf: ".join(" ",@arraydiceout)."\n\t".$logs2plot." ".$outputprefix."_".$bamkey." ".join(" ",@arraydiceout). "\n";

  foreach my $fidice (@arraydiceToZip){
    $stringToPrint.= "\n".$fidice.".gz: ".$outputprefix."_".$bamkey.".dice.txt\n\t"."gzip ".$fidice."\n";
    push(@arrayTargets,$fidice.".gz");
    push(@arrayFiles,  $fidice.".gz");
  }

  #warn Dumper(@arrayTargets);
  #warn Dumper(@arrayFiles);

}

#todo for each bam a log2plot and log2text for each file independently
#zippe .dice files
print "SHELL := /bin/bash\n\nall:\t".join(" ",@arrayTargets)."\n\nclean:\n\trm -vf ".join(" ",@arrayFiles)."\n\n".$stringToPrint."\n\n";


