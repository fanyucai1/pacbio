#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;

my $PBcR="/home/fanyucai/software/Celera_Assembler/wgs-8.3rc2/Linux-amd64/bin/";
my $sprai="/home/fanyucai/software/sprai/sprai-v0.9.9.23/bin/";
my $blast="/home/fanyucai/software/blast+/ncbi-blast-2.6.0+/bin";
my $java="/home/fanyucai/software/java/jre1.8.0_121/";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
my $shuffle="/home/fanyucai/software/velvet/velvet_1.2.10/contrib/shuffleSequences_fasta/shuffleSequences_fastq.pl";
my $bowtie="/home/fanyucai/software/bowtie/bowtie-1.2/";
my $bowtie2="/home/fanyucai/software/bowtie2/bowtie2-2.3.1/";
my $env="export JAVA_HOME=$java && export PATH=$bowtie:$bowtie2:\$JAVA_HOME/bin:\$PATH && export CLASSPATH=.:\$JAVA_HOME/lib/dt.jar:\$JAVA_HOME/lib/tools.jar && export LD_LIBRARY_PATH=/home/fanyucai/software/tbb/tbb-2017_U5/build/linux_intel64_gcc_cc4.4.7_libc2.12_kernel2.6.32_release/:/home/fanyucai/software/glibc/glibc-v2.14/lib/:\$LD_LIBRARY_PATH";
my $smrtanalysis="/home/Softwares/smrtanalysis_2.3.0/current/etc/setup.sh";

my ($subreads,$genomesize,$minlen,$prefix,$config,$outdir,$queue,$pe1,$pe2);
$minlen||=1000;
$outdir||=getcwd;
$pe1||="0";
$queue||="fat";
GetOptions(
    "r:s"=>\$subreads,
    "g:s"=>\$genomesize,
    "m:s"=>\$minlen,
    "p:s"=>\$prefix,
    "o:s"=>\$outdir,
    "a:s"=>\$pe1,
    "b:s"=>\$pe2,
    "q:s"=>\$queue,
           );

sub usage{
    print qq{
This script will run PBcR assembly.
usage:
perl $0 -r subreads.fastq -g 4500000 -p Ecoli -o /path/to/directory -q fat
    or
perl $0 -r subreads.fastq -g 4500000 -p Ecoli -o /path/to/directory -a pe1.fastq -b pe2.fastq -q fat
options:
-r          subreads (fastq from Pacbio)
-g          genome size
-m          minlength using analysis(default:500)
-p          prefix of output
-o          outputdirectory(default:$outdir)
-a          NGS data(fastq 5 reads)
-b          NGS data(fastq 3 reads)
-q          which queue you want(default:fat)

Email:fanyucai1\@126.com
2017.6.15
    };
    exit;
}
if(!$subreads||!$genomesize||!$prefix)
{
    &usage();
}
`echo "merSize=14">$outdir/$prefix.spec`;
`echo "gridEngine=PBS">>$outdir/$prefix.spec`;
`echo "useGrid = 1">>$outdir/$prefix.spec`;

open(PBCR,">$outdir/PBcR.sh");
if($pe1 eq "0")
{
    print PBCR "source $smrtanalysis && $env && cd $outdir && $PBcR/PBcR -length $minlen -partitions 300 -l $prefix -s $prefix.spec -fastq $subreads genomeSize=$genomesize";
}
else
{
    print PBCR "perl $shuffle $pe1 $pe2 Illumina.fq\n";
    print PBCR "source $smrtanalysis && $env && cd $outdir && $PBcR/fastqToCA -libraryname illumina -technology illumina -type sanger -innie -reads Illumina.fq > illumina.frg\n";
    print PBCR "source $smrtanalysis && $env && cd $outdir && $PBcR/PBcR -length $minlen -partitions 300 -l $prefix -s $prefix.spec -fastq $subreads genomeSize=$genomesize illumina.frg"; 
}
my $line=`wc -l $outdir/PBcR.sh`;
chomp($line);

`perl $qsub --lines $line --queue fat $outdir/PBcR.sh`;
