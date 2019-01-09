#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use File::Basename;

my($long,$outdir,$contig,$identify,$ppn);
my $LAST="/home/fanyucai/software/last/last-921";
my $BWA="/home/fanyucai/software/bwa/bwa-0.7.12/";
my $GapCloser="/home/fanyucai/software/GapCloser";
my $SNAP="/home/fanyucai/software/SNAP/snap-0.15.4-linux/";
my $samtools="/home/fanyucai/software/samtools/samtools-v1.4/bin";
my $redundans="/home/fanyucai/software/redundans/redundans/redundans.py";
my $env="export PATH=$LAST/scripts:$LAST/src:$BWA:$GapCloser:$SNAP:$samtools:$redundans:\$PATH && export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64:\$LD_LIBRARY_PATH ";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
#FastaIndex\pyScaf install use pip
$identify||=0.51;
$outdir||=getcwd;
my $overlap||=0.66;
my $nogap||=0;
GetOptions(          
    "l:s"=>\$long,      
    "c:s"=>\$contig,
    "o:s"=>\$outdir,
    "i:s"=>\$identify,
    "overlap:s"=>\$overlap,
    "nogap:s"=>\$nogap,
           );
sub usage{
    print qq{
Redundans takes as input assembled contigs and returns scaffolded homozygous genome assembly.     
usage:
Long reads:
perl $0 -l subreads.fq.gz -c contigs.fna -o $outdir

options:
-l                  long reads
-o                  output directory
-c                  FASTA file with contigs / scaffolds
-i                  min. identity [0.51]
-overlap            min. overlap  [0.66]
-nogap              nogapclosing[1],gap-closing[0]

Email:fanyucai1\@126.com
2018.1.31
    };
    exit;
}
if(!$long ||!$contig)
{
    &usage();
}
$redundans.=" -i $long -f $contig -o $outdir/redundans -t 20 ";
$redundans .=" --identity $identify --overlap $overlap ";
if($nogap==1)
{
   $redundans.=" --nogapclosing ";
}

system "echo '$env && python $redundans '>$outdir/redundans.sh";
`perl $qsub --ppn 15 $outdir/redundans.sh`;
