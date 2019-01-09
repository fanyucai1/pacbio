#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;

my($outdir,$ref,$input,$node,$ppn,$queue);
my $quiver="/smrtlinks/smrtlink/smrtcmds/bin/quiver";
my $pbalign="/smrtlinks/smrtlink/smrtcmds/bin/pbalign";
my $samtools="/smrtlinks/smrtlink/smrtcmds/bin/samtools";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
$outdir||=getcwd;
$ppn||=25;
$queue||="all";
GetOptions(
    "o:s"=>\$outdir,
    "r:s"=>\$ref,
    "i:s"=>\$input,
    "ppn:s"=>\$ppn,
    "queue:s"=>\$queue,
           );
sub usage{
    print qq{
Beacuse the enviroment PATH,you must run in node of fat02.
This script will GenomeicConsensus(https://github.com/PacificBiosciences/GenomicConsensus).
usage:
perl $0 -r contigs.fa -i subreads.bam -o /path/to/directory/ -ppn 25

options:
-r          the contig sequence
-i          subreads file(BAM)
-o          output directory
-ppn        the cpu you want(default:25)
-queue      which queue you will run

Email:fanyucai1\@126.com
2017.7.13
    };
    exit;
}
if(!$ref||!$input)
{
    &usage();
}
system "mkdir -p $outdir";
open(SH,">$outdir/quiver.sh");
#Using bamtools a list of subread files contained in a fofn can be converted to a single bam file.
#https://github.com/PacificBiosciences/PacBioFileFormats/wiki/BAM-recipes
####################################Beacuse the enviroment PATH,you must run in node of fat02.
print SH "$pbalign --algorithm blasr --nproc 30 $input $ref $outdir/aligned_subreads.bam\n";
print SH "$samtools faidx $ref\n";
print SH "$quiver $outdir/aligned_subreads.bam  -j 30 --algorithm=arrow -r $ref -o $outdir/variants.gff -o $outdir/consensus.fasta -o $outdir/consensus.fastq\n";
my $lines=`wc -l $outdir/quiver.sh`;
chomp($lines);
`perl $qsub --queue $queue --ppn $ppn --lines $lines $outdir/quiver.sh`;
