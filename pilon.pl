#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use Cwd;
use Getopt::Long;
use File::Basename;

my $samtools="/home/fanyucai/software/samtools/samtools-v1.4/bin/samtools";
my $bwa="/home/fanyucai/software/bwa/bwa-0.7.12/bwa";
my $pilon="/home/fanyucai/software/pilon/pilon-1.22-2.jar";
my $picard="/home/fanyucai/software/picard/picard.jar";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
my $index="/home/fanyucai/script/index.pl";
my ($ref,$a,$b,$outdir,$queue,$diploid,$lines);
$queue||="fat";
$outdir||=getcwd;
GetOptions(
    "r:s"=>\$ref,
    "o:s"=>\$outdir,
    "a:s"=>\$a,
    "b:s"=>\$b,
    "q:s"=>\$queue,
    "d:s"=>\$diploid,
           );

sub usage{
    print qq{
Details:
Paired-end reads from Illumina platform were aligned to the assembly using bwa mem, followed with duplication removal using Picard tools.Base-pair correction of the assembly was performed using Pilon. Pilon mostly corrected single insertions and deletions in regions enriched with homopolymer.Contigs or scaffolds shorter than 10 kb were excluded from the overall analysis to avoid results from spurious misassembly.

Usage:
perl $0 -a pe1.fq -b pe2.fq -o $outdir -r ref.fa -q fat -d t
Options:
-r              the fasta file(froce)
-a              5 read(froce)
-b              3 read(froce)
-o              outputdirectory(default:$outdir)
-q              queue run(default:fat)
-d              diploid or not(true or false,force)

Email:fanyucai1\@126.com
2017.6.2

[1]Seo J S, Rhie A, Kim J, et al. De novo assembly and phasing of a Korean human genome[J]. Nature, 2016, 538(7624): 243-7.
[2]Bickhart D M, Rosen B D, Koren S, et al. Single-molecule sequencing and chromatin conformation capture enable de novo reference assembly of the domestic goat genome[J]. Nature genetics, 2017, 49(4): 643.
[3]Tørresen O K, Star B, Jentoft S, et al. An improved genome assembly uncovers prolific tandem repeats in Atlantic cod[J]. BMC genomics, 2017, 18(1): 95.
[4]Walker B J, Abeel T, Shea T, et al. Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement[J]. PloS one, 2014, 9(11): e112963.
 };
    exit;
}
if(!$ref ||!$a||!$b||!$diploid)
{
    &usage();
}
my $base=basename $ref;
open(PILON,">$outdir/pilon.sh");
#index
system "perl $index -ref $ref";

#mapping
print PILON "cd $outdir && $bwa mem -M -t 8 $ref $a $b >$outdir/$base.sam\n";

#sam->sort_bam
print PILON "cd $outdir && $samtools sort -@ 8 -o $outdir/$base.bam $outdir/$base.sam && rm $outdir/$base.sam\n";

#Add read groups, coordinate sort and index using AddOrReplaceReadGroups
print PILON "java -Xmx32g -XX:ParallelGCThreads=10 -jar $picard AddOrReplaceReadGroups INPUT=$outdir/$base.bam OUTPUT=$outdir/$base.addRG.bam RGID=$base RGLB=$base RGPL=illumina RGPU=$base RGSM=$base SORT_ORDER=coordinate CREATE_INDEX=true TMP_DIR=$outdir/tmp && rm $outdir/$base.bam\n";

#marking duplication
print PILON "java -Xmx32g -XX:ParallelGCThreads=10 -jar $picard MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=$outdir/$base.addRG.bam OUTPUT=$outdir/$base.markdup.bam METRICS_FILE=$outdir/$base\_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=$outdir/tmp && rm $outdir/$base.addRG.bam\n";

#index bam file
print PILON "cd $outdir && $samtools index $outdir/$base.markdup.bam\n";

#pilon(https://github.com/skoren/PilonGrid/blob/master/pilonParallelSGE.sh)
#Supplemental Material for "Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation"
if($diploid=~/t/)
{
    print PILON "cd $outdir && java -jar -Xmx200g -XX:ParallelGCThreads=10 $pilon --threads 20 --vcf --fix bases,local --genome $ref --bam $base.markdup.bam --outdir $outdir --output $base.pilon --vcf --changes --diploid\n";
}
else
{
    print PILON "cd $outdir && java -jar -Xmx200g -XX:ParallelGCThreads=10 $pilon --threads 20 --vcf --fix bases,local --genome $ref --bam $base.markdup.bam --outdir $outdir --output $base.pilon --vcf --changes\n";
}
$lines=`wc -l $outdir/pilon.sh`;
chomp($lines);
`perl $qsub --queue $queue --lines $lines $outdir/pilon.sh`;
