#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;

my $minimap="/home/fanyucai/software/minimap/minimap-master/";
my $miniasm="/home/fanyucai/software/miniasm/miniasm-master/";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
my $racon="/home/fanyucai/software/racon/racon/bin";
my $env="export LD_LIBRARY_PATH=/home/Softwares/gcc-5.2.0/lib64/:\$LD_LIBRARY_PATH";

my($fq,$outdir,$queue,$thread);
$queue||="fat";
$outdir||=getcwd;
$thread||=10;
sub usage{
    print qq{
This script will assembly genome using miniasm,minimap and racon only PacBio reads.
URL:
https://github.com/lh3/miniasm
https://github.com/isovic/racon

usage:
perl $0 -i reads.fq -o /path/to/directory/ -q fat
options:
-i              raw reads(fastq)
-o              outputdirectory(default:$outdir)
-t              thread(default:10)
-q          specify the queue to use(default:fat)

Email:fanyucai1\@126.com
2017.6.19
    };
    exit;
}
GetOptions(
    "i:s"=>\$fq,       
    "o:s"=>\$outdir,       
    "q:s"=>\$queue,
    "t:s"=>\$thread,
           );
if(!$fq||!$queue)
{
    &usage();
}

system "mkdir -p $outdir";
#################################################################1
open(AS1,">$outdir/minimap.sh");#mapping
print AS1 "$minimap/minimap -Sw5 -L100 -m0 -t$thread $fq $fq | gzip -1 > $outdir/reads.paf.gz\n";
if(! -e "$outdir/reads.paf.gz")
{
    `perl $qsub --queue $queue $outdir/minimap.sh`;
}
print "step1:minimap(mapping) has done\n";
#################################################################2
open(AS2,">$outdir/miniasm.sh");#assembly
print AS2 "$miniasm/miniasm -f $fq $outdir/reads.paf.gz > $outdir/contigs.gfa \n";
if(! -e "$outdir/contigs.gfa")
{
   `perl $qsub --queue $queue $outdir/miniasm.sh`;
}
print "step2 miniasm(assembly) has done\n";
#################################################################3
open(AS3,">$outdir/racon.sh");
`awk \'/\^S/ {print \">\"\$2\"\\n\"\$3}\' $outdir/contigs.gfa | fold > $outdir/contigs.gfa.fasta`;#convert gfa to fasta
print AS3 "$minimap/minimap -t$thread $outdir/contigs.gfa.fasta $fq > $outdir/contigs.paf\n";#mapping contig to subreads
print AS3 "$env && $racon/racon $fq $outdir/contigs.paf $outdir/contigs.gfa.fasta $outdir/out_consensus.fasta\n";#Consensus
my $line=`wc -l $outdir/racon.sh`;
chomp($line);
if(! -e "$outdir/out_consensus.fasta")
{
    `perl $qsub --queue $queue --lines $line $outdir/racon.sh`;
}
print "step3 Racon(Consensus) has done\n";
