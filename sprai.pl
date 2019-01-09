#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;

my ($subread,$outdir,$genomesize,$split,$minlen,$thread,$queue,$illumina,$use);
my $wgs="/home/fanyucai/software/Celera_Assembler/wgs-8.3rc2/Linux-amd64/bin";
my $blast="/home/fanyucai/software/blast+/ncbi-blast-2.6.0+/bin/";
my $sprai="/home/fanyucai/software/sprai/sprai-v0.9.9.23/bin/";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
my $glibc="export LD_LIBRARY_PATH=/home/fanyucai/software/glibc/glibc-v2.14/lib/:\$LD_LIBRARY_PATH";
$minlen||=500;
$outdir||=getcwd;
$queue||="big";
$use||=1;
GetOptions(
    "r:s"=>\$subread,
    "o:s"=>\$outdir,
    "z:s"=>\$genomesize,
    "min:s"=>\$minlen,
    "q:s"=>\$queue,
    "i:s"=>\$illumina,
           );
sub usage{
    print qq{
    This script you will run sprai to correct reads.
    usage:
    perl $0 -r all.subreads.fastq -z 5000000 -min 500 -o $outdir
    options:
    -r              subreads in fastq(fasta) format(force)
    -z              genome size(force)      
    -min            the subreads longer than or equal to this value will be corrected(defualt:500)
    -o              output directory(default:$outdir)
    -q              which queue you run(default:big01)
    
Email:fanyucai1\@126.com
2017.7.21
    };
    exit;
}
if(!$genomesize||!$subread||!$queue)
{
    &usage();
}
open(EC,">$outdir/ec.spec");
print EC "input_for_database ",$subread,"\n";
print EC "estimated_genome_size ",$genomesize,"\n";
print EC "estimated_depth 0","\n";
print EC "partition 20","\n";
print EC "evalue 1e-50\n";
print EC "trim 42\n";
print EC "ca_path $wgs\n";
print EC "word_size 18\n";
print EC "max_target_seqs 50\n";
print EC "min_len_for_query ",$minlen,"\n";
print EC "sprai_path ",$sprai,"\n";
print EC "blast_path ",$blast,"\n";
print EC "num_threads 20","\n";
print EC "use_one_subread 1\n";
print EC "direct_vote 0\n";
print EC "copy_blastdb 1\n";

#Correct errors
system "echo 'cd $outdir && perl $sprai/ezez_vx1.pl $outdir/ec.spec -ec_only'>sprai.sh";
my $lines=`wc -l $outdir/sprai.sh`;
chomp($lines);
`perl $qsub --queue $queue --lines $lines $outdir/sprai.sh`;

