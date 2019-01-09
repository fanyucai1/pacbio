#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;

my $OrganelleRef_PBA="/home/fanyucai/software/Organelle_PBA/Organelle_PBA/OrganelleRef_PBA";
my $sprai="export SPRAI_PATH=/home/fanyucai/software/sprai/sprai-v0.9.9.23/bin";
my $samtools="export SAMTOOLS_PATH=/home/fanyucai/software/samtools/samtools-v1.4/bin/";
my $blastn="export BLAST_PATH=/home/fanyucai/software/blast+/ncbi-blast-2.6.0+/bin/";
my $sspace_long="export SSPACELONG_PATH=/home/fanyucai/software/SSPACE/SSPACE-LongRead_v1-1/";#and blasr in this directory
my $bedtools="export BEDTOOLS_PATH=/home/fanyucai/software/bedtools/bedtools2/bin/";
my $CA="export CA_PATH=/home/fanyucai/software/Celera_Assembler/wgs-8.3rc2/Linux-amd64/bin/";
my $blasr="export BLASR_PATH=/home/Softwares/smrtanalysis_2.3.0/current/analysis/bin/";
my $glibc="export LD_LIBRARY_PATH=/home/fanyucai/software/glibc/glibc-v2.14/lib:\$LD_LIBRARY_PATH";
my $env="$glibc && $blasr && $samtools && $sprai && $blastn && $CA && $sspace_long && $bedtools";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";

my($input,$ref,$outdir,$type,$ppn);
$outdir||=getcwd;
$type||="fastq";
$ppn||=20;
GetOptions(
  "i:s"=>\$input,
  "r:s"=>\$ref,
  "o:s"=>\$outdir,
  "t:s"=>\$type,
  "ppn:s"=>\$ppn,
);
sub usage{
    print qq{
This script will perform a de-novo PacBio assemblies of any organelle (chloroplast or mitochondrial genomes) using several programs using Pacboio.
usage:
perl $0 -i pacbio.fasta -t fasta -r reference.fa -o $outdir
        or
perl $0 -i pacbio.fastq -t fastq -r reference.fa -o $outdir
options:
-i<input_pacbio>        input PacBio subreads
-t <input_type>         input type (fasta or fastq; default=fastq)
-r<fasta_reference>     organelle reference genome, fasta format
-o<output_dir>          output directory
-ppn                    cpu number(default:20)

Email:fanyucai1\@126.com
2017.9.22       
};
    exit;
}
if(!$type||!$input||!$ref)
{
    &usage();
}
system "mkdir -p $outdir/";
if($type=~/fastq/)
{
    system "echo 'cd $outdir && $env && $OrganelleRef_PBA -i $input -r $ref -o $outdir/ -b \'-nproc=40\' -s \'num_threads=40\' -l \'t=20\' '>$outdir/Organelle.sh";
}
else
{
    system "echo 'cd $outdir && $env && $OrganelleRef_PBA -t fasta -i $input -r $ref -o $outdir/ -b \'-nproc=40\' -s \'num_threads=40\' -l \'t=20\' '>$outdir/Organelle.sh";
}

system "cd $outdir && perl $qsub --ppn $ppn $outdir/Organelle.sh";

print "This process has done.";
