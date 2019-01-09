#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin);

my($ref,$read,$outdir,$line);
$outdir||=getcwd;
my $ngmlr="/home/fanyucai/software/ngmlr/ngmlr-0.2.3/ngmlr";
my $sniffles="/home/fanyucai/software/Sniffles/Sniffles-master/bin/sniffles-core-1.0.5/sniffles";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
my $samtools="/home/fanyucai/software/samtools/samtools-v1.4/bin/samtools";
GetOptions(
    "ref:s"=>\$ref,
    "read:s"=>\$read,
    "o:s"=>\$outdir, 
           );

if(!$read||!$ref)
{
    print "perl $0 -ref reference.fa -read pacbio.fq(fa) -o $outdir";
    exit;
}
open(SH,">$outdir/run.sh");
#ngmlr(https://github.com/philres/ngmlr)
print SH "$ngmlr -t 8 -r $ref -q $read -o $outdir/mapped.sam\n";
#sam2bam
print SH "$samtools view $outdir/reslut.sam -o $outdir/mapped.bam\n";
print SH "$samtools sort $outdir/mapped.bam -o $outdir/mapped_sort.bam\n";
#sniffles(https://github.com/fritzsedlazeck/Sniffles/wiki)
print SH "$sniffles -m $outdir/mapped.sort.bam -v $outdir/output.vcf\n";


$line=`wc -l $outdir/run.sh`;
chomp($line);
`perl $qsub --lines $line $outdir/run.sh`;