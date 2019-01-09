#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use Cwd;
use Bio::SeqIO;
use Bio::Seq;
use File::Basename;

my($genomesize,$outdir,$subreads);
$outdir||=getcwd;
my $bamtools="bamtools";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";

GetOptions(
    "g:s"=>\$genomesize,
    "o:s"=>\$outdir,
    "i:s"=>\$subreads,       
           );
sub usage{
    print qq{
This scipt will return a read length which could be used by falcon.
usage:
perl $0 -i subreads.bam -o $outdir -g 45000000

options:
-i          input subreads(bam,fa)
-g          genome size
-o          outdir default($outdir)

Email:fanyucai1\@126.com
2018.3.9
    };
    exit;
}
if(!$genomesize||!$subreads)
{
    &usage();
}
my $seqin;
if($subreads=~/bam$/i)
{
    system "echo '$bamtools convert -format fasta -in $subreads -out $outdir/subreads_tmp.fa'>$outdir/convert.sh";
    system "perl $qsub $outdir/convert.sh";
    $seqin= Bio::SeqIO->new( -format => 'Fasta', -file => "$outdir/subreads_tmp.fa");
}
else
{
   $seqin= Bio::SeqIO->new( -format => 'Fasta', -file => "$subreads");
}

my (@seqlen,%hash);
while((my $seqobj=$seqin->next_seq()))
{
    my $len=$seqobj->length();
    for(my $i=500;$i<=10000;$i+=500)
    {
        if($len>=$i)
        {
            $hash{$i}+=$len;
        }
    }
}
foreach my $key(keys %hash)
{
    $hash{$key}=int($hash{$key}/$genomesize);
}
open(OUT,">$outdir/Depth.xls");
print OUT "Length","\t","Depth(X)","\n";
for(my $i=500;$i<=10000;$i+=500)
{
    if($hash{$i}>=30)
    {
        print OUT $i,"\t",$hash{$i},"\n";
    }
}

