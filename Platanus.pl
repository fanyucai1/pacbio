#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use Cwd;

my (@mp,$outdir,$prefix,@pi,@mi,@pe,$queue,$ppn);
my $Platanus="/home/fanyucai/software/Platanus/platanus";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
$outdir||=getcwd;
$ppn||=20;
$queue||="fat";
GetOptions(
    "pe:s{1,}"=>\@pe,
    "mp:s{1,}"=>\@mp,
    "o:s"=>\$outdir,       
    "p:s"=>\$prefix,
    "pi:s{1,}"=>\@pi,
    "mi:s{1,}"=>\@mi,
    "ppn:s"=>\$ppn,
    "queue:s"=>\$queue,
           );
sub usage{
    print qq{
This script will use Platanus to assemble genome.
Reference:Efficient de novo assembly of highly heterozygous genomes from whole-genome shotgun short reads
usage:
perl $0 -pe lib1.1.fq lib1.2.fq lib2.1.fq lib2.2.fq -mp lib3.1.fq lib3.2.fq lib4.1.fq lib4.2.fq -p prefix -pi 300 500 -mi 1500 2000
options:
-pe         pe library files(fastq)
-pi         pe library insert length

-mp         mp library files(fastq)
-mi         mp library insert length

-p          prefix output
-o          output directory
-ppn        number of threads(default:20)
-queue      which queue you will run(default:fat)
   
Email:fanyucai1\@126.com
2018.1.19
    };
    exit;
}
if(!@pe ||!$outdir||!$prefix||!@pi)
{
    &usage();
}
system "mkdir -p $outdir";
open(AS,">$outdir/assembly.sh");
my $string;
for(my $i=0;$i<=$#pe;$i++)
{
    $string.=" $pe[$i] ";
}
if(@mp)
{
    for(my $i=0;$i<=$#mp;$i++)
    {
        $string.=" $mp[$i] ";
    }
}
print AS "$Platanus assemble -o $outdir/$prefix -k 32 -f $string -t 20 -s 5 -m 100\n";
`perl $qsub --queue $queue --ppn $ppn $outdir/assembly.sh`;

open(SC,">$outdir/scaffold.sh");
my $k=0;
$string="";
my $insert;
for(my $i=0;$i<=$#pe;$i=$i+2)
{
    $k++;
    $string.=" -IP$k $pe[$i] $pe[$i+1] ";
}
my $n=0;
for(my $i=0;$i<=$#pi;$i++)
{
    $n++;
    my $l=$pi[$i]-35;
    if($pi[$i]<=1000)
    {
        $insert.=" -d$n 35 -a$n $pi[$i] -n$n $l ";     
    }
}
if(@mp)
{
    for(my $i=0;$i<=$#mp;$i=$i+2)
    {
        $k++;
        $string.=" -OP$k $mp[$i] $mp[$i+1] ";
    }
    for(my $i=0;$i<=$#mi;$i++)
    {
        $n++;
        if($mi[$i]<=1500)
        {
            my $l=$mi[$i]-150;
            $insert.=" -d$n 150 -a$n $mi[$i] -n$n $l ";
        }
        if($mi[$i]<=2500 && $mi[$i]>=1500)
        {
            my $l=$mi[$i]-350;
            $insert.=" -d$n 350 -a$n $mi[$i] -n$n $l ";
        }
        if($mi[$i]>2500 && $mi[$i]<=3500)
        {
            my $l=$mi[$i]-400;
            $insert.=" -d$n 400 -a$n $mi[$i] -n$n $l ";
        }
        if($mi[$i]>3500 && $mi[$i]<=5000)
        {
            my $l=$mi[$i]-600;
            $insert.=" -d$n 600 -a$n $mi[$i] -n$n $l ";
        }
        if($mi[$i]>5000 && $mi[$i]<=8000)
        {
            my $l=$mi[$i]-800;
            $insert.=" -d$n 800 -a$n $mi[$i] -n$n $l ";
        }
        if($mi[$i]>8000)
        {
            my $l=$mi[$i]-1000;
            $insert.=" -d$n 1000 -a$n $mi[$i] -n$n $l ";
        }
    }
}
print SC "$Platanus scaffold -o $outdir/$prefix -b $outdir/$prefix\_contigBubble.fa -c $outdir/$prefix\_contig.fa $string $insert -t $ppn\n";
`perl $qsub --queue $queue --ppn $ppn $outdir/scaffold.sh`;