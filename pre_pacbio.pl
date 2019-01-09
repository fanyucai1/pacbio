#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename;
use Getopt::Long;
use Cwd;
use List::Util qw(max min sum);
my $R="/home/fanyucai/software/R/R-v3.2.0/bin/Rscript";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
my $env="export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:\$LD_LIBRARY_PATH && export PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/bin/:\$PATH";
my $bam2fasta="/smrtlinks/smrtlink/smrtcmds/bin/bam2fasta";
my ($subreads,$outdir,$bin,$prefix,$type,$bam);
$outdir||=getcwd;
$bin||=1000;
$type||="";
GetOptions(
    "s:s"=>\$subreads,
    "o:s"=>\$outdir,
    "b:s"=>\$bin,
    "p:s"=>\$prefix,
    "t:s"=>\$type,
    "bam:s"=>\$bam,
           );
sub usage{
    print qq{
 This script will statistics subreads from PacBio.
 usage:
 perl $0 -s subreads.fasta -o $outdir -b 1000 -p prefix -t ISO-seq
            or
 perl $0 -bam subreads.bam -o $outdir -b 1000 -p prefix
 options:
 -s         subreads (fasta)
 -o         outputdirectory(default:$outdir)
 -b         plot length distribution of bin size(default:1000)
 -p         prefix
 -bam       subreads bam
 -t         if you defined this (ISO-seq) will plot special
 
 Email:fanyucai1\@126.com
 2017.7.20
    };
    exit;
}
if(!$subreads)
{
    &usage();
}
system "mkdir -p $outdir/";
if($subreads)
{
    open(FA,"$subreads");
}
if($bam)
{
    system "echo 'cd $outdir && $bam2fasta -o $outdir/subreads -u $bam'>$outdir/bam2fasta.sh";
    `perl $qsub $outdir/bam2fasta.sh`
}
my ($name,%hash,$num);
while(<FA>)
{
    chomp;
    if($_=~/^\>/)
    {
        $name=$_;
        $num++;
    }
    else
    {
        $hash{$name}.=$_;
    }
}
###########################
my (%dis,$num500,$sum500,$number,@all,%hash2,$other);
open(OUT,">$outdir/$prefix\_lenth_dis.txt");
print OUT "Read_length\tCounts\n";
$number=-1;
$other=0;
TT:foreach my $key(keys %hash)
{
    $number++;
    $all[$number]=length($hash{$key});
    if(length($hash{$key})>=500)
    {
        $num500++;
        $sum500+=length($hash{$key});
    }
    if($all[$number]>15000)
    {
        $other++;
    }
    for(my $i=1;$i<=15;$i++)
    {
        if($all[$number]<=$i*1000 && $all[$number]>($i-1)*1000)
        {
            my $string=1000*($i-1)+1;
            $string.="_";
            $string.=1000*$i;
            $hash2{$string}++;
            next TT;
        }
    }
}
foreach my $key(sort {$a <=> $b} keys %hash2)
{
    print OUT $key,"\t",$hash2{$key},"\n";
}
print OUT "(>15000)\t",$other,"\n";
##########################
my $seqmin=min (@all);
my $seqmax=max (@all);
my $seqsum=sum (@all);
my $mean=$seqsum/$num;


system "echo '#!$R
a=read.table(\"$outdir/$prefix\_lenth_dis.txt\",header=T,sep=\"\\t\")
library(ggplot2)
p=ggplot(a,aes(x=factor(Read_length,order=T,levels=Read_length),y=Counts))+geom_bar(stat= \"identity\", fill = \"steelblue\")
p=p+geom_text(aes(label =Counts), vjust = 0.2, size = 2.5, position = position_dodge(0.9))
p=p+xlab(\"Read length (bp)\")+theme_bw()+theme(panel.grid =element_blank(),axis.text.x=element_text(angle=45,vjust=0.5))
png(\"$outdir/$prefix\_length_dis.png\",res=300,width=2500,height=2500)
p
dev.off()
pdf(\"$outdir/$prefix\_length_dis.pdf\",width=9,height=9)
p
dev.off()
'>$outdir/$prefix\_len_dis.Rscript";
system "$env && $R $outdir/$prefix\_len_dis.Rscript";
###################################
open(ST,">$outdir/$prefix\_stat.xls");
$sum500=$sum500/$num500;
print ST "#Total_reads\tTotal_base(bp)\tMin\tMax\tMean\n";
print ST "$num\t$seqsum\t$seqmin\t$seqmax\t$mean\n";
