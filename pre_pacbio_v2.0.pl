#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename;
use Getopt::Long;
use Cwd;
use Bio::SeqIO;
use Bio::Seq;
use List::Util qw(max min sum);
my $R="/home/fanyucai/software/R/R-v3.2.0/bin/Rscript";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
my $env="export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:\$LD_LIBRARY_PATH && export PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/bin/:\$PATH";
my $bamtools="/smrtlinks/smrtlink/smrtcmds/bin/bamtools";

my (@subreads,$outdir,@prefix,@bam);
$outdir||=getcwd;
GetOptions(
    "s:s{1,}"=>\@subreads,
    "o:s"=>\$outdir,
    "p:s{1,}"=>\@prefix,
    "bam:s{1,}"=>\@bam,
           );
sub usage{
    print qq{
 This script will statistics subreads from PacBio.
 usage:
 perl $0 -s m1.subreads.fasta m2.subreads.fasta -o $outdir -p m1 m2
            or
 perl $0 -bam m1.subreads.bam m2.subreads.bam  -o $outdir -p m1 m2
 options:
 -s         subreads (fasta):several split by space
 -o         output directory(default:$outdir)
 -p         output prefix:several split by space
 -bam       subreads bam:several split by space
 
 Email:fanyucai1\@126.com
 2018.1.19
    };
    exit;
}
if(!@prefix)
{
    &usage();
}
system "mkdir -p $outdir";
open(HIST,">$outdir/hist.Rscript");
open(STAT,">$outdir/length_stat.xls");
print STAT "#Sample_ID\tTotal_reads\tTotal_base(bp)\tMin_length\tMax_length\tMean_length\n";
if(@bam)
{
    open(BAM,">$outdir/bam2fa.sh");
    for(my $i=0;$i<=$#prefix;$i++)
    {
        print BAM "cd $outdir && $bamtools convert -format fasta -in $bam[$i] -out $outdir/$prefix[$i].fasta\n";
    }
    system "perl $qsub $outdir/bam2fa.sh";
}
for(my $i=0;$i<=$#prefix;$i++)
{
     $subreads[$i]="$outdir/$prefix[$i].fasta";
}
if(@subreads)
{
    for(my $i=0;$i<=$#prefix;$i++)
    {
        open(LEN,">$outdir/$prefix[$i].len_dis.txt");
        my $seqin= Bio::SeqIO->new( -format => 'Fasta', -file => "$subreads[$i]");
        my (@seqlen,$total,$num);
        while((my $seqobj=$seqin->next_seq()))
        {
            $num++;
            push @seqlen,$seqobj->length();
            print LEN $seqobj->length(),"\n";
        }
        my $mean=sum(@seqlen)/$num;
        print STAT $prefix[$i],"\t",$num,"\t",sum(@seqlen),"\t",min(@seqlen),"\t",max(@seqlen),"\t",$mean,"\n";
        print HIST "
        #!$R
        ################################
        a=read.table(\"$outdir/$prefix[$i].len_dis.txt\",header=F)
        library(ggplot2)
        p=ggplot(a,aes(a[,1]))+geom_histogram(binwidth =1000,fill=\"steelblue\",colour=\"black\")+ylab(\"Count\")+xlab(\"Length(bp)\")+ggtitle(\"$prefix[$i]:The distribution of sequence\")
        p=p+xlim(0, 15000)
        pdf(\"$outdir/$prefix[$i].pdf\",width=12,height=8)
        p
        dev.off()
        png(\"$outdir/$prefix[$i].png\",width=1500,height=1200,res=300)
        p
        dev.off()
        ";
    }
}
###############################################
system "cd $outdir && $R $outdir/hist.Rscript";
system "rm -rf $outdir/*.len_dis.txt";
