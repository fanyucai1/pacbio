#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
my $pbsuite="/home/fanyucai/software/PBSuite/PBSuite_15.8.24/";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
my $blasr="/home/Softwares/smrtanalysis_2.3.0/install/smrtanalysis_2.3.0.140936/analysis/bin/";
my $python="/home/fanyucai/software/python/Python-v2.7.9/bin/python";#pip install networkx
my $split="/home/fanyucai/software/FASTA_Splitter/fasta-splitter.pl";
my($contig,$outdir,$reads,$queue,$lines,$part,$ppn,$minGap,$minReads);
$outdir||=getcwd;
$queue||="all";
$part||=20;
$ppn=5;
$minGap||=10;
$minReads=2;
GetOptions(
    "s:s"=>\$contig,
    "o:s"=>\$outdir,
    "r:s"=>\$reads,
    "q:s"=>\$queue,
    "part:s"=>\$part,
    "ppn:s"=>\$ppn,
    "minGap:s"=>\$minGap,
           );
sub usage{
    print qq{
This script will build the scaffold using subreads.
links:https://sourceforge.net/p/pb-jelly/wiki/Home/?#058c
usage:
perl $0 -s scaffold.fa -o $outdir/pbjelly -r /path/to/subreads.fasta -q all -part 20 -ppn 5

options:
-c              scaffold sequence
-r              subreads(fasta)
-o              out directory(default:$outdir/pbjelly)
-q              which queue you run this shell script(default:all)
-part           Divide into <N> parts(default:20)
-ppn            cpu number(default:5)
-minGap         Minimum number of consecutive Ns to be considered a gap(default:10)
-minReads       Minimum number of reads required to fill a gap(default:2)

reference:
1:set minreads=2 #Improved maize reference genome with single-molecule technologies
2:set minGap=10   #The genome sequence of the Antarctic bullhead notothen reveals evolutionary adaptations to a cold environment
Email:fanyucai1\@126.com
2017.12.20
    };
    exit;
}
if(!$contig||!$reads)
{
    &usage();
}
system "mkdir -p $outdir/pbjelly/";
system "mkdir -p $outdir/pbjelly/data";
system "mkdir -p $outdir/pbjelly/data/reference";
system "mkdir -p $outdir/pbjelly/data/reads";
if(! -e "outdir/pbjelly/data/reference/contigs.fasta")
{
    system "ln -s $contig $outdir/pbjelly/data/reference/contigs.fasta";
}
#######################################################################split subreads
open(RUN1,">$outdir/pbjelly/fasta_spliter.sh");
print RUN1 "awk \'{print \$1}\' $reads >$outdir/pbjelly/subreads.fasta\n";
print RUN1 "perl $split --n-parts $part --out-dir $outdir/pbjelly/data/reads $outdir/pbjelly/subreads.fasta && rm $outdir/pbjelly/subreads.fasta\n";
$lines=`wc -l $outdir/pbjelly/fasta_spliter.sh`;
chomp($lines);
`perl $qsub --lines $lines $outdir/pbjelly/fasta_spliter.sh`;
##########################################################################prepare the subread and reference
open(PR,">$outdir/pbjelly/pre.sh");
#remove the " RQ=*" from the fastq read names, you can resume processing at the 'extraction' stage (https://sourceforge.net/p/pb-jelly/discussion/pbjtiks/thread/82457eda/)
print PR "source $pbsuite/setup.sh && fakeQuals.py $outdir/pbjelly/data/reference/contigs.fasta $outdir/pbjelly/data/reference/contigs.qual && cd $outdir/pbjelly/data/reference/ && sawriter contigs.fasta.sa contigs.fasta\n";
my @array=glob("$outdir/pbjelly/data/reads/*fasta");
foreach my $key(@array)
{
    my @array1=split(/.fasta/,$key);
    print PR "source $pbsuite/setup.sh && fakeQuals.py $key $array1[0].qual\n";
}
`perl $qsub --queue $queue $outdir/pbjelly/pre.sh`;
##########################################################################get Protocol.xml
open(XML,">$outdir/pbjelly/Protocol.xml");
print XML "<jellyProtocol>
    <reference>$outdir/pbjelly/data/reference/contigs.fasta</reference>
    <outputDir>$outdir/pbjelly/</outputDir>
    <cluster>
    <command>echo \'\${CMD}\'|qsub -q $queue -l nodes=1:ppn=$ppn -N pbjelly.\${JOBNAME} -o \${STDOUT} -e \${STDERR}</command>
    <nJobs>$part</nJobs>
    </cluster>
    <blasr>-minMatch 8 -sdpTupleSize 8 -minPctIdentity 75 -bestn 1 -nCandidates 10 -maxScore -500 -noSplitSubreads -nproc 20</blasr>
    <input baseDir=\"$outdir/pbjelly/data/reads/\">\n";
foreach my $key(@array)
{
    my $file=basename $key;
    print XML "<job>$key</job>\n";
    
}
print XML "</input>
    </jellyProtocol>";
my $run=0;
###################################################################################################setup
open(PB1,">$outdir/pbjelly/pbjelly1.sh");
print PB1 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB1 "$python $pbsuite/bin/Jelly.py setup $outdir/pbjelly/Protocol.xml -x \"--minGap=$minGap\"\n";
if(! -e "$outdir/setup.log")
{
        `sh $outdir/pbjelly/pbjelly1.sh`;
        sleep(15);
        &check("setup");
}
#####################################################################################################mapping
open(PB2,">$outdir/pbjelly/pbjelly2.sh");
print PB2 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB2 "$python $pbsuite/bin/Jelly.py mapping $outdir/pbjelly/Protocol.xml -x \"--nproc=20\"\n";
if(! -e "$outdir/mapping.log")
{
        `sh $outdir/pbjelly/pbjelly2.sh`;
        sleep(15);
        &check("mapping");
}
#####################################################################################################support
open(PB3,">$outdir/pbjelly/pbjelly3.sh");
print PB3 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB3 "$python $pbsuite/bin/Jelly.py support $outdir/pbjelly/Protocol.xml\n";
if( -e "$outdir/mapping.log")
{
    if(! -e "$outdir/support.log")
    {
        `sh $outdir/pbjelly/pbjelly3.sh`;
                sleep(15);
        &check("support");
    }
}
#####################################################################################################extraction
open(PB4,">$outdir/pbjelly/pbjelly4.sh");
print PB4 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB4 "$python $pbsuite/bin/Jelly.py extraction $outdir/pbjelly/Protocol.xml\n";
if(-e "$outdir/support.log")
{
    if(! -e "$outdir/extraction.log")
    {
        `sh $outdir/pbjelly/pbjelly4.sh`;
                sleep(15);
        &check("extraction");
    }
}
#####################################################################################################assembly
open(PB5,">$outdir/pbjelly/pbjelly5.sh");
print PB5 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB5 "$python $pbsuite/bin/Jelly.py assembly $outdir/pbjelly/Protocol.xml -x \"--nproc=20\"\n";
if(-e "$outdir/extraction.log")
{
    if(! -e "$outdir/assembly.log")
    {
        `sh $outdir/pbjelly/pbjelly5.sh`;
            sleep(15);
        &check("assembly");
    }
}
#####################################################################################################output
open(PB6,">$outdir/pbjelly/pbjelly6.sh");
print PB6 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB6 "$python $pbsuite/bin/Jelly.py output $outdir/pbjelly/Protocol.xml -x \"--minReads=$minReads\"\n";
if(-e "$outdir/assembly.log")
{
    if(! -e "$outdir/output.log")
    {
        `sh $outdir/pbjelly/pbjelly6.sh`;
        sleep(15);
        &check("output");
    }
}

sub check
{
    my $process=$_[0];
    my $user = `whoami`;
    chomp $user;
    while(1)
    {
        my $qstat=`qstat -u $user`;
        my @jobs = split("\n", $qstat);
        my $m=0;
        foreach my $key(@jobs)
        {
            if($key=~/pbjelly/i)
            {
                $m=1;
            }
        }
        last if($m==0);
        sleep(10);
    }
    system "echo '$process done'>$outdir/$process.log";
}
