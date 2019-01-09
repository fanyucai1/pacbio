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
my($contig,$outdir,$reads,$queue,$lines);
$outdir||=getcwd;
$queue||="all";
GetOptions(
    "s:s"=>\$contig,
    "o:s"=>\$outdir,
    "r:s"=>\$reads,
    "q:s"=>\$queue,
           );
sub usage{
    print qq{
This script will build the scaffold using subreads.
links:https://sourceforge.net/p/pb-jelly/wiki/Home/?#058c
usage:
perl $0 -s scaffold.fa -o $outdir/pbjelly -r /path/to/subreads.fasta -q all 

options:
-c              scaffold sequence
-r              subreads(fasta)
-o              outdirectory(default:$outdir/pbjelly)
-q              which queue you run this shell script(default:all)
-split          Divide into <N> parts(default:20)

Email:fanyucai1\@126.com
2017.7.18
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
system "ln -s $contig $outdir/pbjelly/data/reference/contigs.fasta";

##########################################################################
open(PR,">$outdir/pbjelly/pre.sh");
#remove the " RQ=*" from the fastq read names, you can resume processing at the 'extraction' stage (https://sourceforge.net/p/pb-jelly/discussion/pbjtiks/thread/82457eda/)
print PR "awk \'{print \$1}\' $reads >$outdir/pbjelly/data/reads/subreads.fasta\n";
print PR "source $pbsuite/setup.sh && fakeQuals.py $outdir/pbjelly/data/reference/contigs.fasta $outdir/pbjelly/data/reference/contigs.qual\n";
print PR "source $pbsuite/setup.sh && fakeQuals.py $outdir/pbjelly/data/reads/subreads.fasta $outdir/pbjelly/data/reads/subreads.qual\n";
print PR "cd $outdir/pbjelly/data/reference/ && source $pbsuite/setup.sh && sawriter contigs.fasta.sa contigs.fasta\n";
$lines=`wc -l $outdir/pbjelly/pre.sh`;
chomp($lines);
`perl $qsub --lines $lines --queue $queue $outdir/pbjelly/pre.sh`;
##########################################################################
open(XML,">$outdir/pbjelly/Protocol.xml");
print XML "<jellyProtocol>
    <reference>$outdir/pbjelly/data/reference/contigs.fasta</reference>
    <outputDir>$outdir/pbjelly/</outputDir>
    <blasr>-minMatch 8 -sdpTupleSize 8 -minPctIdentity 75 -bestn 1 -nCandidates 10 -maxScore -500 -noSplitSubreads -nproc=40</blasr>
    <input baseDir=\"$outdir/pbjelly/data/reads/\">
    <job>subreads.fasta</job>
   </input>
    </jellyProtocol>";
###################################################################################################
open(PB1,">$outdir/pbjelly/pbjelly1.sh");
print PB1 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB1 "$python $pbsuite/bin/Jelly.py setup $outdir/pbjelly/Protocol.xml -x \"--minGap=1\"\n";
`$qsub $outdir/pbjelly/pbjelly1.sh`;


open(PB2,">$outdir/pbjelly/pbjelly2.sh");
print PB2 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB2 "$python $pbsuite/bin/Jelly.py mapping $outdir/pbjelly/Protocol.xml -x \"--nproc=20\"\n";
`$qsub --ppn 20  $outdir/pbjelly/pbjelly2.sh`;


open(PB3,">$outdir/pbjelly/pbjelly3.sh");
print PB3 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB3 "$python $pbsuite/bin/Jelly.py support $outdir/pbjelly/Protocol.xml\n";
`$qsub $outdir/pbjelly/pbjelly3.sh`;

open(PB4,">$outdir/pbjelly/pbjelly4.sh");
print PB4 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB4 "$python $pbsuite/bin/Jelly.py extraction $outdir/pbjelly/Protocol.xml\n";
`$qsub $outdir/pbjelly/pbjelly4.sh`;

open(PB5,">$outdir/pbjelly/pbjelly5.sh");
print PB5 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB5 "$python $pbsuite/bin/Jelly.py assembly $outdir/pbjelly/Protocol.xml -x \"--nproc=20\"\n";
`$qsub --ppn 20 $outdir/pbjelly/pbjelly5.sh`;


open(PB6,">$outdir/pbjelly/pbjelly6.sh");
print PB6 "export PATH=$pbsuite/bin:$blasr:\$PATH && source $pbsuite/setup.sh && ";
print PB6 "$python $pbsuite/bin/Jelly.py output $outdir/pbjelly/Protocol.xml\n";
`$qsub $outdir/pbjelly/pbjelly6.sh`;
