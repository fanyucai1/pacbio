#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use Cwd;
my $env="export LD_LIBRARY_PATH=/home/Softwares/gcc-5.2.0/lib64/:\$LD_LIBRARY_PATH";
my $MECAT="/home/fanyucai/software/MECAT/Linux-amd64/bin/";
my $qsub="/home/fanyucai/software/qsub/qsub-pbs.pl";
my ($input,$outdir,$genomesize,$threads,$num,$queue,$a,$l);
$threads||=20;
$queue||="fat";
$a||=2000;
$l=5000;
$num||=100;
$outdir||=getcwd;
GetOptions(
    "i:s"=>\$input,
    "o:s"=>\$outdir,
    "g:s"=>\$genomesize,
    "t:s"=>\$threads,
    "a:s"=>\$a,
    "l:s"=>\$l,
    "n:s"=>\$num,
    "q:s"=>\$queue,
           );
sub usage{
    print qq{
This scirpt will assembly the Pacbio data use MECAT.
usage:
perl $0 -i input -o $outdir -g 4700000 -n 50 -t 16 -q fat
options:
-i          input file(fastq or fasta)
-o          dir of output
-t          thread number(default:20)
-g          the genome size
-a          minimum overlap size(default:2000)
-l          minimum length of corrected sequence(default:5000)
-n          Number of candidates considered for gapped extension:
                For GS(genome size) < 20M, NC should be set as 200;
                For GS>20M and GS<200M; NC should be set as 100;
                For GS>200M, NC should be set as 50.
-q          run shell in queue(default:fat)

Email:fanyucai1\@126.com
2017.5.27
    };
    exit;
}

if(!$input ||!$queue||!$genomesize)
{
    &usage();
}
system "mkdir -p $outdir";
open(RUN1,">$outdir/MECAT1.sh");#Step 1, using mecat2pw to detect overlapping candidates
print RUN1 "cd $outdir && $env && $MECAT/mecat2pw -j 0 -n $num -d $input -o mecat.pm.can -w $outdir/ -t $threads\n";

open(RUN2,">$outdir/MECAT2.sh");#Step 2, correct the noisy reads based on their pairwise overlapping candidates.
print RUN2 "cd $outdir &&  $env && $MECAT/mecat2cns -a $a -l $l -i 0 -t $threads mecat.pm.can $input corrected_mecat\_filtered.fasta\n";

open(RUN3,">$outdir/MECAT3.sh");#Step 3, extract the longest 25X corrected reads
system "mkdir -p $outdir/25x";
system "mkdir -p $outdir/30x";
system "mkdir -p $outdir/35x";
system "mkdir -p $outdir/40x";
system "mkdir -p $outdir/50x";
system "mkdir -p $outdir/60x";
system "mkdir -p $outdir/all";
print RUN3 "cd $outdir/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/extract_sequences corrected_mecat\_filtered.fasta $outdir/25x/corrected_mecat\_25x.fasta $genomesize 25\n";
print RUN3 "cd $outdir/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/extract_sequences corrected_mecat\_filtered.fasta $outdir/30x/corrected_mecat\_30x.fasta $genomesize 30\n";
print RUN3 "cd $outdir/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/extract_sequences corrected_mecat\_filtered.fasta $outdir/40x/corrected_mecat\_40x.fasta $genomesize 40\n";
print RUN3 "cd $outdir/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/extract_sequences corrected_mecat\_filtered.fasta $outdir/50x/corrected_mecat\_50x.fasta $genomesize 50\n";
print RUN3 "cd $outdir/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/extract_sequences corrected_mecat\_filtered.fasta $outdir/80x/corrected_mecat\_60x.fasta $genomesize 60\n";


open(RUN4,">$outdir/MECAT4.sh");#Step 4, assemble the longest 25X corrected reads using mecat2cacu
print RUN4 "cd $outdir/25x/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/mecat2canu -trim-assemble -p mecat -d mecat genomeSize=$genomesize ErrorRate=0.02 maxMemory=200 maxThreads=$threads useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected corrected_mecat\_25x.fasta.fasta\n";
print RUN4 "cd $outdir/30x/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/mecat2canu -trim-assemble -p mecat -d mecat genomeSize=$genomesize ErrorRate=0.02 maxMemory=200 maxThreads=$threads useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected corrected_mecat\_30x.fasta.fasta\n";
print RUN4 "cd $outdir/40x/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/mecat2canu -trim-assemble -p mecat -d mecat genomeSize=$genomesize ErrorRate=0.02 maxMemory=200 maxThreads=$threads useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected corrected_mecat\_40x.fasta.fasta\n";
print RUN4 "cd $outdir/50x/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/mecat2canu -trim-assemble -p mecat -d mecat genomeSize=$genomesize ErrorRate=0.02 maxMemory=200 maxThreads=$threads useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected corrected_mecat\_50x.fasta.fasta\n";
print RUN4 "cd $outdir/60x/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/mecat2canu -trim-assemble -p mecat -d mecat genomeSize=$genomesize ErrorRate=0.02 maxMemory=200 maxThreads=$threads useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected corrected_mecat\_60x.fasta.fasta\n";
print RUN4 "cd $outdir/all/ &&  export PATH=$MECAT:\$PATH && $env && $MECAT/mecat2canu -trim-assemble -p mecat -d mecat genomeSize=$genomesize ErrorRate=0.02 maxMemory=200 maxThreads=$threads useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected corrected_mecat\_filtered.fasta\n";

=head
system "perl $qsub --queue $queue $outdir/MECAT1.sh";
system "perl $qsub --queue $queue $outdir/MECAT2.sh";
system "ln -s $outdir/corrected_mecat\_filtered.fasta $outdir/25x/";
system "ln -s $outdir/corrected_mecat\_filtered.fasta $outdir/30x/";
system "ln -s $outdir/corrected_mecat\_filtered.fasta $outdir/40x/";
system "ln -s $outdir/corrected_mecat\_filtered.fasta $outdir/50x/";
system "ln -s $outdir/corrected_mecat\_filtered.fasta $outdir/60x/";
system "perl $qsub --queue $queue $outdir/MECAT3.sh";
system "perl $qsub --queue $queue $outdir/MECAT4.sh";
=cut
