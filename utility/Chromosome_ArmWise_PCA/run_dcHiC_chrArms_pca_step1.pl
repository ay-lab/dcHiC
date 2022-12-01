#!/usr/bin/env perl

use File::Basename;
use Getopt::Long;

my $DCHIC = "../../dchicf.r";
chomp $DCHIC; 

 GetOptions(
 'cpu=s' => \my $cpu,
 'mem=s' => \my $mem,
 'hrs=s' => \my $hrs,
 'pfx=s' => \my $pfx,
 'cmd=s' => \my $cmd,
 'cen=s' => \my $cen,
 'pexcl=s' => \my $pexcl,
 'qexcl=s' => \my $qexcl
);

sub help {
  print("
HELP:
\t--cpu=Number of CPUs to be used for the job
\t--mem=Amount of memory to be used for the job (in GB)
\t--hrs=Amount of wall time requested for the job (in hours)
\t--pfx=Prefix of job name
\t--cmd=input.txt file
\t--pexcl=Exclude chromosomes from parm pca analysis (provide multiple chromosome names as e.g. chr21,chr22)
\t--qexcl=Exclude chromosomes from qarm pca analysis (provide multiple chromosome names as e.g. chr21,chr22)
\t--cen=centromere file (Downlod centromere file from https://www.ncbi.nlm.nih.gov/grc/human)\n\n");
  exit;
}

sub error {
  my $msg = $_[0];
  print ("$msg\n");
  help();
}

chomp ($cpu, $mem, $hrs, $pfx, $cmd, $cen, $pexcl, $qexcl);

$cpu eq "" ?  error("Please provide CPU (--cpu) to be used! existing") : print ("CPU = $cpu\n");
$mem eq "" ?  error("Please provide MEMORY (--mem) to be used! existing") : print ("MEMORY = $mem\n");
$hrs eq "" ?  error("Please provide HOURS (--hrs) to be used! existing") : print ("HOURS = $hrs\n");
$pfx eq "" ?  error("Please provide job PREFIX (--pfx) name to be used! existing") : print ("PREFIX = $pfx\n");
$cmd eq "" ?  error("Please provide job COMMAD (--cmd) file to be used! existing") : print ("COMMAND = $cmd\n");
$cen eq "" ?  error("Please provide job CENTROMERE (--cen) file to be used! existing") : print ("CENTROMERE = $cen\n");

my %pexclude;
my %qexclude;
my @pexcludeChrs = split(/,/,$pexcl);
my @qexcludeChrs = split(/,/,$qexcl);

map {
 chomp $_;
 $pexclude{$_} = 1;
} @pexcludeChrs;

map {
 chomp $_;
 $qexclude{$_} = 1;
} @qexcludeChrs;


my %centromere;
open(CEN, $cen);
while (my $line = <CEN>) {
  chomp $line;
  if ($. > 1) {
    my $chr   = "chr".(split(/\s+/,$line))[1];
    my $start = (split(/\s+/,$line))[2];
    my $end   = (split(/\s+/,$line))[3];
    $centromere{$chr}{'start'} = $start;
    $centromere{$chr}{'end'}   = $end;
    print ("$chr\t$centromere{$chr}{'start'}\t$centromere{$chr}{'end'}\n");
  }
}
close CEN;

-d "parm" ? print ("parm folder exists\n") : system("mkdir -p parm/data");
-d "qarm" ? print ("qarm folder exists\n") : system("mkdir -p qarm/data");

open(PINPUT, ">parm/input_parm.txt");
open(QINPUT, ">qarm/input_qarm.txt");
open(FILE, $cmd);
  while (my $line = <FILE>) {
    chomp $line;
    print ("$line\n");
    my $mat = (split(/\s+/,$line))[0];
    my $bed = (split(/\s+/,$line))[1];
    my $rep = (split(/\s+/,$line))[2];
    my $spl = (split(/\s+/,$line))[3];


    my %index2coord;
    my @parm_coord;
    my @qarm_coord;
    open(POUT, ">parm/data/$rep\_parm.bed");
    open(QOUT, ">qarm/data/$rep\_qarm.bed");
    open(BED, $bed);
    while (my $bed_line = <BED>) {
      chomp $bed_line;
      my $chr   = (split(/\s+/,$bed_line))[0];
      my $start = (split(/\s+/,$bed_line))[1];
      my $end   = (split(/\s+/,$bed_line))[2];
      my $index = (split(/\s+/,$bed_line))[3];
      $index2coord{$index}{'chr'}   = $chr;
      $index2coord{$index}{'start'} = $start;
      $index2coord{$index}{'end'}   = $end;

      if ($end < $centromere{$chr}{'start'}) {
        push (@parm_coord, $bed_line);
      } elsif ($start > $centromere{$chr}{'end'}) {
        push (@qarm_coord, $bed_line);
      }
    }
    close BED;

    %armwisecov;
    open(MAT, $mat);
    while (my $mat_line = <MAT>) {
      chomp $mat_line;
      my $indexA = (split(/\s+/,$mat_line))[0];
      my $indexB = (split(/\s+/,$mat_line))[1];
      my $chrA   = $index2coord{$indexA}{'chr'};
      my $chrB   = $index2coord{$indexB}{'chr'};
      my $startA = $index2coord{$indexA}{'start'};
      my $startB = $index2coord{$indexB}{'start'};
      my $endA   = $index2coord{$indexA}{'end'};
      my $endB   = $index2coord{$indexB}{'end'};

      if ($chrA eq $chrB) {
        if ($endA < $centromere{$chrA}{'start'} && $endB < $centromere{$chrA}{'start'}) {
          $armwisecov{'parm'}{$chrA} += 1;
        } elsif ($startA > $centromere{$chrA}{'end'} && $startB > $centromere{$chrA}{'end'}) {
          $armwisecov{'qarm'}{$chrA} += 1;
        }
      }
    }
    close MAT;

    foreach my $pcoord (@parm_coord) {
      chomp $pcoord;
      my $chr   = (split(/\s+/,$pcoord))[0];
      my $start = (split(/\s+/,$pcoord))[1];
      my $end   = (split(/\s+/,$pcoord))[2];
      my $index = (split(/\s+/,$pcoord))[3];
      if ($armwisecov{'parm'}{$chr} >= 50 && $pexclude{$chr} eq "") {
        print POUT ("$pcoord\n"); 
      } 
    }

    foreach my $qcoord (@qarm_coord) {
      chomp $qcoord;
      my $chr   = (split(/\s+/,$qcoord))[0];
      my $start = (split(/\s+/,$qcoord))[1];
      my $end   = (split(/\s+/,$qcoord))[2];
      my $index = (split(/\s+/,$qcoord))[3];
      if ($armwisecov{'qarm'}{$chr} >= 50 && $qexclude{$chr} eq "") {
        print QOUT ("$qcoord\n"); 
      } 
    }

    close POUT;
    close QOUT;

    my $mat_abs_path = `readlink -f $mat`;
    chomp $mat_abs_path;

    system("ln -fs $mat_abs_path ./parm/data/");
    system("ln -fs $mat_abs_path ./qarm/data/");
   
    print PINPUT ("./data/".basename($mat_abs_path),"\t./data/$rep\_parm.bed\t$rep\t$spl\n");
    print QINPUT ("./data/".basename($mat_abs_path),"\t./data/$rep\_qarm.bed\t$rep\t$spl\n");
  }
  close FILE; 
  close PINPUT;
  close QINPUT;

my @job_param ="#!/bin/bash -ex
#PBS -l nodes=1:ppn=$cpu
#PBS -l mem=$mem\GB
#PBS -l walltime=$hrs:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
hostname
TMPDIR=/scratch\n";

chdir "parm";
$pwd = `pwd`;
chomp $pwd;
open(OUT, ">parm_$pfx\_job.sh");
print OUT @job_param;
print OUT "cd $pwd\n\n";
print OUT ("Rscript $DCHIC --file input_parm.txt --pcatype cis --pc 4 --dirovwt T\n");
print OUT ("Rscript $DCHIC --file input_parm.txt --pcatype select --pc 4 --dirovwt T --genome hg38\n");
close OUT;
chdir "../";


chdir "qarm";
$pwd = `pwd`;
chomp $pwd;
open(OUT, ">qarm_$pfx\_job.sh");
print OUT @job_param;
print OUT "cd $pwd\n\n";
print OUT ("Rscript $DCHIC --file input_qarm.txt --pcatype cis --pc 4 --dirovwt T\n");
print OUT ("Rscript $DCHIC --file input_qarm.txt --pcatype select --pc 4 --dirovwt T --genome hg38\n");
close OUT;
chdir "../";
