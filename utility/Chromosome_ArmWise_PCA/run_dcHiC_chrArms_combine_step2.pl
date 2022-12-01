#!/usr/bin/env perl

use File::Basename;
use Getopt::Long;

 GetOptions(
 'cmd=s' => \my $cmd,
);

sub help {
  print("
HELP:
\t--cmd=input.txt file\n\n");
  exit;
}

sub error {
  my $msg = $_[0];
  print ("$msg\n");
  help();
}

sub findMax {
  my $fileName = $_[0];
  my $dirName  = $_[1];
  my $max = 0;
  open(PCFILE, "$dirName/$fileName");
  while (my $value = <PCFILE>) {
    chomp $value;
    my $v = (split(/\s+/,$value))[3];
    if (abs($v) > $max) {
      $max = abs($v);
    } 
  }
  close PCFILE;
  return $max;
}

sub writeFile {
  my $chrName = $_[0];
  my $dirName = $_[1];
  my $repName = $_[2];
  my $maxParm = $_[3];
  my $append  = $_[4];

  system("mkdir -p $repName\_pca/intra_pca/$repName\_mat");
  if ($append == 0) {
    -e "$repName\_pca/intra_pca/$repName\_mat/$chrName" ? system("$repName\_pca/intra_pca/$repName\_mat/$chrName") : 0;
    open(OUT, ">$repName\_pca/intra_pca/$repName\_mat/$chrName");
  } elsif ($append == 1) {
    open(OUT, ">>$repName\_pca/intra_pca/$repName\_mat/$chrName");
  }

  open(PCFILE, "$dirName/$chrName");
  while (my $bdg = <PCFILE>) {
    chomp $bdg;
    my $c = (split(/\s+/,$bdg))[0];
    my $s = (split(/\s+/,$bdg))[1];
    my $e = (split(/\s+/,$bdg))[2];
    my $v = (split(/\s+/,$bdg))[3];
    
    my $scaled = sprintf("%.5f", $v/$maxParm);
    print OUT ("$c\t$s\t$e\t$scaled\n");
  }
  close PCFILE;
  close OUT;
}

chomp ($cpu, $mem, $hrs, $pfx, $cmd, $cen);

$cmd eq "" ?  error("Please provide job COMMAD (--cmd) file to be used! existing") : print ("COMMAND = $cmd\n");

open(FILE, $cmd);
  while (my $line = <FILE>) {
    chomp $line;
    print ("$line\n");
    my $mat = (split(/\s+/,$line))[0];
    my $bed = (split(/\s+/,$line))[1];
    my $rep = (split(/\s+/,$line))[2];
    my $spl = (split(/\s+/,$line))[3];

    my @parm_pcs = `ls parm/$rep\_pca/intra_pca/$rep\_mat/*.pc.bedGraph`;
    my @qarm_pcs = `ls qarm/$rep\_pca/intra_pca/$rep\_mat/*.pc.bedGraph`;

    foreach my $parm_pc (@parm_pcs) {
      chomp $parm_pc;
      my $max_parm = findMax(basename($parm_pc), dirname($parm_pc));
      writeFile(basename($parm_pc), dirname($parm_pc), $rep, $max_parm, 0);
    }

    foreach my $qarm_pc (@qarm_pcs) {
      chomp $qarm_pc;
      my $max_qarm = findMax(basename($qarm_pc), dirname($qarm_pc));
      writeFile(basename($qarm_pc), dirname($qarm_pc), $rep, $max_qarm, 1);
    }
  }
  close FILE; 
