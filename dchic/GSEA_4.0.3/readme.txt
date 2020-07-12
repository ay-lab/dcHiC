=======================
GSEA Binary Distribution
=======================

Prerequisites:

Java 11 (http://openjdk.java.net).  This is bundled with our distributions.
Not compatible with Java 8, 9, 10.


Instructions:

1. Download and unzip the distribution file to a directory of your choice.

2. Note that the Java command line has become more complex with Java 11 compared to
Java 8.  You can start GSEA with one of the following scripts; this is the recommended 
way to launch.  Some of these may not be present depending on the distribution you
downloaded.  You might have to make the script executable (chmod a+x gsea.sh).


gsea.bat       (for Windows)
gsea.sh        (for Linux and macOS)
gsea_hidpi.sh  (for Linux with HiDPI screens)
gsea-cli.sh    (for Linux and macOS command line usage)
gsea-cli.bat   (for Windows command line usage)
gsea.command   (for macOS, double-click to start)

The bat and shell scripts are configured to start GSEA with 4GB of
memory.  This is a reasonable default for most machines.  If you are
working with very large datasets you can increase the amount of memory
available to GSEA by editing the first line of the startup script.
Specifically change the value of the "-Xmx" parameter.  For example,
to start GSEA with 8 GB of memory, change the value

   -Xmx4g

to

   -Xmx8g

Use gsea-cli.sh (Linux, Mac) or gsea-cli.bat (Windows) to run GSEA at the command line,
For example:
     $ gsea-cli.sh GSEA [parameters]
Or, more generally
     $ gsea-cli.sh [operationName] [parameters]
     
Where [operationName] is one of GSEA, GSEAPreranked, CollapseDataset, Chip2Chip, or LeadingEdgeTool
for the chosen operation and [parameters] are the corresponding parameters for that operation.
Use of the Command feature in the UI is the best way to discover the available parameters 
for each operation.

Alternatively, to start the GSEA Desktop UI execute the following from the command line:

     $ java --module-path=modules -Xmx4g @gsea.args --patch-module=jide.common=lib/jide-components-3.7.4.jar:lib/jide-dock-3.7.4.jar:lib/jide-grids-3.7.4.jar --module=org.gsea_msigdb.gsea/xapps.gsea.GSEA

As noted above, use of the scripts above instead is strongly recommended.