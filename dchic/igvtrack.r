htmlheader <- function() {
cat ("<!DOCTYPE html>
<html lang=\"en\">
<head>
    <meta charset=\"utf-8\">
    <meta content=\"IE=edge\" http-equiv=\"X-UA-Compatible\">
    <meta content=\"width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no\" name=\"viewport\">
    <meta content=\"\" name=\"description\">
    <meta content=\"\" name=\"",paste0(Sys.info()["user"][[1]]),"\">
    <link href=https://igv.org/web/release/2.5.5/examples/img/favicon.ico rel=\"igv icon\">
    <title>IGV - Dev</title>
    <script src=\"https://cdn.jsdelivr.net/npm/igv@2.5.5/dist/igv.min.js\"></script>
</head>

<body>
<h1><u>dcHiC: Differential Compartment calling from Hi-C data</u></h1>
<p style=\"color: #C86400\"><b>A compartment</b></p>
<p style=\"color: #0064C8\"><b>B compartment</b></p>
<p style=\"color: #99FFCC\"><b>Mdist: Mahalanobis distance score to represent outlierness of the bin</b></p>
<p style=\"color: #009900\"><b>log10Padj: P adjusted value of corresponding Mdist score</b></p>
<p style=\"color: #994C00\"><b>dZsc: Distance Zscore that represent biological variability as compared to technical variability</b></p>
<div id=\"igvDiv\" style=\"padding-top: 10px;padding-bottom: 10px; border:1px solid lightgray\"></div>\n")
}

scriptbody <- function(genome,start=F,end=F) {

  if (start) {
  cat ("<script type=\"text/javascript\">
    document.addEventListener(\"DOMContentLoaded\", function () {
      var igvDiv = document.getElementById(\"igvDiv\");
      var options = {
        locus: '19:49301000-49305700',
        genome:",paste0("\"",genome,"\""),",
	tracks: [\n")
  }
  if (end) {
  cat ("
       ]
     }; igv.createBrowser(igvDiv, options)
                .then(function (browser) {
                    console.log(\"Created IGV browser\")});
   })
  </script>
 </body>
</html>\n")
  }
}

htmlbody <- function(file, folder, genome) {
 
  dir.create(paste0(folder,"/data"))
  scriptbody(genome,start=T)
  
  for(i in 1:nrow(file)) {
    if (file$group[i] == "compartment") {

      compbdg <- read.table(file$file[i], h=T, as.is=T, comment.char = "")
      colnames(compbdg)[1] <- "chr"

      ## Write the PC coordinates
      for(j in 4:(ncol(compbdg)-4)) {
cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
track name=\"",colnames(compbdg)[j]," PC1\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20 plotType=\"points\"\n",
file=paste0(folder,"/data/",colnames(compbdg)[j],".PC.bedGraph"))

      write.table(compbdg[,c(1:3,j)], file=paste0(folder,"/data/",colnames(compbdg)[j],".PC.bedGraph"), row.names=F, col.names=F,sep="\t", quote=F, append=T)
      R.utils::gzip(paste0(folder,"/data/",colnames(compbdg)[j],".PC.bedGraph"), paste0(folder,"/data/",colnames(compbdg)[j],".PC.bedGraph.gz"), overwrite=T)
      cmd <- paste0("create_datauri ",folder,"/data/",colnames(compbdg)[j],".PC.bedGraph.gz")
      datauri <- system(cmd, wait=T, intern=T)
        cat ("
             {
                name: ",paste0("\'",colnames(compbdg)[j],".PC',"),"
                url: ",paste0("\"",datauri,"\""),",
                indexed: false,
                format: \"bedGraph\"
             },\n")       
      }
      ## Write Mdist values
trackmaster <- paste0(colnames(compbdg)[4:(ncol(compbdg)-4)],collapse="_vs_")
cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
track name=\"",trackmaster," Mdist\" description=\"BedGraph format\" visibility=full color=153,255,204 priority=20 plotType=\"points\"\n",
file=paste0(folder,"/data/",trackmaster,".Mdist.bedGraph"))

      write.table(compbdg[,c(1:3,(ncol(compbdg)-3))], file=paste0(folder,"/data/",trackmaster,".Mdist.bedGraph"), row.names=F, col.names=F,sep="\t", quote=F, append=T)
      R.utils::gzip(paste0(folder,"/data/",trackmaster,".Mdist.bedGraph"), paste0(folder,"/data/",trackmaster,".Mdist.bedGraph.gz"), overwrite=T)    
      cmd <- paste0("create_datauri ",folder,"/data/",trackmaster,".Mdist.bedGraph.gz")
      datauri <- system(cmd, wait=T, intern=T)
        cat ("
            {
               name: ",paste0("\'",trackmaster,".Mdist',"),"
               url: ",paste0("\"",datauri,"\""),",
               indexed: false,
               format: \"bedGraph\"
            },\n")
  
     ## Write Padjusted values
cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
track name=\"",trackmaster," log10Padj\" description=\"BedGraph format\" visibility=full color=0,153,0 priority=20 plotType=\"points\"\n",
file=paste0(folder,"/data/",trackmaster,".log10Padj.bedGraph"))    

     compbdg$padj <- -log10(compbdg$padj)
     compbdg$padj[compbdg$padj == "Inf"] <- max(compbdg$padj[compbdg$padj != "Inf"]) + 1
     write.table(compbdg[,c(1:3,ncol(compbdg))], file=paste0(folder,"/data/",trackmaster,".log10Padj.bedGraph"), row.names=F, col.names=F,sep="\t", quote=F, append=T)
     R.utils::gzip(paste0(folder,"/data/",trackmaster,".log10Padj.bedGraph"), paste0(folder,"/data/",trackmaster,".log10Padj.bedGraph.gz"), overwrite=T)
     cmd <- paste0("create_datauri ",folder,"/data/",trackmaster,".log10Padj.bedGraph.gz")
     datauri <- system(cmd, wait=T, intern=T)
        cat ("
             {
                name: ",paste0("\'",trackmaster,".log10Padj',"),"
                url: ",paste0("\"",datauri,"\""),",
                indexed: false,
                format: \"bedGraph\"
             },\n")

     ## Write dZsc values
cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
track name=\"",trackmaster," dZsc\" description=\"BedGraph format\" visibility=full color=153,76,0 altColor=254,254,254 priority=20 plotType=\"points\"\n",
file=paste0(folder,"/data/",trackmaster,".dZsc.bedGraph"))

     write.table(compbdg[,c(1:3,(ncol(compbdg)-2))], file=paste0(folder,"/data/",trackmaster,".dZsc.bedGraph"), row.names=F, col.names=F,sep="\t", quote=F, append=T)
     R.utils::gzip(paste0(folder,"/data/",trackmaster,".dZsc.bedGraph"), paste0(folder,"/data/",trackmaster,".dZsc.bedGraph.gz"), overwrite=T)
     cmd <- paste0("create_datauri ",folder,"/data/",trackmaster,".dZsc.bedGraph.gz")
     datauri <- system(cmd, wait=T, intern=T) 
        cat ("
             {
                name: ",paste0("\'",trackmaster,".dZsc',"),"
                url: ",paste0("\"",datauri,"\""),",
                indexed: false,
                format: \"bedGraph\"
             }")
      if (i < nrow(file)) { 
       cat (",\n") 
      }
    } else {

        bdg <- read.table(file$file[i],h=F)
        trackname <- file$name[i]
cat("# locus chr19:49302001-49304701\n# refGene encodeRegions\n# zero-based, half-open coords
track name=\"",trackname,"\" description=\"BedGraph format\" visibility=full color=255,102,178 altColor=102,102,255 priority=20 plotType=\"points\"\n",
file=paste0(folder,"/data/",trackname,".bedGraph"))
        write.table(bdg, file=paste0(folder,"/data/",trackname,".bedGraph"),row.names=F, col.names=F,sep="\t", quote=F, append=T)
        R.utils::gzip(paste0(folder,"/data/",trackname,".bedGraph"), paste0(folder,"/data/",trackname,".bedGraph.gz"), overwrite=T)
        cmd <- paste0("create_datauri ",folder,"/data/",trackname,".bedGraph.gz")
        datauri <- system(cmd, wait=T, intern=T)
        cat ("
             {
                name: ",paste0("\'",trackname,"',"),"
                url: ",paste0("\"",datauri,"\""),",
                indexed: false,
                format: \"bedGraph\"
             }")
        if (i < nrow(file)) {
          cat (",\n")
        }
    }
  } 
  scriptbody(genome,end=T) 
}


## Need pip install igv-reports
## Require a visualization.txt file with following coloumns
## file    name    group
## _full_compartment_details.bedGraph A.PC1 compartment [This file is a must and the group should be compartment]
## RNA-Seq_signal.bedGraph	A.RNAseq	RNAseq	[These are optional]
##
## @genome version
## @folder name of the html and data files to be saved

args = commandArgs(trailingOnly=TRUE) 
file <- read.table(args[length(args)], h=T, as.is=T)
genome <- args[length(args)-1]
#file   <- read.table("visualization.txt",h=T,as.is=T)
#genome <- "mm10"
folder <- "vizIGV"
print("Visualization Details")
print(file)

if (!file.exists(folder)) {
  dir.create(folder)
}

sink(paste0(folder,"/igv.html"))
htmlheader()

htmlbody(file, folder, genome)
sink()
