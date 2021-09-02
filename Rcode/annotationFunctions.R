require(Homo.sapiens)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(Gviz)
require(gtools)

getGenesLocations <- function (geneSymbolsSEL, sortByChrom=TRUE, csvFileName=NULL){
  if (length(geneSymbolsSEL)>0){
    anotacs<- select(Homo.sapiens, 
                   keys=geneSymbolsSEL, columns="ENTREZID", keytype="SYMBOL")
    transcriptCoords2 <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,
                              keys = anotacs$ENTREZID,
                              columns=c('GENEID', 'TXNAME', 'TXID', 'TXCHROM', 'TXSTART', 'TXEND' ),
                              keytype="GENEID")            
    names(transcriptCoords2) <- c('ENTREZID', 'TXNAME', 'TXID', 'TXCHROM', 'TXSTART', 'TXEND')
    anotacs2<- merge(anotacs, transcriptCoords2, 'ENTREZID')
    anotacs3a <- anotacs2[!duplicated(anotacs2$ENTREZID),]
    anotacs3b <- na.omit(anotacs3a)
    if (sortByChrom){
      anotacs3 <- anotacs3b[order(anotacs3b$TXCHROM), ]
    }else{
     anotacs3 <- anotacs3b
    }
    if(!is.null(csvFileName)) write.csv(anotacs3, file=csvFileName, row.names = FALSE)
  }else{
    anotacs3=NULL
  }
  return(anotacs3)
}

plotGenesInChroms <- function (transcriptCoords, plotsFilename, minbase, maxbase, islandData, dnaseData){
  if(!is.null(transcriptCoords)){
    anotacs4<-transcriptCoords[,c('TXCHROM', 'TXSTART', 'TXEND')]
    anotacs4<-anotacs4[complete.cases(anotacs4),]
  
    #change column names to read it into GenomicRanges
    colnames(anotacs4)<-c("chromosome","start","end")
    genRangList<-makeGRangesFromDataFrame(anotacs4) #if we keep geneid we add keep.extra.columns=TRUE inside function
  
    axisT<-GenomeAxisTrack()
  
    data <- read.table(paste("dades","cytoBandIdeo.txt", sep="/"), header=F, sep="\t")
    colnames(data) <-c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')
    
  
    pdf(file=plotsFilename, width= 8, height = 12)
    #draw chromosomes (it has to be done 1 by 1)
    for (i in 1:length(unique(anotacs4$chromosome))){
      chr<-mixedsort(unique(anotacs4$chromosome))[i]
      #draw axis list from our genomic ranges object containing the gene positions
      #read each gene per chromosome
      genList<-AnnotationTrack(anotacs4, name = "Genes", genome ="hg19", chromosome = chr,  stacking ="dense", col= "#5E2366", fill= "#5E2366")
      ideoT<-IdeogramTrack(chromosome = chr,genome="hg19", bands=data,  showId=FALSE)
      # CpG island track
      islandData2 <- islandData[seqnames(islandData) == chr &  (start(islandData) >= minbase & end(islandData) <= maxbase)]
      islandTrack <- AnnotationTrack(range=islandData2, genome="hg19", name="CpG Islands", 
                                 chromosome=chr)
      # DNaseI hypersensitive site data track
      dnaseTrack <- DataTrack(range=dnaseData, genome="hg19", name="DNAseI", 
                            type="gradient", chromosome=chr)
      plotTracks(list(ideoT,axisT,genList, dnaseTrack, islandTrack), sizes=c(1,2,2,1,20), main=chr, cex.main=1,littleTicks = TRUE, showTitle=TRUE)
    }
    dev.off()
  }
}

