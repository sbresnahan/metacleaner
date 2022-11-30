#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
library(taxonomizr)

if(args[[1]]=="NONE"){
  prepareDatabase("accessionTaxa.sql",tmpDir=args[[4]])
  args[[1]] <- paste(getwd(),"accessionTaxa.sql",sep="/")
}

if(args[[6]]%in%c("1","3")){
  badseqids_query <- read.table(args[[2]],sep="\t",header=T)[,1]
  badQ_taxaId <- accessionToTaxa(badseqids_query,args[[1]])
  badQ_taxonomy <- data.frame(getTaxonomy(badQ_taxaId,args[[1]]))
  badQ_taxonomy$taxaID <- row.names(badQ_taxonomy)
  row.names(badQ_taxonomy) <- NULL
  badQ_taxonomy$qseqid <- badseqids_query
  
  badseqids_subject <- read.table(args[[2]],sep="\t",header=T)[,2]
  badS_taxaId <- accessionToTaxa(badseqids_subject,args[[1]])
  badS_taxonomy <- data.frame(getTaxonomy(badS_taxaId,args[[1]]))
  badS_taxonomy$taxaID <- row.names(badS_taxonomy)
  row.names(badS_taxonomy) <- NULL
  badS_taxonomy$sseqid <- badseqids_subject
  
  badseqids <- data.frame(qseqid=badseqids_query,sseqid=badseqids_subject)
}

if(args[[6]]%in%c("1","2")){
  goodseqids_query <- read.table(args[[3]],sep="\t",header=T)[,1]
  goodQ_taxaId <- accessionToTaxa(goodseqids_query,args[[1]])
  goodQ_taxonomy <- data.frame(getTaxonomy(goodQ_taxaId,args[[1]]))
  goodQ_taxonomy$taxaID <- row.names(goodQ_taxonomy)
  row.names(goodQ_taxonomy) <- NULL
  goodQ_taxonomy$qseqid <- goodseqids_query
  
  goodseqids_subject <- read.table(args[[3]],sep="\t",header=T)[,2]
  goodS_taxaId <- accessionToTaxa(goodseqids_subject,args[[1]])
  goodS_taxonomy <- data.frame(getTaxonomy(goodS_taxaId,args[[1]]))
  goodS_taxonomy$taxaID <- row.names(goodS_taxonomy)
  row.names(goodS_taxonomy) <- NULL
  goodS_taxonomy$sseqid <- goodseqids_subject
  
  newbadids <- list()
  for(i in 1:length(goodQ_taxonomy$qseqid)){
    if(goodQ_taxonomy[i,args[[5]]]!=goodS_taxonomy[i,args[[5]]]){
      newbadids <- append(newbadids, goodQ_taxonomy$qseqid[i])
    }
  }
  newbadids <- unlist(newbadids)
  print(paste(">>> ",paste(length(newbadids)," accession(s) filtered",sep="")),sep="")
  goodseqids <- data.frame(qseqid=goodseqids_query,sseqid=goodseqids_subject)
  goodseqids <- goodseqids[!goodseqids$qseqid%in%newbadids,]
  row.names(goodseqids) <- NULL
  if(length(newbadids)>0){
    newbadout <- data.frame(badseqid=unlist(newbadids),
                            tophit=goodseqids[goodseqids$qseqid%in%newbadids,"sseqid"],
                            reason=rep.int(paste("Wrong",args[[5]],sep=" "),
                                          length(newbadids)))
    oldbadids <- read.table(args[[2]],sep="\t",header=T)
    newbadids_out <- rbind(oldbadids,newbadout)
    write.table(newbadids_out,args[[2]],sep="\t",header=T,row.names=F,quote=F)
  }
}

goodout <- function(){
  goodtax_out <- goodQ_taxonomy[,c(9,4,5,6,7)]
  goodtax_out$species <- gsub(" ", "_", goodtax_out$species)
  goodtax_out <- goodtax_out[goodtax_out$qseqid%in%goodseqids$qseqid,]
  write.table(goodtax_out,paste(args[[4]],paste(args[[7]],"_clean.tax",sep=""),sep="/"),
              sep="\t",col.names=F,row.names=F,quote=F)
}

badout <- function(){
  badseqids_out <- read.table(args[[2]],sep="\t",header=T)
  
  badQ_taxonomy$queryTaxa <- paste(badQ_taxonomy$superkingdom,
                               badQ_taxonomy$phylum,
                               badQ_taxonomy$class,
                               badQ_taxonomy$order,
                               badQ_taxonomy$family,
                               badQ_taxonomy$genus,
                               badQ_taxonomy$species,sep="|")
  badQ_taxonomy$queryTaxa <- gsub(" ", "_", badQ_taxonomy$queryTaxa)
  badQ_merge <- badQ_taxonomy[,c("qseqid","queryTaxa")]
  names(badQ_merge) <- c("badseqid","queryTaxa")
  badseqids_out <- merge(x=badseqids_out,y=badQ_merge, 
                         by="badseqid", all.x=TRUE)
  
  badS_taxonomy$subjectTaxa <- paste(badS_taxonomy$superkingdom,
                                     badS_taxonomy$phylum,
                                     badS_taxonomy$class,
                                     badS_taxonomy$order,
                                     badS_taxonomy$family,
                                     badS_taxonomy$genus,
                                     badS_taxonomy$species,sep="|")
  badS_taxonomy$subjectTaxa <- gsub(" ", "_", badS_taxonomy$subjectTaxa)
  badS_merge <- badS_taxonomy[,c("sseqid","subjectTaxa")]
  names(badS_merge) <- c("tophit","subjectTaxa")
  badseqids_out <- merge(x=badseqids_out,y=badS_merge, 
                         by="tophit", all.x=TRUE)
  
  badseqids_out <- badseqids_out[,c("badseqid","queryTaxa","tophit","subjectTaxa","reason")]
  write.table(badseqids_out,args[[2]],sep="\t",col.names=T,row.names=F,quote=F)
}

if(args[[6]]=="1"){
  goodout()
  badout()
}

if(args[[6]]=="2"){
  goodout()
}

if(args[[6]]=="3"){
  badout()
}