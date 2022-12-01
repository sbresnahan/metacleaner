#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
library(taxonomizr)

if(args[[1]]=="NONE"){
  prepareDatabase("accessionTaxa.sql",tmpDir=args[[4]])
  args[[1]] <- paste(getwd(),"accessionTaxa.sql",sep="/")
}


# Check through badseqids, recover qseqids where sseqid is same filterlevel
badseqids_df <- read.table(args[[2]],sep="\t",header=T)

badseqids_query <- badseqids_df$badseqid
badQ_taxaId <- accessionToTaxa(badseqids_query,args[[1]])
badQ_taxonomy <- data.frame(getTaxonomy(badQ_taxaId,args[[1]]))
badQ_taxonomy$taxaID <- row.names(badQ_taxonomy)
row.names(badQ_taxonomy) <- NULL
badQ_taxonomy$qseqid <- badseqids_query

badseqids_subject <- badseqids_df$tophit
badS_taxaId <- accessionToTaxa(badseqids_subject,args[[1]])
badS_taxonomy <- data.frame(getTaxonomy(badS_taxaId,args[[1]]))
badS_taxonomy$taxaID <- row.names(badS_taxonomy)
row.names(badS_taxonomy) <- NULL
badS_taxonomy$sseqid <- badseqids_subject

badseqids_confirmed <- badseqids_df[grepl(args[[7]], badseqids_df$reason),"badseqid"]
badseqids_check <- badseqids_df[grepl(args[[8]], badseqids_df$reason),"badseqid"]
badQ_taxonomy_check <- badQ_taxonomy[badQ_taxonomy$qseqid%in%badseqids_check,]

newgoodids <- list()
for(i in 1:length(badQ_taxonomy_check$qseqid)){
  c <- T
  x <- badQ_taxonomy[i,args[[5]]]
  if(is.na(x)){c <- F}
  y <- badS_taxonomy[i,args[[5]]]
  if(is.na(y)){c <- F}
  if((x==y)&(c==T)){
    newgoodids <- append(newgoodids, badQ_taxonomy_check$qseqid[i])
  }
}
newgoodids <- unlist(newgoodids)
print(paste(">>> ",paste(length(newgoodids),
                         paste(" flagged accession(s) were correct ",args[[5]],sep=""),
                         sep="")),sep="")

goodseqids_df <- read.table(args[[3]],sep="\t",header=T)
newgoodids_df <- badseqids_df[badseqids_df$badseqid%in%newgoodids,]
newgoodids_df <- newgoodids_df[!newgoodids_df$badseqid%in%badseqids_confirmed,]
newgoodids_df$reason <- NULL
names(newgoodids_df) <- names(goodseqids_df)
badseqids_df <- badseqids_df[!badseqids_df$badseqid%in%newgoodids_df$goodseqid,]
goodseqids_df <- rbind(goodseqids_df,newgoodids_df)


# Check through goodseqids, remove qseqids where sseqid is not same filterlevel
goodseqids_query <- goodseqids_df$goodseqid
goodQ_taxaId <- accessionToTaxa(goodseqids_query,args[[1]])
goodQ_taxonomy <- data.frame(getTaxonomy(goodQ_taxaId,args[[1]]))
goodQ_taxonomy$taxaID <- row.names(goodQ_taxonomy)
row.names(goodQ_taxonomy) <- NULL
goodQ_taxonomy$qseqid <- goodseqids_query

goodseqids_subject <- goodseqids_df$tophit
goodS_taxaId <- accessionToTaxa(goodseqids_subject,args[[1]])
goodS_taxonomy <- data.frame(getTaxonomy(goodS_taxaId,args[[1]]))
goodS_taxonomy$taxaID <- row.names(goodS_taxonomy)
row.names(goodS_taxonomy) <- NULL
goodS_taxonomy$sseqid <- goodseqids_subject

newbadids <- list()
for(i in 1:length(goodQ_taxonomy$qseqid)){
  c <- F
  x <- goodQ_taxonomy[i,args[[5]]]
  if(is.na(x)){c <- T}
  y <- goodS_taxonomy[i,args[[5]]]
  if(is.na(y)){c <- T}
  if((x!=y)||(c==T)){
    newbadids <- append(newbadids, goodQ_taxonomy$qseqid[i])
  }
}
newbadids <- unlist(newbadids)
print(paste(">>> ",paste(length(newbadids)+length(badseqids_df$badseqid),
                         " total accession(s) filtered",sep="")),sep="")

if(length(newbadids)>0){
  newbadout <- goodseqids_df[goodseqids_df$goodseqid%in%newbadids,]
  newbadout$reason <- rep.int(paste("Wrong",args[[5]],sep=" "),length(newbadids))
  names(newbadout) <- names(badseqids_df)
  badseqids_df <- rbind(badseqids_df,newbadout)
  write.table(badseqids_df,args[[2]],sep="\t",row.names=F,quote=F)
}

goodseqids_df <- goodseqids_df[!goodseqids_df$goodseqid%in%badseqids_df$badseqid,]
row.names(goodseqids_df) <- NULL


# Save files
goodtax_out <- goodQ_taxonomy[,c(9,4,5,6,7)]
goodtax_out$species <- gsub(" ", "_", goodtax_out$species)
goodtax_out <- goodtax_out[goodtax_out$qseqid%in%goodseqids_df$goodseqid,]
write.table(goodtax_out,paste(args[[4]],paste(args[[6]],"_clean.tax",sep=""),sep="/"),
            sep="\t",col.names=F,row.names=F,quote=F)

badQ_taxonomy <- rbind(badQ_taxonomy,goodQ_taxonomy[goodQ_taxonomy$qseqid%in%badseqids_df$badseqid,])
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
badseqids_df <- merge(x=badseqids_df,y=badQ_merge, 
                       by="badseqid", all.x=TRUE)

badS_taxonomy <- rbind(badS_taxonomy,goodS_taxonomy[goodS_taxonomy$sseqid%in%badseqids_df$tophit,])
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
badseqids_df <- merge(x=badseqids_df,y=badS_merge, 
                       by="tophit", all.x=TRUE)
badseqids_df <- badseqids_df[!duplicated(badseqids_df),]

badseqids_df <- badseqids_df[,c("badseqid","queryTaxa","tophit","subjectTaxa","reason")]
write.table(badseqids_df,args[[2]],sep="\t",col.names=T,row.names=F,quote=F)
