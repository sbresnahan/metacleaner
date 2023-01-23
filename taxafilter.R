#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
library(taxonomizr)

if(args[[1]]=="NONE"){
  prepareDatabase("accessionTaxa.sql",tmpDir=args[[4]])
  args[[1]] <- paste(getwd(),"accessionTaxa.sql",sep="/")
}


# Check through badseqids
badseqids_all <- read.table(args[[2]],sep="\t",header=T)
badseqids_df <- badseqids_all[grep(args[[8]],badseqids_all$reason),]
badseqids_confirmed <- badseqids_all[grepl(args[[7]], badseqids_all$reason),"badseqid"]

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

## recover qseqids where sseqid is same filterlevel
newgoodids <- list()
for(i in 1:length(badQ_taxonomy$qseqid)){
  c <- T
  x <- badQ_taxonomy[i,args[[5]]]
  if(is.na(x)){c <- F}
  y <- badS_taxonomy[i,args[[5]]]
  if(is.na(y)){c <- F}
  if((x==y)&(c==T)){
    newgoodids <- append(newgoodids, badQ_taxonomy$qseqid[i])
  }
}
newgoodids <- unlist(newgoodids)
print(paste0(">>> ",paste0(length(newgoodids),
                         paste0(" flagged accession(s) were correct ",args[[5]]))))
goodseqids_df <- read.table(args[[3]],sep="\t",header=T)
newgoodids_df <- badseqids_df[badseqids_df$badseqid%in%newgoodids,]
newgoodids_df$reason <- NULL
names(newgoodids_df) <- names(goodseqids_df)
badseqids_df <- badseqids_all[!badseqids_all$badseqid%in%newgoodids_df$goodseqid,]
goodseqids_df <- rbind(goodseqids_df,newgoodids_df)
row.names(goodseqids_df) <- NULL


# Check through goodseqids
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

## Remove qseqids with any hits against arg_addfilter
if(!args[[9]]=="NONE"){
  newbadids <- list()
  aflevel <- strsplit(args[[9]],",")[[1]][1]
  af <- strsplit(args[[9]],",")[[1]][2]
  q <- unique(goodQ_taxonomy$qseqid)
  
  for(i in 1:length(q)){
    x <- goodS_taxonomy[goodQ_taxonomy$qseqid==q[i],aflevel]
    if(any(!x%in%af)){
      newbadids <- append(newbadids, q[i])
    }
  }
  
  newbadids <- unlist(newbadids)
  print(paste0(">>> ",paste0(length(newbadids),
                           paste0(" accession(s) filtered for wrong ",args[[9]]))))

 if(length(newbadids)>0){
   newbadout <- goodseqids_df[goodseqids_df$goodseqid%in%newbadids,]
   newbadout$reason <- rep.int(paste("Wrong",aflevel,sep=" "),length(newbadout$goodseqid))
   names(newbadout) <- names(badseqids_df)
   badseqids_df <- rbind(badseqids_df,newbadout)
   write.table(badseqids_df,args[[2]],sep="\t",row.names=F,quote=F)
   goodseqids_df <- goodseqids_df[!goodseqids_df$goodseqid%in%badseqids_df$badseqid,]
   row.names(goodseqids_df) <- NULL
 }
}

## remove qseqids where sseqid is not same filterlevel
newbadids <- list()
q <- unique(goodQ_taxonomy$qseqid)
reason <- list()
for(i in 1:length(q)){
  c <- F
  x <- goodQ_taxonomy[goodQ_taxonomy$qseqid==q[i],args[[5]]]
  if(any(is.na(x))){
    c <- T
    reason <- append(reason,rep.int("No taxa info for query",length(x)))
  }
  if(c==F){
    y <- unique(goodS_taxonomy[goodQ_taxonomy$qseqid==q[i],args[[5]]])
    y[is.na(y)] <- "FILTER"
    if(length(y)>1){
      c <- T
      reason <- append(reason,rep.int("Matches to multiple genera",length(x)))
    }else{
      if(!y%in%x){
        c <- T
        reason <- append(reason,rep.int("Match to wrong genera",length(x)))
      }
    }
  }
  if(c==T){newbadids <- append(newbadids, q[i])}
}
newbadids <- unlist(newbadids)
print(paste0(">>> ",paste0(length(newbadids),
                           paste0(" accession(s) filtered for wrong ",args[[5]]))))
if(length(newbadids)>0){
  newbadout <- goodseqids_df[goodseqids_df$goodseqid%in%newbadids,]
  newbadout$reason <- unlist(reason)
  names(newbadout) <- names(badseqids_df)
  badseqids_df <- rbind(badseqids_df,newbadout)
  write.table(badseqids_df,args[[2]],sep="\t",row.names=F,quote=F)
}
goodseqids_df <- goodseqids_df[!goodseqids_df$goodseqid%in%badseqids_df$badseqid,]
row.names(goodseqids_df) <- NULL


# Save files
## Cleaned tax
goodtax_out <- unique(goodseqids_df$goodseqid)
taxaId <- accessionToTaxa(goodtax_out,args[[1]])
taxonomy <- getTaxonomy(taxaId,args[[1]])
Final.df <- cbind(taxaId ,goodtax_out, taxaId, taxonomy)
write.table(Final.df,paste(args[[4]],paste0(args[[6]],"_clean.tax"),sep="/"),
            row.names=F,col.names=F,sep=",")

## badseqids
badseqids_df <- badseqids_df[!duplicated(badseqids_df),]
badseqids_query <- badseqids_df$badseqid
badQ_taxaId <- accessionToTaxa(badseqids_query,args[[1]])
badQ_taxonomy <- data.frame(getTaxonomy(badQ_taxaId,args[[1]]))
badQ_taxonomy$taxaID <- row.names(badQ_taxonomy)
row.names(badQ_taxonomy) <- NULL
badQ_taxonomy$qseqid <- badseqids_query
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

badseqids_subject <- badseqids_df$tophit
badS_taxaId <- accessionToTaxa(badseqids_subject,args[[1]])
badS_taxonomy <- data.frame(getTaxonomy(badS_taxaId,args[[1]]))
badS_taxonomy$taxaID <- row.names(badS_taxonomy)
row.names(badS_taxonomy) <- NULL
badS_taxonomy$sseqid <- badseqids_subject
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
badseqids_df <- badseqids_df[order(badseqids_df$badseqid),]
write.table(badseqids_df,args[[2]],sep="\t",col.names=T,row.names=F,quote=F)
