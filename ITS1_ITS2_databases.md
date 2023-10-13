# Example of contrasting sequence databases (plant vs non-plant ITS1/ITS2 regions)

## "Bad" sequences (for removal of non-plant sequences from query database):

Downloaded list of sequence accessions from NCBI:
- ((((ITS1) OR 5.8S) OR 28S) OR ITS2) NOT Embryophyta[Organism] AND ("0"[SLEN] : "10000"[SLEN]) 

Retrieved fasta sequences using [reutils](https://cran.r-project.org/web/packages/reutils/index.html) in R:
```
# Ensure working directory path is set, and there is a subdirectory called "fasta" here
wd <- setwd("working/directory/path")
dir.create(file.path(wd, "fasta"), showWarnings = FALSE)

# Download fasta sequences in chunks of 5000 and write to individual file for each chunk
library(reutils)
accessions <- read.table("nonplant_accessions.seq.txt",sep="\t",header=F)[,c(1)]
chunks <- split(accessions,cut(seq_along(accessions),5000,labels = FALSE))
for(i in 1:length(chunks)){
  print(i)
  fetch <- efetch(unlist(chunks[i]), db="nuccore", rettype = "fasta", retmode = "text",strand=1)
  write(content(fetch), file = paste(paste("fasta/",i,sep=""),".fasta",sep=""))
  cat("\014")
}
```

Combined chunks to single fasta file:
```
find . -maxdepth 1 -type f -name '*.fasta' -print0 | sort -zV | xargs -0 cat > /Users/sbresnahan/Desktop/NCBI/ITS1_ITS2_nonplant_NCBI_seqs.fasta
```

Removed spaces and cleaned up sequence headers:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_ITS2_nonplant_NCBI_seqs.fasta | sed '/^$/d' > ITS1_ITS2_nonplant_NCBI.fasta  
cat ITS1_ITS2_nonplant_NCBI.fasta  | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_ITS2_nonplant.fasta   
```

Combined and removed duplicates using seqkit:
```
conda install -c biobuilds fastx-toolkit  
seqkit rmdup -s < ITS1_ITS2_nonplant.fasta > ITS1_ITS2_nonplant_collapsed.fasta  
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_ITS2_nonplant_collapsed.fasta | sed '/^$/d' > ITS1_ITS2_nonplant_collapsed_fixed.fasta  
cat ITS1_ITS2_nonplant_collapsed_fixed.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_ITS2_nonplant_database.fasta
```

## "Good" sequences (for removal of "mislabeled" sequences from query database):

Downloaded list of sequence accessions from NCBI:
- ((((ITS1) OR 5.8S) OR 28S) OR ITS2) AND Embryophyta[Organism] AND ("0"[SLEN] : "10000"[SLEN]) 

Retrieved fasta sequences using [reutils](https://cran.r-project.org/web/packages/reutils/index.html) in R:
```
# Ensure working directory path is set, and there is a subdirectory called "fasta" here
wd <- setwd("working/directory/path")
dir.create(file.path(wd, "fasta"), showWarnings = FALSE)

# Download fasta sequences in chunks of 5000 and write to individual file for each chunk
library(reutils)
accessions <- read.table("plant_accessions.seq.txt",sep="\t",header=F)[,c(1)]
chunks <- split(accessions,cut(seq_along(accessions),5000,labels = FALSE))
for(i in 1:length(chunks)){
  print(i)
  fetch <- efetch(unlist(chunks[i]), db="nuccore", rettype = "fasta", retmode = "text",strand=1)
  write(content(fetch), file = paste(paste("fasta/",i,sep=""),".fasta",sep=""))
  cat("\014")
}
```

Combined chunks to single fasta file:
```
find . -maxdepth 1 -type f -name '*.fasta' -print0 | sort -zV | xargs -0 cat > /Users/sbresnahan/Desktop/NCBI/ITS1_ITS2_plant_NCBI_seqs.fasta
```

Removed spaces and cleaned up sequence headers:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_ITS2_plant_NCBI_seqs.fasta | sed '/^$/d' > ITS1_ITS2_plant_NCBI.fasta   
cat ITS1_ITS2_plant_NCBI.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_ITS2_plant.fasta      
```

Combined and removed duplicates using fastx-toolkit:
```  
seqkit rmdup -s < ITS1_ITS2_plant.fasta > ITS1_ITS2_plant_collapsed.fasta   
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_ITS2_plant_collapsed.fasta | sed '/^$/d' > ITS1_ITS2_plant_collapsed_fixed.fasta   
cat ITS1_ITS2_plant_collapsed_fixed.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_ITS2_plant_database.fasta
```
