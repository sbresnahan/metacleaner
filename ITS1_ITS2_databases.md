# Example of contrasting sequence databases (plant vs non-plant ITS1/ITS2 regions)

## "Bad" sequences (for removal of non-plant sequences from query database):

Downloaded sequences manually from NCBI:
- ((((ITS1) OR 5.8S) OR 28S) OR ITS2) NOT Embryophyta[Organism] AND ("0"[SLEN] : "10000"[SLEN]) 

Removed spaces and cleaned up sequence headers:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_ITS2_nonplant_NCBI.fasta.txt | sed '/^$/d' > ITS1_ITS2_nonplant_NCBI.fasta  
cat ITS1_ITS2_nonplant_NCBI.fasta  | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_ITS2_nonplant.fasta   
rm ITS1_ITS2_nonplant_NCBI.fasta 
```

Combined and removed duplicates using seqkit:
```
conda install -c biobuilds fastx-toolkit  
seqkit rmdup -s < ITS1_ITS2_nonplant.fasta > ITS1_ITS2_nonplant_collapsed.fasta  
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_ITS2_nonplant_collapsed.fasta | sed '/^$/d' > ITS1_ITS2_nonplant_collapsed_fixed.fasta  
cat ITS1_ITS2_nonplant_collapsed_fixed.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_ITS2_nonplant_database.fasta
rm ITS1_ITS2_nonplant_collapsed.fasta 
rm ITS1_ITS2_nonplant_collapsed_fixed.fasta
```

## "Good" sequences (for removal of "mislabeled" sequences from query database):

Downloaded sequences manually from NCBI:
- ((((ITS1) OR 5.8S) OR 28S) OR ITS2) AND Embryophyta[Organism] AND ("0"[SLEN] : "10000"[SLEN]) 

Removed spaces and cleaned up sequence headers:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_ITS2_plant_NCBI.fasta.txt | sed '/^$/d' > ITS1_ITS2_plant_NCBI.fasta   
cat ITS1_ITS2_plant_NCBI.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_ITS2_plant.fasta   
rm ITS1_ITS2_plant_NCBI.fasta     
```

Combined and removed duplicates using fastx-toolkit:
```  
seqkit rmdup -s < ITS1_ITS2_plant.fasta > ITS1_ITS2_plant_collapsed.fasta   
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_ITS2_plant_collapsed.fasta | sed '/^$/d' > ITS1_ITS2_plant_collapsed_fixed.fasta   
cat ITS1_ITS2_plant_collapsed_fixed.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_ITS2_plant_database.fasta
rm ITS1_ITS2_plant_collapsed.fasta 
rm ITS1_ITS2_plant_collapsed_fixed.fasta
```
