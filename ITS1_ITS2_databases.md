# Example of contrasting sequence databases (plant vs non-plant ITS1/ITS2 regions)

## Step 1: removal of non-plant sequences from database:

Downloaded sequences manually from NCBI:
- ITS1[All Fields] NOT "Embryophyta"[Organism]
- ITS2[All Fields] NOT "Embryophyta"[Organism]

Removed spaces and cleaned up sequence headers:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1.fasta.txt | sed '/^$/d' > ITS1_nonplant_NCBI.fasta  
cat ITS1_nonplant_NCBI.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_nonplant.fasta  
rm ITS1_nonplant_NCBI.fasta  
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS2.fasta.txt | sed '/^$/d' > ITS2_nonplant_NCBI.fasta  
cat ITS2_nonplant_NCBI.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS2_nonplant.fasta  
rm ITS2_nonplant_NCBI.fasta  
```

Combined and removed duplicates using seqkit:
```
conda install -c biobuilds fastx-toolkit  
cat ITS1_nonplant.fasta ITS2_nonplant.fasta > ITS1_ITS2_nonplant.fasta  
seqkit rmdup -s < ITS1_ITS2_nonplant.fasta > ITS1_ITS2_nonplant_collapsed.fasta  
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_ITS2_nonplant_collapsed.fasta | sed '/^$/d' > ITS1_ITS2_nonplant_collapsed_fixed.fasta  
cat ITS1_ITS2_nonplant_collapsed_fixed.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_ITS2_nonplant_database.fasta  
```

## Step 2: removal of "mislabeled" sequences from database:

Downloaded sequences manually from NCBI:
- ITS1[All Fields] AND "Embryophyta"[Organism]
- ITS2[All Fields] AND "Embryophyta"[Organism]

Removed spaces and cleaned up sequence headers:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_plant_NCBI.fasta.txt | sed '/^$/d' > ITS1_plant_NCBI.fasta   
cat ITS1_plant_NCBI.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_plant.fasta   
rm ITS1_plant_NCBI.fasta   
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS2_plant_NCBI.fasta.txt | sed '/^$/d' > ITS2_plant_NCBI.fasta   
cat ITS2_plant_NCBI.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS2_plant.fasta   
rm ITS2_plant_NCBI.fasta   
```

Combined and removed duplicates using fastx-toolkit:
```
cat ITS1_plant.fasta ITS2_plant.fasta > ITS1_ITS2_plant.fasta   
seqkit rmdup -s < ITS1_ITS2_plant.fasta > ITS1_ITS2_plant_collapsed.fasta   
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ITS1_ITS2_plant_collapsed.fasta | sed '/^$/d' > ITS1_ITS2_plant_collapsed_fixed.fasta   
cat ITS1_ITS2_plant_collapsed_fixed.fasta | perl -pe 's/^>gi\|\d+\|.*\|(.*)\|.*/>$1/' | sed '/^>/ s/ .*//' > ITS1_ITS2_plant_database.fasta   
```
