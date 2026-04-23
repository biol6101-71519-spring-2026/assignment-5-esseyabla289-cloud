#DESCRIPTION OF DATASET
The dataset is made up of 4 genes with 10 taxa each. 
#taxa.txt - used a script to get the taxa and genes from NCBI
Homo sapiens - Human
Pan troglodytes - Chimpanzee
Mus musculus - House mouse
Canis lupus familiaris - Domestic Dog
Bos taurus - Cow
Felis catus - Domestic cat
Equus caballus - Horse
Oryctoglagus cuniculus - European rabbit
Sus scrofa - Domestic pig
Rattus norvegicus - Brown rat
#genes
COI (Cytochrome c Oxidase)  - 9 taxa sequences found
CYTB (Cytochrome B) - 9 taxa sequences found
BRCA1 (Breast Cancer gene 1) - 6 taxa sequences found
RAG2 (Recombination Activating Gene 2) - 6 taxa sequences found
#Best-Fit
COI - TIM2+F+G4
CYTB - GTR+F+R2
BRCA1 - HKY
RAG2 - TN+F+I
Due to the differences in nucleotide submissions across the different genes
shows the variation in model complexity for each gene
#SUMMARY OF KEY RESULTS
1. Coalescent species tree - Astral
The coalescent-based species was implied using Astral from four individual gene trees.
The Nodal support on the trees is expressed as local posterior probabilities (LPP).
Therefore deep nodes receieve low support (LPP=0.33) while shallower trees receive 
higher support (LPP = 0.67)
2. Supermatrix tree - IQ-TREE
The supermatrix or concatenated tree was estimated by IQ-TREE2. It used a partition
model and applied it to the 4 genes. Bootstrap support values were set from low-1-10
(for the deep nodes) to high-89-96(for the shallow nodes). The concatenated tree formed
a well-supported early-diverging clades. 
#COMPARISON OF SPECIES COALESCENT VS SUPERMATRIX
Comparing the Coalescent vs Supermatrix methods and these yielded different results,
based on the Robinson-Foulds (RF) distance. The normalised RF distance between 
coalescent and the concatenated trees was 0.7143 which was indicative of a difference
in results between the two methods. Gene trees show a lot of agreement with the 
coalescent tree than the concatenated tree. 
#GENE TREE CONCORDANCE AND DISCORDANCE
Gene concordance factor(gCF) analysis was performed using IQ-TREE which mapped all four
gene unto the ASTRAL species tree. Results reveal discordance across all internal nodes
It showed COI conflicted with the outergroup related nodes while CYTB conflicted with 
internal nodes. BRAC1 and RAG2 conflicted with the deepest internal groups.
The fact that thee genes were only 4 could be a factor in clade evaluation
Every resolved clade had a gCF of 0.00 indicating no single tree supported
bipartition in the species tree indicating that none of the trees contained 
a specific branch from the species tree created. It was observed that the 
nuclear genes (RAG2 and BRCA1) were concordant with the coalescent trees
with RF =0 but the mitochondrial genes (COI and CYTB) had some discordance 
(RF, COI= 0.1429 & 0.8571, CYTB= 0.1429).
#EVALUATION OF APPROPRAITE METHOD
With the results obtained, the coalescent approach using ASTRAL seems better
and appropraite for this dataset. The coalescent framework treats each genes
seperately therefore the discordance seen in concatenated tree will not
be seen with the coalescent. With the coalescent, the gene trees are individually
closer to the coalescent tree. Due to the discordance seen with concatenated
tree, it is better to use the coalescent for this dataset. 

