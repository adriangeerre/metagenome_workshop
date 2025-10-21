# Workflow: Metagenome Analysis

#### Software

* R 4.4.3
* python 3.13.9
* jupyter 1.1.1
    * jupyterlab 4.4.7
    * notebook 7.4.5
* r-irkernel 1.3.2
* r-vegan 2.7_2
* bioconductor:
    * ancombc 2.8.0
    * phyloseq 1.50.0
    * microbiome 1.28.0

### Create conda environment

We will start by creating a conda environment with the software required for this tutorial. Instead of conda, we will use mamba to install the packages. There are no differences in the usage of mamba in respect to conda, except that mamba is faster.

*Suggestion:* Install mamba in the base environment of your conda installation. This would make mamba visible in all your environments. Nonetheless, do not install any other package in the base environment.

```
mamba create -n WF5 anaconda::jupyter conda-forge::r-irkernel bioconda::bioconductor-ancombc bioconda::bioconductor-phyloseq conda-forge::r-vegan bioconda::bioconductor-microbiome
```

### Create folders

```
mkdir -p Tutorial/Metagenomics
cd Tutorial/Metagenomics
mkdir Data/
```

### Download data

We will need different files:

* Metadata (Required)
* Isolate counts (Required)
* Taxonomy (Optional)
* Phylogeny (Optional)

```
# Enter folder
cd Data/
mkdir -p tmp
cd tmp

# Download data (3.2 GB)
curl -o wf5_data.tar.gz https://zenodo.org/records/15656403/files/SSC_Community_2025.tar.gz?download=1

# Decompress
tar -xvzf wf5_data.tar.gz

# Extract data for this tutorial
mv Isolate_tables/Original/SSC_norm.tsv ../OTU_table.tsv
mv SSC_taxonomy_GTDB.tsv ../taxa_table.tsv
mv SSC_R2_metadata.tsv ../metadata.tsv

# Delete rest of data
cd ../
rm -r tmp/
```

### Start jupyter notebook within an interactive session

The idea is to run R inside a jupyter notebook to obtain good resources and have the flexibility to analyze the data with visualizations. For this workshop we will refeer to session ([https://genome.au.dk/docs/software-specific/#:~:text=Jupyter%20Notebook/Lab%C2%B6](GenomeDK's interactive jupyter notebook)) session. Remember that you _<jupyter-env>_ is WF5. I recommend to use Visual Studio Code for all this process (folder organization, cluster connection in dual terminal, and notebook analysis).

From now one, each section could be a separate code block.

### Load libraries

```
# Libraries
library(phyloseq)
library(vegan)
library(microbiome)
```

### Prepare and load data into phyloseq

```
# Load data
otus <- read.table("Data/OTU_table.tsv", header=T)
taxa <- read.table("Data/taxa_table.tsv", header=T)
metadata <- read.table("Data/metadata.tsv", header=T)
```

```
# Define row names: OTUs
rownames(otus) <- otus$isolate
otus <- otus[,-1]

# Define row names: Taxa
rownames(taxa) <- taxa$isolate
taxa <- taxa[,-1]

# Define row names: Metadata
rownames(metadata) <- metadata$sample_id
metadata <- metadata[,-1]
```

```
# Make data contain same information
taxa <- taxa[rownames(taxa) %in% rownames(otus),]
metadata <- metadata[rownames(metadata) %in% colnames(otus),]
otus <- otus[,colnames(otus) %in% rownames(metadata)]
```

```
# Transform data into phyloseq class
otus <- otu_table(as.matrix(otus), taxa_are_rows = TRUE)
taxa <- tax_table(as.matrix(taxa))
metadata <- sample_data(metadata)
```

```
# Create phyloseq object
data <- phyloseq(otus, taxa, metadata)
```

### Subset data




### Diversity estimates

Alpha: intra sample richness based on species and abundance
Beta: inter samples richenes based on species (Presence-Absence)

```
# Alpha diversity
alpha(data, index="all")

# KS test

```

```
# Create empty matrix
beta <- matrix(nrow=ncol(data@otu_table), ncol=ncol(data@otu_table))
rownames(beta) <- colnames(data@otu_table)
colnames(beta) <- colnames(data@otu_table)

# Compute pair-wise beta diversity values (takes long time!!)
for (i in colnames(beta)) {
    for (j in colnames(beta)) {
        if (i == j) {
            beta[i,j] <- 1
        } else {
            beta[i,j] <- divergence(abundances(data)[,i], abundances(data)[,j], method="bray")
        }
    }
}

```



### Principal Coordinates analysis (PCoA) / NMDS

### Sample clustering (+ optimal clustering threshold - GAP)

### Run Adonis and Betadisper

Adonis, also known as PERMANOVA (permutation-based multivariate analysis of variance), is a statistical test which allows us to test if two or more samples have a different composition (e.g., bacterial composition). More information, in the [https://rdrr.io/rforge/vegan/man/adonis.html](vegan) package page.


### Normalize data

### Core Microbiome

### Diferential Abundance
