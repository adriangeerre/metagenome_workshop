# Workflow: Metagenome Analysis

#### Software

* R 4.4.1
* jupyter 1.1.1
    * jupyterlab 4.4.7
    * notebook 7.4.5
* r-irkernel 1.3.2
* ancombc 2.8.0
* phyloseq 1.50.0
* r-vegan 2.7_2
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
```

### Download data

We will need different files:

* Metadata (Required)
* Isolate counts (Required)
* Taxonomy (Optional)
* Phylogeny (Optional)

To avoid an long workshop, I am going to use the example data from ([https://joey711.github.io/phyloseq/Example-Data.html](Phyloseq)).

### Start jupyter notebook within an interactive session

The idea is to run R inside a jupyter notebook to obtain good resources and have the flexibility to analyze the data with visualizations. For this workshop we will refeer to session ([https://genome.au.dk/docs/software-specific/#:~:text=Jupyter%20Notebook/Lab%C2%B6](GenomeDK's interactive jupyter notebook)) session. Remember that you _<jupyter-env>_ is WF5. I recommend to use Visual Studio Code for all this process (folder organization, cluster connection in dual terminal, and notebook analysis).

From now one, each section could be a separate code block.

### Load libraries

```
# Libraries
library(phyloseq)
library(vegan)
```

### Load Example data

```
# Example data
data(enterotype)
enterotype
```


### Load our own data (clean)

```
# Create phyloseq object
data = phyloseq(otu_table, tax_table, sample_data, tree)
```

### Identify outliers

### Subset data

### Diversity estimates

### Principal Coordinates analysis (PCoA) / NMDS

### Sample clustering (+ optimal clustering threshold - GAP)

### Run Adonis, Betadisper

### Permanova

### Normalize data

### Core Microbiome

### Diferential Abundance
