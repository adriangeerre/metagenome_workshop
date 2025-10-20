# Workflow: Genome Assembly

#### Software

* fastqc 0.12.1
* multiqc 1.27.1
* spades 4.0.0
* quast 5.3.0
* checkm 1.2.3
* python 3.12.8

### Create conda environment

We will start by creating a conda environment with the software required for this tutorial. Instead of conda, we will use mamba to install the packages. There are no differences in the usage of mamba in respect to conda, except that mamba is faster.

*Suggestion:* Install mamba in the base environment of your conda installation. This would make mamba visible in all your environments. Nonetheless, do not install any other package in the base environment.

```
mamba create -n WF1 bioconda::spades bioconda::quast bioconda::checkm-genome bioconda::fastqc bioconda::multiqc
``` 

### Create folder

Let's create a folder to run the tutorial.

```
mkdir -p Tutorial/Assembly
cd Tutorial/Assembly
mkdir -p Data/Reads
```

### Download sequencing data

We need to download the paired-end Illumina sequencing reads from ENA (Bioproject: PRJEB37696), which occupies approx. 1.5 GB.

```
# Enter folder
cd Data/Reads

# Download reads from ENA
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR593/ERR5936217/LjRoot135_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR593/ERR5936217/LjRoot135_R2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR593/ERR5936086/LjNodule214_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR593/ERR5936086/LjNodule214_R2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR593/ERR5936286/LjRoot221_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR593/ERR5936286/LjRoot221_R2.fastq.gz

# Back to the initial folder
cd ../../

```

### Start interactive session

The code-box below activates an interactive session in GenomeDK. This session would allow us to work inside the cluster without the need of using a script. Be aware that the <PROJECT> tag needs to be modified for your working project. Also, the number of threads (-n) and the amount of memory (--mem) could be ajusted.

```
# Start a session on GenomeDK
srun --account <PROJECT> -n 8 --mem 24GB --time=12:00:00 --pty /bin/bash

# Activate Environment
conda activate WF1
```

### Quality control reads

```
# Create output folder
mkdir -p readsQC

# Run fastqc
fastqc -t 8 Data/Reads/*.fastq.gz

# Move results
mv Data/Reads/*.html readsQC/
mv Data/Reads/*.zip readsQC/

# Run multiqc
multiqc readsQC/

# Move results
mv multiqc_* readsQC/
```

Once multiqc has finished, we can download the results into our personal computers (Filezilla, ). By downloading the HTML file, we can check the quality across all the sample reads in one go.

### Assembly genomes

Next, we will assemble the genomes. This process takes several hours, speed up could be achieved by incresing the number of cores and memory.

```
# Create output folder
mkdir -p SPAdes

# Define input directory
indir="Data/Reads"

# Run SPAdes
for name in $(ls ${indir} | grep "R1" | sed 's/_R1.fastq.gz//g'); do
	spades.py -1 ${indir}/${name}_R1.fastq.gz -2 ${indir}/${name}_R2.fastq.gz --careful -t 8 -m 24 -o SPAdes/${name}
done
```

The argument "--careful" is not mandatory, you could play with different arguments to obtain a different assembly.

### Quality control assemblies

Once we have obtain the assemblies, we can check their quality. Quast will provide an overview of the genome (number of contigs, N50, etc) whereas CheckM would provide information about taxonomy, completness and contamination.

```
# Create output folder
mkdir -p assemblyQC

# Quast
for name in $(ls ${indir} | grep "R1" | sed 's/_R1.fastq.gz//g'); do
	quast -o assemblyQC/${name}/quast -t 8 SPAdes/${name}/contigs.fasta
done

# CheckM
for name in $(ls ${indir} | grep "R1" | sed 's/_R1.fastq.gz//g'); do
	checkm lineage_wf -x fasta -t 8 SPAdes/${name} assemblyQC/${name}/checkm --reduced_tree
	checkm qa assemblyQC/${name}/checkm/lineage.ms assemblyQC/${name}/checkm -f assemblyQC/${name}/checkm/results.tsv --tab_table -t 8
done
```

*Note:* Since spades produces multiple fasta file, checkm would contain several rows in the output.