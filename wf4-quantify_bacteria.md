# Workflow: Quantify Bacteria

#### Software

* bowtie2 2.5.4
* samtools 1.22.1
* ncurses 6.5
* salmon 1.10.3
* ncbi-dataset 18.9.0
* sra-tools 3.2.1
* python 3.13.7
* biopython 1.85

### Create conda environment

We will start by creating a conda environment with the software required for this tutorial. Instead of conda, we will use mamba to install the packages. There are no differences in the usage of mamba in respect to conda, except that mamba is faster.

*Suggestion:* Install mamba in the base environment of your conda installation. This would make mamba visible in all your environments. Nonetheless, do not install any other package in the base environment.

```
mamba create -n WF4 bioconda::bowtie2 bioconda::samtools anaconda::ncurses bioconda::salmon conda-forge::ncbi-datasets-cli bioconda::sra-tools anaconda::biopython
```

### Create folders

```
mkdir -p Tutorial/QuanBac
cd Tutorial/QuanBac
mkdir -p Data/FNA Data/Fastq
```

### Download data

The idea for this specific workshop is to analyze metagenomics data obtained after the inoculation of a bacterial culture collection or SynCom into a host. Therefore, prior to this we need to have the genome assemblies of the bacteria and the host. Moreover, we need to obtain metagenomics data obtained from a host (e.g., human, plant, etc).In our case, we will be working with three bacteria extracted from the roots of _Lotus japonicus_, and one metagenomics sample.

*Important*: We will not perform the quality control of the metagenomics data. However, this is a crucial step required for the correct analysis of high-throughput biological data. Remember that you have tools like *FastQC* and *MultiQC* for this purpose.

```
# Enter folder
cd Data/FNA

# Download and extract host genome from NCBI
datasets download genome accession GCF_012489685.1 --include genome
unzip ncbi_dataset.zip
cp ncbi_dataset/data/GCF_012489685.1/GCF_012489685.1_*.fna  LjGifu_v1.2.fna
rm -r ncbi_dataset ncbi_dataset.zip md5sum.txt README.md

# Download and extract high-quality genomes from NCBI
datasets download genome accession GCA_964253635.1 --include genome # LjRoot135
unzip ncbi_dataset.zip
cp ncbi_dataset/data/GCA_964253635.1/GCA_964253635.1_*.fna LjRoot135.fna
rm -r ncbi_dataset ncbi_dataset.zip md5sum.txt README.md

datasets download genome accession GCF_964252735.1 --include genome # LjNodule214
unzip ncbi_dataset.zip
cp ncbi_dataset/data/GCF_964252735.1/GCF_964252735.1_*.fna LjNodule214.fna
rm -r ncbi_dataset ncbi_dataset.zip md5sum.txt README.md

datasets download genome accession GCA_964253145.1 --include genome # LjRoot206
unzip ncbi_dataset.zip
cp ncbi_dataset/data/GCA_964253145.1/GCA_964253145.1_*.fna LjRoot206.fna
rm -r ncbi_dataset ncbi_dataset.zip md5sum.txt README.md

# Enter folder
cd ../Fastq

# Download and extract meteganomics data
prefetch SRX25212494
fastq-dump --split-files SRR29709995.sra
rm -r SRR29709995/
ls | grep *.fastq | parallel -j 2 "gzip {}"

# Exit data folder
cd ../../
```

If you want to extend the number of metagenomics samples you can explore the SRA repository belonging to the _PRJNA1131994_ BioProject.

### Start interactive session

The code-box below activates an interactive session in GenomeDK. This session would allow us to work inside the cluster without the need of using a script. Be aware that the <PROJECT> tag needs to be modified for your working project. Also, the number of threads (-n) and the amount of memory (--mem) could be ajusted.

```
# Start a session on GenomeDK
srun --account <PROJECT> -n 8 --mem 16GB --time=8:00:00 --pty /bin/bash

# Activate Environment
conda activate WF4
```

### Remove host reads

```
# Create host index (>30 minutes)
mkdir -p Index
bowtie2-build Data/FNA/LjGifu_v1.2.fna IndexBW2/LjGifu --threads 8

# Map reads against host
mkdir -p HostMap
bowtie2 -x IndexBW2/LjGifu -1 Data/Fastq/SRR29709996_1.fastq.gz -2 Data/Fastq/SRR29709996_2.fastq.gz -p 8 -S HostMap/aln.sam &> bowtie2.log
```

Before continuing the analysis, it is advisable to read the log from the bowtie2 map. There, we will obtain an overview of how many reads where host or potential bacteria out of the total.

_Note_: Do not mind if samtools says "libtinfow.so.6" and "libncursesw.so.6" is missing (.symver is missing). Just verify that the files are produced. You can verify that the libraries are present using `ldd $(which samtools)`.

```
# Extract bacteria reads
samtools view -b -f 4 HostMap/aln.sam > HostMap/unmapped.bam # subset unmapped reads
samtools collate HostMap/unmapped.bam HostMap/unmapped.collate
samtools fastq -1 HostMap/clean_reads_R1.fastq -2 HostMap/clean_reads_R2.fastq -s HostMap/clean_reads_leftover.fastq HostMap/unmapped.collate.bam

# Gzip fastqs
ls -d HostMap/* | grep fastq | parallel -j 3 "gzip {}"
```

To free space, we can remove the intermediate files and keep the resulting fastqs.

### Quantify bacteria

Once we have the metagenomics reads clean from host, we can quantify the _relative_ abudance of bacteria. First, we need to create a single fasta file with the one-line bacterial genomes. In other words, if a genome has more than a chromosome we need to paste it together using N's to create spacers. To simplify this task, here is an script in BioPython:

```
vim paste_genome.py

    ----- Include in script the content below -----

    # Libraries
    import sys
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    # Arguments
    infasta = sys.argv[1]
    outfasta = sys.argv[2]
    header = sys.argv[3]

    # Load genome
    contigs = list(SeqIO.parse(infasta, "fasta"))

    # Join sequences with 100 Ns in between
    concat = Seq("")
    for i, record in enumerate(contigs):
        if i > 0:
            concat += Seq("N" * 100)
        concat += record.seq

    # Create a new SeqRecord
    genome = SeqRecord(
        concat,
        id=header,
        description=""
    )

    # Write to output FASTA with line wrapping at 80 characters
    SeqIO.write(genome, outfasta, "fasta")

    ----- End of script -----

```

Let's paste LjNodule214

```
# Paste genome
python paste_genome.py Data/FNA/LjNodule214.fna Data/FNA/LjNodule214.concat.fna LjNodule214

# Create single fasta
cat Data/FNA/LjRoot135.fna Data/FNA/LjNodule214.concat.fna Data/FNA/LjRoot206.fna > Data/bacterial_genomes.fna

# Verify fasta content (Optional!)
grep "^>" Data/bacterial_genomes.fna
```

Now we can proceed with the bacteria quantification:

```
# Create folders
mkdir -p IndexSalmon BactMap

# Create bacteria index
salmon index -t Data/bacterial_genomes.fna -i IndexSalmon/bact_index -k 31 -p 8

# Map reads against host
salmon quant -i IndexSalmon/bact_index -l IU -1 HostMap/clean_reads_R1.fastq.gz -2 HostMap/clean_reads_R2.fastq.gz --validateMappings -o BactMap --minScoreFraction 0.95 --writeMappings=BactMap/aln.sam --writeUnmappedNames --threads 8
```

Once finished, our results are available in the files `quant.sf`. The column _NumReads_ represent the number of reads mapped to each bacteria (read more about salmon in their github page: [https://combine-lab.github.io/salmon/](https://combine-lab.github.io/salmon/)). Following this, the data can be taken into R or Python for visualization.
