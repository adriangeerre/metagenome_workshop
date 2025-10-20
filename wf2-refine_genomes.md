# Workflow: Refine Culture Collection

#### Software

* python 3.10.14
* biopython 1.84
* pandas 2.2.2
* bowtie2 2.5.4
* samtools 1.21
* bedtools 2.31.1
* qiime2 2024.10.1
* seqkit 2.10.0
* gtdb-tk 2.4.1

### Create conda environment

We will start by creating a conda environment with the software required for this tutorial. Instead of conda, we will use mamba to install the packages. There are no differences in the usage of mamba in respect to conda, except that mamba is faster.

*Suggestion:* Install mamba in the base environment of your conda installation. This would make mamba visible in all your environments. Nonetheless, do not install any other package in the base environment.

In this case we will create a conda environment for Qiime2, and then, we would add some extra packages (if needed).

```
# Qiime2: Amplicon distribution
mamba env create -n WF2 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
```

As extra package, we will install SeqKit, Bedtools, and GTDB-tk. Other packages required, such as Bowtie2, Samtools, Biopython, or Pandas; are installed by Qiime2.

```
# Activate Environment
conda activate WF2

# Install packages
mamba install bioconda::seqkit bioconda::bedtools bioconda::gtdbtk
```

Importantly, GTDB-tk requires to download their database (information displayed at bottom of installation log). The database (_gtdbtk-r226_data.tar.gz_) weigths 131.73G. This step is a bottleneck of the workflow, but it is not required until the last step. Thus, one can download the database in the background while ẃorking on the first steps.

```
# Download GTDB-tk automatically
download-db.sh &
```

If you find issues downloading the GTDB-tk database, please consider downloading the database splitted in multiple parts and follow the steps in the [HOWTO.txt](https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/split_package/HOWTO.txt) file.

### Create folders

Let's create a folder to run the tutorial.

```
mkdir -p Tutorial/RefineCC
cd Tutorial/RefineCC
mkdir src silvadb checkm_results qiime2
```

### Download SILVA Database

Here, we would like to download the database corresponding to the 18S/16S ribosomal RNAs, these correspond to the Small subunit, also known as SSU. The LSU corresponds to the 23/28S ribosomal RNAs. SILVA provides sequences and taxonomy for Qiime2, however, the latest version is not available. Instead, we can download the database directly from qiime ([pre-formatted silva db](https://docs.qiime2.org/2024.10/data-resources/)).

```
# Download qiime2 pre-formatted silva database
wget https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza
```

Remember to verify, checking the MD5 on qiime2 website, and move the database into a folder:

```
# Verify downloaded files
md5sum silva-138-99-*

# and, if good, organize folder!!
mv silva-138-99-* silvadb/
```

Detail about this specific SILVA database release could be obtained from `https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/README.txt`.

### Start interactive session

The code-box below activates an interactive session in GenomeDK. This session would allow us to work inside the cluster without the need of using a script. Be aware that the <PROJECT> tag needs to be modified for your working project. Also, the number of threads (-n) and the amount of memory (--mem) could be ajusted.

```
# Start a session on GenomeDK
srun --account <PROJECT> -n 8 --mem 16GB --time=8:00:00 --pty /bin/bash

# Activate Environment
conda activate WF2
```

### Explore CheckM results

We will first gather all CheckM results, obtained from WF1, into a single table.

```
# Merge results
ls ../Assembly/assemblyQC/ | xargs -I {} sh -c 'cat ../Assembly/assemblyQC/{}/checkm/results.tsv | grep -v "^Bin" | sed "s/^/{}\t/g"' > checkm_results/results_merged.tsv
```

Filter table to show only SPAdes assembly results.

```
# Select contigs assembly
word="contigs"
awk -v word=${word} '($2 == word)' checkm_results/results_merged.tsv > checkm_results/results_merged_clean.tsv
```

Split table on good and bad assemblies.

```
# Define parameters
min_comple=95
max_contam=5

# Filter: Good
awk -v min_comple=${min_comple} -v max_contam=${max_contam} '($14 >= min_comple) && ($15 <= max_contam)' checkm_results/results_merged_clean.tsv > checkm_results/passed_assemblies.tsv

# Filter: Bad
awk -v min_comple=${min_comple} -v max_contam=${max_contam} '($14 < min_comple) && ($15 > max_contam)' checkm_results/results_merged_clean.tsv > checkm_results/failed_assemblies.tsv
```

After dividing the data, we can explore the failed assemblies to determine the potential source of error.

For example:

- Remove small contigs not containing genes (e.g., <1000 bps)
- Visualize the coverage versus the GC content per contig
- Extract the 16S region (and/or subregions) to determine potential contamination

### Explore Coverage versus GC% per contig

First, we will need to map the reads against the contigs. Let's create a script to process the assemblies.

```
vim src/compute_coverage.sh


	----- Include in script the content below -----

	#!/bin/bash

	# Variables
	fasta=${1}
	R1=${2}
	R2=${3}
	prefix=${4}
	indexdir=${5}
	outdir=${6}
	threads=${7}

	# Create folders
	mkdir -p ${indexdir} ${outdir}

	# Create index
	bowtie2-build --threads ${threads} ${fasta} ${indexdir}/${prefix} &> /dev/null

	# Map with Bowtie2
	bowtie2 -x ${indexdir}/${prefix} -1 ${R1} -2 ${R2} -S ${outdir}/${prefix}.sam --threads ${threads} 2> ${outdir}/bowtie2.log

	# SAM to BED
	samtools view -bS ${outdir}/${prefix}.sam -@ ${threads} > ${outdir}/${prefix}.bam
	samtools sort -o ${outdir}/${prefix}.sort.bam -O bam ${outdir}/${prefix}.bam -@ ${threads}
	samtools index -b ${outdir}/${prefix}.sort.bam ${outdir}/${prefix}.sort.bai -@ ${threads}

	# BAM to BED
	bedtools genomecov -d -ibam ${outdir}/${prefix}.sort.bam > ${outdir}/${prefix}.cov

	# Remove intermediate files
	rm ${outdir}/${prefix}.sam ${outdir}/${prefix}.bam

	----- End of script -----


```

Once the script is ready, we can process the assemblies.

```
ls ../Assembly/SPAdes/ | xargs -I {} sh -c 'echo {}; bash src/compute_coverage.sh ../Assembly/SPAdes/{}/contigs.fasta ../Assembly/Data/Reads/{}_R1.fastq.gz ../Assembly/Data/Reads/{}_R2.fastq.gz {} coverage/{}/index coverage/{} 8'
```

Once the process is finished we should obtain a file ended in ".cov" for each of the assemblies. Now, we can proceed to compute the GC content per assembly. Thus, we need a new script to process the assemblies.

```
vim src/compute_gc_content.py


	----- Include in script the content below -----

	# Libraries
	import os
	import sys
	import pandas as pd
	from Bio import SeqIO, SeqUtils

	# Define arguments
	infasta = sys.argv[1]
	outfile = sys.argv[2]

	# Create folder
	dir = os.path.dirname(outfile)
	if not os.path.isdir(dir) and dir != '': os.makedirs(dir)

	# Load fasta
	fasta = SeqIO.to_dict(SeqIO.parse(infasta, "fasta"))

	# Compute length and GC per contig
	values = {i: (len(j.seq), SeqUtils.gc_fraction(j)*100) for i,j in fasta.items()}

	# Values to dataframe and save into table
	df = pd.DataFrame(values.values(), index=values.keys(), columns=['Contig_length','GC_perc'])

	df.to_csv(outfile, sep="\t")

	----- End of script -----


```

```
ls ../Assembly/SPAdes/ | xargs -I {} sh -c 'python src/compute_gc_content.py ../Assembly/SPAdes/{}/contigs.fasta gc_perc/{}.tsv'
```

Finally, we can visualize the GC content versus coverage.

```
vim src/plot_gc_vs_cov.py


	----- Include in script the content below -----

	# Libraries
	import os
	import sys
	import pandas as pd
	import seaborn as sns
	import matplotlib.pyplot as plt

	# Define arguments
	cov = sys.argv[1]
	gc = sys.argv[2]
	outfile = sys.argv[3]

	# Create folder
	dir = os.path.dirname(outfile)
	if not os.path.isdir(dir) and dir != '': os.makedirs(dir)

	# Load data
	cov = pd.read_table(cov, header=None).rename(columns = {0: 'Contig', 1: 'Position', 2: 'Coverage'})
	gc = pd.read_table(gc).rename(columns = {'Unnamed: 0': 'Contig'})

	# Transform data
	mean_cov = cov[['Contig','Coverage']].groupby('Contig').mean()

	# Merge data
	df = pd.merge(mean_cov, gc, on="Contig")

	# Plot
	sns.scatterplot(df, x="GC_perc", y="Coverage", size="Contig_length", edgecolor="black", alpha=0.5)
	plt.savefig(outfile)

	----- End of script -----


```

```
ls ../coverage/ | xargs -I {} sh -c 'python src/plot_gc_vs_cov.py coverage/{}/{}.cov gc_perc/{}.tsv gc_vs_cov/{}.pdf' # Change extension for different format (e.g, png)
```

### Taxonomic identification with Qiime2 using SILVA

We would like to taxonomically identify the 16S Ribosomal RNAs from our genomes using Qiime2. First, we need extract the 16S region/s from our assemblies. For that, we require a target template (i.e., a guiding 16S sequence) for Blastn. We will download E. coli's 16S ribosomal DNA.

```
wget "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NR_024570.1&rettype=fasta" -O Ecoli_16S.fasta
```

Next, let's create a script to extract the 16S region from the assemblies.

```
vim src/extract_16S.sh


	----- Include in script the content below -----

	#!/bin/bash

	# Variables
	assembly=${1}
	target=${2}
	prefix=${3}
	outdir=${4}

	# Create folder
	mkdir -p ${outdir}

	# Blast
	blastn -query ${assembly} -subject ${target} -strand both -outfmt "6 std qseq" > ${outdir}/blast_${prefix}.tsv

	# Save blast hits as fasta
	cut -f 13 ${outdir}/blast_${prefix}.tsv | sed "s/^/>${prefix}\\x0A/g" | sed 's/-//g' > ${outdir}/tmp_${prefix}

	# Number sequences
	counter=1
	while IFS= read -r line; do
		if [[ $line == ">"* ]]; then
			new_header=">${prefix}.${counter}"
			((counter++))
			echo "$new_header"
		else
			echo "$line"
		fi
	done < ${outdir}/tmp_${prefix} > ${outdir}/${prefix}.pre-clean.fasta

	# SeqKit
	seqkit rmdup -s ${outdir}/${prefix}.pre-clean.fasta > ${outdir}/${prefix}.fasta

	# Remove intermediate files
	rm ${outdir}/tmp_${prefix}

	----- End of script -----


```

Once the script is ready, we can process the assemblies.

```
# Loop assemblies
ls ../Assembly/SPAdes/ | xargs -I {} sh -c 'echo {}; bash src/extract_16S.sh ../Assembly/SPAdes/{}/contigs.fasta Ecoli_16S.fasta {} qiime2/{} 8'
```

Finally, we will classify the obtained sequences using Qiime2.

```
vim src/classify_16S.sh


	----- Include in script the content below -----

	#!/bin/bash

	# Variables
	fasta=${1}
	database=${2}
	taxonomy=${3}
	prefix=${4}
	outdir=${5}
	threads=${6}

	# Load fasta into qiime2 format
	qiime tools import --type 'FeatureData[Sequence]' --input-path ${fasta} --output-path ${outdir}/${prefix}.qza

	# Classify 16S
	qiime feature-classifier classify-consensus-vsearch --i-query ${outdir}/${prefix}.qza --i-reference-reads ${database} --i-reference-taxonomy ${taxonomy} --p-threads ${threads} --output-dir ${outdir}/Blast_16S --o-classification ${outdir}/${prefix}_16S.taxonomy.qza
	qiime tools export --input-path ${outdir}/${prefix}_16S.taxonomy.qza --output-path ${outdir}
	mv ${outdir}/taxonomy.tsv ${outdir}/${prefix}_16S.taxonomy.tsv

	----- End of script -----

```

```
ls qiime2/ | xargs -I {} sh -c 'echo {}; bash src/classify_16S.sh qiime2/{}/{}.fasta silvadb/silva-138-99-seqs.qza silvadb/silva-138-99-tax.qza {} qiime2/{}'
```

Once finished, we would need to explore the results in order to identify potential source/s of contamination across our assemblies.

**Important:** double-check in qiime2 if the species taxonomic level could be trusted. Otherwise, limit your analysis to the genus taxonomic level.

### Run GTDB-tk

Once we have obtained a final set of genomes to work with, we can obtain a taxonomic identification using GTDB-tk. However, GTDB-tk requires an input folder containing the set of genomes. Here, I will assume that the folder is named as `culture_collection`.

```
# Run GTDB-tk after organizing the genomes in a folder
gtdbtk classify_wf --genome_dir culture_collection --out_dir gtdbtk -x fasta --prefix cultcoll --cpus 8 --pplacer_cpus 8 --write_single_copy_genes --keep_intermediates
```

**Important**: If you have problems running GTDB-tk, verify that the database and the software version are compatible. We experience issues with versions below 2.4.1 using databases below version 226, and viceversa. Apparently, newer databases do not include the FastANI database (potential solution: --skip_ani_screen).



### Extra: 

This steps might require larger memory than the previously defined in the srun command. I recommend to double the amount of memory.

**Classify 16S subregions**

```
vim src/subset_silvadb.sh


	----- Include in script the content below -----

	#!/bin/bash

	# Variables
	database=${1}
	taxonomy=${2}
	primer1=${3}
	primer2=${4}
	prefix=${5}
	outdir=${6}
	threads=${7}

	# Create folder
	mkdir -p ${outdir}

	# Subset database (classifier)
	qiime feature-classifier extract-reads --i-sequences ${database} --p-f-primer ${primer1} --p-r-primer ${primer2} --p-min-length 300 --p-max-length 500 --o-reads ${outdir}/${prefix}.qza --p-n-jobs ${threads}

	# Train classifier
	qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ${outdir}/${prefix}.qza --i-reference-taxonomy ${taxonomy} --o-classifier ${outdir}/${prefix}.qza

	----- End of script -----


```

```
vim src/classify_16S_subregion.sh


	----- Include in script the content below -----

	#!/bin/bash

	# Variables
	fasta=${1}
	database=${2}
	primer1=${3}
	primer2=${4}
	prefix=${5}
	outdir=${6}
	threads=${7}

	# Create folder
	mkdir -p ${outdir}

	# Load fasta into qiime2 format
	qiime tools import --type 'FeatureData[Sequence]' --input-path ${fasta} --output-path ${outdir}/${prefix}.qza

	# Extract region
	qiime feature-classifier extract-reads --i-sequences ${outdir}/${prefix}.qza --p-f-primer ${primer1} --p-r-primer ${primer2} --p-min-length 250 --p-max-length 800 --o-reads ${outdir}/${prefix}.qza
	qiime tools export --input-path ${outdir}/${prefix}.qza --output-path ${outdir}/${prefix}

	# Move file
	mv ${outdir}/${prefix}/dna-sequences.fasta ${outdir}/${prefix}.pre-clean.fasta

	# Taxonomy
	qiime feature-classifier classify-sklearn --i-classifier ${database} --i-reads ${outdir}/${prefix}.qza --o-classification ${outdir}/${prefix}.taxonomy.qza
	qiime tools export --input-path ${outdir}/${prefix}.taxonomy.qza --output-path ${outdir}/${region}
	mv ${outdir}/${region}/taxonomy.tsv ${outdir}/${prefix}.taxonomy.tsv
	rmdir ${outdir}/${region}

	# Clean fasta
	seqkit rmdup -s ${outdir}/${prefix}.pre-clean.fasta > ${outdir}/${prefix}.fasta

	----- End of script -----


```

```
# V3-V4

primer1="CCTACGGGNGGCWGCAG"
primer2="GACTACHVGGGTATCTAATCC"

bash src/subset_silvadb.sh silvadb/silva-138-99-seqs.qza silvadb/silva-138-99-tax.qza ${primer1} ${primer2} silva-138-99-seqs-v3v4 silvadb 8

ls ../Assembly/SPAdes/ | xargs -I {} sh -c "echo {}; bash src/classify_16S_subregion.sh ../Assembly/SPAdes/{}/contigs.fasta silvadb/silva-138-99-seqs-v3v4.qza $primer1 $primer2 {}_v3v4 qiime2/{} 8"

```

```
# V5-V7

primer1="AACMGGATTAGATACCCKG"
primer2="ACGTCATCCCCACCTTCC"

bash src/subset_silvadb.sh silvadb/silva-138-99-seqs.qza silvadb/silva-138-99-tax.qza ${primer1} ${primer2} silva-138-99-seqs-v5v7 silvadb 8

ls ../Assembly/SPAdes/ | xargs -I {} sh -c "echo {}; bash src/classify_16S_subregion.sh ../Assembly/SPAdes/{}/contigs.fasta silvadb/silva-138-99-seqs-v5v7.qza $primer1 $primer2 {}_v5v7 qiime2/{} 8"

```
