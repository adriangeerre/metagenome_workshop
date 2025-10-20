# Workflow: Genomes to Pangenome

#### Software

* python 3.13.2
* fastani 1.34
* eggnog-mapper 2.1.12
* orthofinder 2.5.5

### Create conda environment

We will start by creating a conda environment with the software required for this tutorial. Instead of conda, we will use mamba to install the packages. There are no differences in the usage of mamba in respect to conda, except that mamba is faster.

*Suggestion:* Install mamba in the base environment of your conda installation. This would make mamba visible in all your environments. Nonetheless, do not install any other package in the base environment.

```
mamba create -n WF3 anaconda::python bioconda::orthofinder bioconda::fastani bioconda::eggnog-mapper
```

### Create folders

Let's create a folder to run the tutorial.

```
mkdir -p Tutorial/Gen2Pang
cd Tutorial/Gen2Pang
mkdir -p Data/FAA Data/FNA
```

### Download genomes

For this task we have two options:

1. We could use the genomes assembled in the first workflow.
2. We could download the genomes already assembled from the At-Sphere website.

```
# Enter folder
cd Data/FNA

# Download genomes from At-Sphere
wget --no-check-certificate http://www.at-sphere.com/download/assemblies/LjRoot135.fna.gz
wget --no-check-certificate http://www.at-sphere.com/download/assemblies/LjNodule214.fna.gz
wget --no-check-certificate http://www.at-sphere.com/download/assemblies/LjRoot221.fna.gz

# Ungzip files
gzip -d \*.gz
```

### Download proteins

For this task we have two options:

1. We could run software to annotate the genomes assembled in the first workflow, and then use the proteins to run Orthofinder and EggNog. Examples: Prokka or PGAP (not yet explained!!!).
2. If you decided to download the genome, you should download the proteins from the At-Sphere website too.

```
# Enter folder
cd ../FAA

# Download protein fastas from At-Sphere
wget --no-check-certificate http://www.at-sphere.com/download/ORFs_AA/LjRoot135.faa.gz
wget --no-check-certificate http://www.at-sphere.com/download/ORFs_AA/LjNodule214.faa.gz
wget --no-check-certificate http://www.at-sphere.com/download/ORFs_AA/LjRoot221.faa.gz

# Ungzip files
gzip -d \*.gz

# Back to the initial folder
cd ../../
```

### Start interactive session

The code-box below activates an interactive session in GenomeDK. This session would allow us to work inside the cluster without the need of using a script. Be aware that the <PROJECT> tag needs to be modified for your working project. Also, the number of threads (-n) and the amount of memory (--mem) could be ajusted.

```
# Start a session on GenomeDK
srun --account <PROJECT> -n 8 --mem 16GB --time=8:00:00 --pty /bin/bash

# Activate Environment
conda activate WF3
```

### Run FastANI

Once we obtained our genomic dataset, we need to explore the similarity between genomes. Including highly similar or even equal bacterial strains into our analysis could create bias and confound our conclusions. FastANI is a software that performed kmer (genomic sequences of size K) comparison between pairs of genomes reporting a similarity value between 70 and 100%. Values below 70% are not reported, thus considered "unrelated" organisms.

```
# Create output folder
mkdir FastANI

# Run pairwise fastani in a loop
for i in $(ls | grep fna); do
	for j in $(ls | grep fna); do
		fastANI -q ${i} -r ${j} -o FastANI/${i}\_${j}.txt
	done
done
```

The results from FastANI should be followed with data exploration and redundancy cleaning, as described in the next section.

### Explore FastANI results

By exploring FastANI results we might be able to identify highly similar isolates, which might bias our downstream analysis. To explore the data and cluster the isolates, if required, we should implement an script for easy data processing.

```
vim src/create_fastani_table.py

    ----- Include in script the content below -----

    # Libraries
    import os
    import sys
    import glob
    import pandas as pd

    # Define arguments
    indir = sys.argv[1]
    outfile = sys.argv[2]

    # Create folder
    dir = os.path.dirname(outfile)
    if not os.path.isdir(dir) and dir != '': os.makedirs(dir)

    # Load files
    files = glob.glob(f"{indir}/*")

    # Function: Create Table
    def create_table(files):
        # Store values
        vals = {}
        # Loop files
        for i in files:
            # Obtain isolate names
            filename = os.path.basename(i).replace(".txt","").replace(".fna","").split("_")
            a = filename[0]
            b = filename[1]
            # Read file
            try:
                val = pd.read_table(i, header=None)[2].to_list()[0]
            except:
                val = "NA"
            # Append
            if a in vals.keys():
                vals[a][b] = val
            else:
                vals[a] = {b: val}
        # Return
        return vals

    # Create table
    table = create_table(files)

    # Table to dataframe
    df = pd.DataFrame(table)

    # Order dataframe
    df = df.loc[:,df.index]

    # Save dataframe
    df.to_csv(outfile, sep="\t")

    ----- End of script -----


```

Let's create the table.

```
python src/create_fastani_table.py FastANI/ fastani_table.tsv
```

I would recommend to use R/Python to load the table and explore the potential clusters. Examples of similarities thresholds could be 95% (species differentiation), between 97 to 100% (as Vsearch), or 99.99% (used in the SSC experiment).

### Run EggNogMapper

We can perform the annotation of the proteins for the cleaned community. EggNogMapper would provide an annotation per protein from multiple databases, e.g., COG, KEGG, or PFAM.

```
# Run in login server
download_eggnog_data.py &
```

```
# Create output folder
mkdir -p EggNog

# Run eggnog for each proteins
for faa in $(ls -d Data/FAA/\*.faa); do
    name=$(echo ${faa} | sed 's/Data\/FAA\///g' | sed 's/.faa//g')
    if [ ! -f ${name}.emapper.annotations ]; then
        emapper.py -i ${faa} -m diamond --itype proteins --output EggNog/${name} --cpu 4
    fi
done
```

### Run Orthofinder

Finally, we can use a pangenome analysis to obtain a whole view of the proteins and the sharedness across isolates. Also, we could explore specific subsets by taxonomy or any other relevant metadata. This analysis took a month and a half of computation for ~1000 isolates.

```
# Run eggnog for all proteins
orthofinder -f Data/FAA -t 4 -a 2 -n tutorial -o Pangenome
```
