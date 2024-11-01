# CPG DENV  

### Analysing tiled amplicon sequencing data of dengue virus genomes

This repository covers methods for analysing DENV data from both the ONT and Illumina platforms.

The methods described here require working conda environments. They will be created during the Installation and Setup sections below, but first require a working conda installation. One easy method of obtaining this is by following the instructions in the [miniforge](https://github.com/conda-forge/miniforge) repository.  

#### Considerations
I have tried to make this as simple as possible but there are some annoying things to be aware of.  
The easiest way to avoid having to think about these is to follow the installation guides for both the ONT and Illumina analysis tools so that they are both available.  
Cloning the ONT analysis repo (the one which generates the `ont_denv` conda environment) will also give you two utility scripts for the Illumina data:

* [illumina_plotdepth.sh](https://raw.githubusercontent.com/centre-pathogen-genomics/DENV_Amplicon/refs/heads/main/illumina_plotdepth.sh) to plot the depth (of coverage) of each position along the reference  
* [illumina_nextclade.sh](https://raw.githubusercontent.com/centre-pathogen-genomics/DENV_Amplicon/refs/heads/main/illumina_nextclade.sh) to run Nextclade on consensus assemblies for QC and comparison to the reference.

If you are only going to be using Illumina data and don't want to install the repository for ONT data, then you can `wget` or `curl` the scripts directly using the links above.

---
### ONT data  

This uses the [artic](https://github.com/artic-network/fieldbioinformatics) pipeline to assemble genomes, [Nextclade](https://github.com/nextstrain/nextclade) for quality control and comparison to the reference genome, and [csvtk](https://github.com/shenwei356/csvtk) to generate depth (of coverage) plots.   
The BASH script used below to run the pipeline is a modified version of the one written by [Joseph Fauver](https://github.com/josephfauver/DENV_MinION_Script).  

#### Installation and Setup
The commands in the code block below will get everything working by cloning the contents of this repository to your local machine/sever, creating the conda environments, and downloading the Nextclade reference genomes.  
You can install this anywhere, but for the purposes for this documentation we will do it in a new directory called `Tools/` in your home directory (`~`) to make it easy. 
```bash
mkdir ~/Tools/
cd ~/Tools/
git clone https://github.com/centre-pathogen-genomics/DENV_Amplicon.git 
cd DENV_Amplicon
conda create -n ont_denv -c bioconda csvtk medaka=1.11.3 artic nextclade -y
conda activate ont_denv
nextclade dataset get -n community/v-gen-lab/dengue/denv1 -o DENV_Nextclade/DENV1
nextclade dataset get -n community/v-gen-lab/dengue/denv2 -o DENV_Nextclade/DENV2
nextclade dataset get -n community/v-gen-lab/dengue/denv3 -o DENV_Nextclade/DENV3
nextclade dataset get -n community/v-gen-lab/dengue/denv4 -o DENV_Nextclade/DENV4
cd ~
```

#### Expected Input
A single directory (`inputdirectory` in the code block below, but can have any name) with a single fastq.gz file per sample.  

For example:  
```bash
$ ls inputdirectory
sampleA.fastq.gz sampleB.fastq.gz  
sampleC.fastq.gz sampleD.fastq.gz
``` 

#### Running the Pipeline
Make sure the `ont_denv` conda environment is active and use the following command:
```bash
bash ~/Tools/DENV_Amplicon/ont_denv.sh inputdirectory outputdirectory
```

Ensure you run the script with `bash` instead of `sh` - otherwise the step which determines the serotype will not run correctly and will default to DENV1

#### Expected Output
A single directory (`outputdirectory` in the code block above, but can have any name) containing a subdirectory for each sample - the names will be taken from the names of the input fastq.gz files,  
Each sample subdirectory will have loads of different files, following the same naming scheme as the sample subdirectory.  

For example: 
```bash
$ ls outputdirectory
sampleA sampleB sampleC sampleD

$ ls -l outputdirectory/sampleA/
total 26696
-rw-r--r--@ 1 cwwalsh  staff   348649 29 Oct 15:31 sampleA.1.hdf
-rw-r--r--@ 1 cwwalsh  staff    11549 29 Oct 15:31 sampleA.1.vcf
-rw-r--r--@ 1 cwwalsh  staff   371904 29 Oct 15:31 sampleA.2.hdf
-rw-r--r--@ 1 cwwalsh  staff    14600 29 Oct 15:31 sampleA.2.vcf
-rw-r--r--@ 1 cwwalsh  staff   578301 29 Oct 15:30 sampleA.alignreport.er
-rw-r--r--@ 1 cwwalsh  staff   954896 29 Oct 15:30 sampleA.alignreport.txt
-rw-r--r--@ 1 cwwalsh  staff    10941 29 Oct 15:32 sampleA.consensus.fasta
-rw-r--r--@ 1 cwwalsh  staff       86 29 Oct 15:32 sampleA.coverage_mask.txt
-rw-r--r--@ 1 cwwalsh  staff   234820 29 Oct 15:32 sampleA.coverage_mask.txt.1.depths
-rw-r--r--@ 1 cwwalsh  staff   235833 29 Oct 15:32 sampleA.coverage_mask.txt.2.depths
-rw-r--r--@ 1 cwwalsh  staff     9021 29 Oct 15:31 sampleA.depths.pdf
-rw-r--r--@ 1 cwwalsh  staff     9021 29 Oct 15:31 sampleA.depths.txt
-rw-r--r--@ 1 cwwalsh  staff     9021 29 Oct 15:31 sampleA.fail.vcf
-rw-r--r--@ 1 cwwalsh  staff    76971 29 Oct 15:31 sampleA.merged.vcf
-rw-r--r--@ 1 cwwalsh  staff     4371 29 Oct 15:31 sampleA.merged.vcf.gz
-rw-r--r--@ 1 cwwalsh  staff      114 29 Oct 15:31 sampleA.merged.vcf.gz.tbi
-rw-r--r--@ 1 cwwalsh  staff     3561 29 Oct 15:32 sampleA.minion.log.txt
-rw-r--r--@ 1 cwwalsh  staff    21864 29 Oct 15:32 sampleA.muscle.in.fasta
-rw-r--r--@ 1 cwwalsh  staff    21888 29 Oct 15:32 sampleA.muscle.out.fasta
-rw-r--r--@ 1 cwwalsh  staff     7151 29 Oct 15:31 sampleA.pass.vcf.gz
-rw-r--r--@ 1 cwwalsh  staff      114 29 Oct 15:31 sampleA.pass.vcf.gz.tbi
-rw-r--r--@ 1 cwwalsh  staff    10737 29 Oct 15:32 sampleA.preconsensus.fasta
-rw-r--r--@ 1 cwwalsh  staff     1785 29 Oct 15:31 sampleA.primers.vcf
-rw-r--r--@ 1 cwwalsh  staff     1360 29 Oct 15:31 sampleA.primersitereport.txt
-rw-r--r--@ 1 cwwalsh  staff  2716284 29 Oct 15:30 sampleA.primertrimmed.rg.sorted.bam
-rw-r--r--@ 1 cwwalsh  staff       96 29 Oct 15:30 sampleA.primertrimmed.rg.sorted.bam.bai
-rw-r--r--@ 1 cwwalsh  staff  5121033 29 Oct 15:30 sampleA.sorted.bam
-rw-r--r--@ 1 cwwalsh  staff       96 29 Oct 15:30 sampleA.sorted.bam.bai
-rw-r--r--@ 1 cwwalsh  staff  2716512 29 Oct 15:30 sampleA.trimmed.rg.sorted.bam
-rw-r--r--@ 1 cwwalsh  staff       96 29 Oct 15:30 sampleA.trimmed.rg.sorted.bam.bai
```
There will also be a file called `nextclade.tsv` in `outputdirectory` which describes the quality control and serotype assignment results for each sample.  

#### ONT-specific Considerations
The paths to the databases and the model used by medaka to generate the final consensus sequence are hardcoded at the beginning of the script.  
If you need to change these, you can open the `ont_denv.sh` script in your favourite text editor and modify the `SCHEME_DIR` and `MEDAKA_MODEL` varibles respectively.  
The default value of the latter assumes you are using R10.4.1 flow cells, basecalled using the superaccurate model in Dorado.  

There is a potential error that can occur if the model of medaka is too new.   
`medaka 2.x.x` changed the `medaka consensus` command to `medaka inference`, breaking the step which generates the consensus sequence. Make sure you are using `medaka 1.11.3` using the command `medaka --version`. Specifying the medaka version while creating the conda env above should solve this though  
I haven't tested this - but a possible fix would be to change `medaka consensus` to `medaka inference` in the `minion.py` script - located in your conda environment files (eg. `miniforge3/envs/artic/lib/python3.9/site-packages/artic/`).  

---

### Illumina data  

We will use the method described [here](https://github.com/grubaughlab/DENV_pipeline) but with added steps at the end to generate coverage plots and run Nextclade (see the Considerations section near the top for more info).  

#### Installation and Setup
```bash
cd ~/Tools
git clone https://github.com/grubaughlab/DENV_pipeline.git
cd DENV_pipeline
conda env create -f environment.yml -n illumina_denv
conda activate illumina_denv
pip install .
conda install csvtk nextclade -y
nextclade dataset get -n community/v-gen-lab/dengue/denv1 -o DENV_Nextclade/DENV1
nextclade dataset get -n community/v-gen-lab/dengue/denv2 -o DENV_Nextclade/DENV2
nextclade dataset get -n community/v-gen-lab/dengue/denv3 -o DENV_Nextclade/DENV3
nextclade dataset get -n community/v-gen-lab/dengue/denv4 -o DENV_Nextclade/DENV4
cd ~
```
Testing the installation using `denv_pipeline -h` should print the available options.  

#### Optional Step
If you do not have the ONT analysis repository installed, then this is the point at which you should `wget` or `curl` the utility scripts mentioned in the Considerations section near the top of this page.  
Make sure to edit the path to the scripts when using them in the final two steps below.    
If you have the ONT analysis repository installed, then skip this.  
```bash
cd ~/Tools/DENV_pipeline
wget https://raw.githubusercontent.com/centre-pathogen-genomics/DENV_Amplicon/refs/heads/main/illumina_plotdepth.sh
wget https://raw.githubusercontent.com/centre-pathogen-genomics/DENV_Amplicon/refs/heads/main/illumina_nextclade.sh
cd ~
```

#### Expected Input
The expected input is kinda clunky - depending on how your sequencer outputs your FASTQ reads.  There needs to be a single directory (`inputdirectory` in the code block below, but can have any name)  with a subdirectory for each sample.  Each subdirectory should contain forward and reverse read files with `R1` and `R2` somewhere in their filenames. 

If your FASTQ files are all in a single directory, an easy method to get them in the correct structure is:
```bash
ls *R1*fastq.gz | sed 's,_S.*,,' > names
while read i ; do mkdir $i ; done < names
while read i ; do mv $i*R1*fastq.gz $i ; done < names
while read i ; do mv $i*R2*fastq.gz $i ; done < names
rm names
```

Your `inputdirectory` should look something like this:
```bash
$ ls inputdirectory
sampleA sampleB sampleC sampleD

$ ls inputdirectory/*
inputdirectory/sampleA:
sampleA_S1_L001_R1_001.fastq.gz sampleA_S1_L001_R2_001.fastq.gz

inputdirectory/sampleB:
sampleB_S2_L001_R1_001.fastq.gz sampleB_S2_L001_R2_001.fastq.gz

inputdirectory/sampleC:
sampleC_S3_L001_R1_001.fastq.gz sampleC_S3_L001_R2_001.fastq.gz

inputdirectory/sampleD:
sampleD_S4_L001_R1_001.fastq.gz sampleD_S4_L001_R2_001.fastq.gz
```

#### Running the Pipeline
You can run the pipeline with the default settings using the command below. If you want to modify anything, then refer to the source [repository](https://github.com/grubaughlab/DENV_pipeline) for more information.  
```bash
conda activate illumina_denv
denv_pipeline --indir inputdirectory --outdir outputdirectory
```

If you are using modified parameters for every run, then I suggest storing these in a config file and specifying with `--config` for convenience.

#### Expected Output  
A lot of information about the analysis will print to the screen - most of this will also be stored in log files if you need to check anything. The output directory should look like this:
```bash
$ ls outputdirectory
log_files   output_config.yml   results

$ ls outputdirectory/results
alignments                consensus_sequences       low_coverage_calls.csv    summary_all_samples.tsv   variant_plot.pdf          virus_calls.tsv
bam_files                 depth                     low_coverage_consensus    top_virus_all_samples.tsv variants

$ ls outputdirectory/results/consensus_sequences
sampleA.DENV2.10.cons.fa          sampleB.DENV3.10.cons.fa          sampleC.DENV3.10.cons.fa          sampleD.DENV1.10.cons.fa 
```

The consensus sequences will be named `samplename`.`serotype`.`depth`.cons.fa

#### Generating Depth Plots
This can be useful for troubleshooting amplicon dropouts - something that happened with covid as the virus gathered mutations in the primer binding regions.  
Just point the script towards the output directory from the previous step (`outputdirectory` in the example above). It will plot the depth at each position on the genome and store the outputs in `outputdirectory/depth/`.  
```bash
sh ~/Tools/DENV_Amplicon/illumina_plotdepth.sh outputdirectory
```

#### Nextclade
There are lots of useful outputs generated by the pipeline, but if you want to also use Nextclade then run the command below. Just point it at the outputdirectory generated by `denv_pipeline`.
```bash
sh ~/Tools/DENV_Amplicon/illumina_nextclade.sh outputdirectory
```


