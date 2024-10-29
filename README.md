# CPG DENV  

### Analysing tiled amplicon sequencing data of dengue virus genomes

This repository covers methods for analysing DENV data from both the ONT and Illumina platforms.

The methods described here require working conda environments.  
They will be created below, but first require a working conda installation.  
One easy method of obtaining this is by following the instructions in the [miniforge](https://github.com/conda-forge/miniforge) repository.

---
### ONT data  

#### Installation and Setup
Follow the instructions here to get everything working by cloning the contents of this repository to your local machine or sever and creating the conda environments.  

```bash
git clone https://github.com/centre-pathogen-genomics/DENV_Amplicon.git 
cd DENV_Amplicon
conda create -f ont_denv.yml
```

This uses the artic pipeline with reference genomes and BED files for each DENV serotype.  
The BASH script below is a slightly modified version of the one written by [Joseph Fauver](https://github.com/josephfauver/DENV_MinION_Script).
#### Expected Input
A single directory (`inputdirectory` in the code block below, but can have any name) with a single fastq.gz file per sample.  
For example:  
```bash
$ ls inputdirectory
sampleA.fastq.gz sampleB.fastq.gz  
sampleC.fastq.gz sampleD.fastq.gz
``` 
#### Running the Pipeline
Use the following command
```bash
bash ont_denv.sh inputdirectory outputdirectory
```
make sure you run the script with `bash` instead of `sh` - otherwise the step which determines the serotype will not run correctly and will default to DENV1

#### Expected Output
A single directory (`outputdirectory` in the code block above, but can have any name) containing a subdirectory for each sample - the names will be taken from the names of the input fastq.gz files,  
Each sample subdirectory will have loads of different files, following the same naming scheme as the sample subdirectory.  
For example: 
```bash
$ ls outputdirectory
sampleA sampleB sampleC sampleD

$ ls outputdirectory/sampleA/
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

The paths to the database (referred to artic as the scheme), the expected min and max read lengths, and the model used by medaka to genrate the final consensus sequence should be set correctly by default.  
If you need to change these, you can open the `ont_denv.sh` script in your favourite text editor and modify the `SCHEME_DIR`, `MIN_READLEN`, `MAX_READLEN`, and `MEDAKA_MODEL` varibles respectively.  
The latter assumes you are using R10.4.1 flow cells, basecalled using the superaccurate model in Dorado.  

There is a potential error that can occur if the model of medaka is too new. V2 changed the `medaka consensus` command to `medaka inference`, breaking the step which generates the consensus sequence. Make sure you are using `medaka 1.11.3` using the command `medaka --version`.  
I haven't tested this - but a possible fix would be to change `medaka consensus` to `medaka inference` in the `minion.py` script - located in your conda environment files (eg. `miniforge3/envs/artic/lib/python3.9/site-packages/artic/`).  

---
### Illumina data  
We use the method described [here](https://github.com/grubaughlab/DENV_pipeline) but with an added step at the end to generate some plots.

#### Installation and Setup
```bash
git clone https://github.com/grubaughlab/DENV_pipeline.git
cd DENV_pipeline
conda env create -f environment.yml -n illumina_denv
conda activate illumina_denv
pip install .
conda install csvtk -y
```
Testing the installation using `denv_pipeline -h` should print the available options.  

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

Your inputdirectory should look something like this:
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
Just point the script towards the output directory from the previous step (`outputdirectory` in the example above). Tt will plot the depth at each position on the genome and store the outputs in `outputdirectory/depth/`.  
```bash
sh denv_illumina_plotdepth.sh outputdirectory
```



