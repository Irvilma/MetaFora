# MetaFora
### Introduction
The **MetaFora** (Metagenomics For All) pipeline analyzes metagenomic data based on 16S DNA obtained from the application of next-generation sequencing (NGS) platforms. The files should be demultiplexed.

### Pipeline workflow
**Stage 1. Pre-process the raw data**  
- Parse the raw reads  
- Trim based on quality and adaptors  
- Quality control

**Stage 2. Quantification of metagenomic data**   
- Submodule I: Taxonomy profile
- Submodule II: Diversity profile
- Submodule III: Comparative integrative analysis


A schematic of the workflow: 

![](https://github.com/Irvilma/MetaFora/blob/master/Figures/workflow.png)

### Requirements

* [FastQC][5] v.0.11.7 (Andrews, 2010)
* [Java 7 or higher][6] - We used Java 8
* [Trimmomatic][7] v.0.36 (Bolger et al., 2014)
* [PEAR][8] v. 0.9.11 (Zhang et al., 2014)
* [QIIME][9] v.1.9.1 (Caporaso, et al., 2010)
* [Vsearch][10] v2.9.1 (Rognes et al., 2016)
* [The R Project for Statistical Computing][13](Team, 2000)
* [dyplr][14] R package (Wickham et al., 2018)
* [purrr][15] R package(Henry and Wickham, 2018)


### Briefly Manual

**Installing QIIME**

QIIME is a set of python scripts that are called using the terminal.
Authors recommend installing [QIIME][11] through Miniconda package. Here, you can see a guide for Linux OS.

1) Install [Miniconda][12] package. 
2) Create your qiime1 environmental and install QIIME.

$ conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda

*Install the additional components*

Add additional channels for required packages

$ conda config --add channels defaults

$ conda config --add channels r

$ conda config --add channels bioconda

$ conda config --add channels biocore

$ conda config --add channels hcc 
 
Install additional packages

$ conda install --yes r-essentials blast-legacy cd-hit cdbtools cogent ea-utils infernal microbiomeutil muscle mothur r-vegan sourcetracker sortmerna sumaclust rdp-classifier swarm tax2tree

3) Activate your qiime1 environment and test your QIIME installation.

$ source activate qiime1

$ print_qiime_config.py –t


**NOTE:** If you have any problems with the proper functioning of QIIME, you should check the versions of these packages:

NumPy version:	1.15.0

pandas version:	0.23.1

matplotlib version:	1.4.3




**Installing vsearch v2.9.1**

- Download and install vsearch:

$ wget https://github.com/torognes/vsearch/archive/v2.9.1.tar.gz

$ tar xzf v2.9.1.tar.gz

$ cd vsearch-2.9.1

$ ./autogen.sh

$ ./configure

$ make

$ make install  # as root or sudo make install


- Rename vsearch to usearch61:

$ mv vsearch-2.9.1 usearch61


- It is necessary, set executable permissions:

$ chmod a+x usearch61


- Move files to /miniconda3/envs/QIIME1/bin

$ mv usearch61  /$path/miniconda3/envs/QIIME191/bin


### Getting started

* Fastq files have to be included in a folder called “Reads”

* A mapping file should be created including information about samples and metadata. For information about its format file, see the webpage: http://qiime.org/documentation/file_formats.html?highlight=mapping.

### How to use the pipeline

Before starting, you should check that your working directory includes:

- Mapping file
- “Reads” folder
- Pipeline files: metagenomics1.pl, metagenomics2.pl and metagenomics_mod.pm

Run the pipeline in the terminal:

$ perl metafora1.pl -t <data type, PE = Paired-End or SE = Single-End> 
-q <Trimmomatic AVGQUAL value> -tr <Trimmomatic TRAILING value>

$ perl metafora2.pl -m <Filename of mapping file (without extension)> -db <database: gg, RDP or SILVA> -level <taxonomic level (L2, L3…L7)> -count <minimum number of count for retaining a OUT> -samples <minimum sample in which the OUT should be present>

### Results

Some examples of graphics generated using this pipeline:

**Taxonomic profile**

Taxonomy profile from analysis using SILVA database:

![](https://github.com/Irvilma/MetaFora/blob/master/Figures/results-tax-silva.png)

**Diversity profile**

2D representation of PCoA results from the unweighted unifrac, weighted and Bray-Curtis index based on Greengenes database. The "Description" title is referred to samples:

![](https://github.com/Irvilma/MetaFora/blob/master/Figures/figure7.png)


[5]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[6]:https://www.java.com/en/
[7]:http://www.usadellab.org/cms/?page=trimmomatic
[8]:https://www.h-its.org/en/research/sco/software/#NextGenerationSequencingSequenceAnalysis
[9]:http://qiime.org/index.html
[10]:https://github.com/torognes/vsearch
[11]: http://qiime.org/install/install.html
[12]:https://anaconda.org/
[13]: https://www.r-project.org/
[14]:https://dplyr.tidyverse.org/
[15]:https://purrr.tidyverse.org/

