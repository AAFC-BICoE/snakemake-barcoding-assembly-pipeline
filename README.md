# Snakemake Barcoding Assembly Pipeline

Pipeline for processing Illumina sequencing data consisting of COI PCR amplicons . 

1) Trims adapters and bases below <20 quality score [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
2) Assembles trimmed reads [SPAdes](http://cab.spbu.ru/software/spades/)
3) Detects and extracts target contigs
4) Alignssequences to COI reference (Chrysomya putoria (NCBI accession number NC002697) to correct 5'-3' orientation [Mafft](https://mafft.cbrc.jp/alignment/software/) 
5) Submit Sequences to BOLD for Identification [Bold Retriever](https://bold-retriever.readthedocs.io/en/latest/)

### Prerequisites

* [Conda Installation](https://conda.io/docs/user-guide/install/index.html)
* [Snakemake Installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
* Git

Installing Miniconda + Snakemake
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install -c bioconda -c conda-forge snakemake
```

## Getting Started

Within a working directory:

* Copy Snakefile and pipeline_files directory into working directory
* Create a folder named "fastq" that contains Illumina based raw reads in fastq.gz format 

After initializing a conda enviroment containing snakemake, pipeline can be invoked from within working directory 
```
snakemake --use-conda -k
```

## Built With

* [Python](https://www.python.org/doc/) - Programming language
* [Conda](https://conda.io/docs/index.html) - Package, dependency and environment management
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) - Workflow management system
* [BioPython](https://biopython.org/) - Tools for biological computation
* [Mafft](https://mafft.cbrc.jp/alignment/software/) - Multiple sequence alignment
* [Bold Retriever](https://bold-retriever.readthedocs.io/en/latest/) - Automated BOLD Submission
* [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) - Adaptor and quality trimming
* [SPAdes](http://cab.spbu.ru/software/spades/) - De Novo short read assembler

## Copyright
Government of Canada, Agriculture & Agri-Food Canada

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Citations
* Snakemake  
Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.

* Bold Retriever  
Vesterinen, E. J., Ruokolainen, L., Wahlberg, N., Peña, C., Roslin, T., Laine, V. N., Vasko, V., Sääksjärvi, I. E., 
Norrdahl, K., and Lilley, T. M. (2016) What you need is what you eat? Prey selection by the bat Myotis daubentonii.
Molecular Ecology, 25(7), 1581–1594. doi:10.1111/mec.13564

* Mafft  
Nakamura, Yamada, Tomii, Katoh 2018 (Bioinformatics 34:2490–2492) 
Parallelization of MAFFT for large-scale multiple sequence alignments. 
(describes MPI parallelization of accurate progressive options)

* SPAdes  
Nurk S. et al. (2013) Assembling Genomes and Mini-metagenomes from Highly Chimeric Reads. In: Deng M., Jiang R., 
Sun F., Zhang X. (eds) Research in Computational Molecular Biology. RECOMB 2013. Lecture Notes in Computer Science, 
vol 7821. Springer, Berlin, Heidelberg

* BBTools  
Brian-JGI (2018) BBTools is a suite of fast, multithreaded bioinformatics tools designed for analysis of DNA and RNA 
sequence data.https://jgi.doe.gov/data-and-tools/bbtools/ 

* FASTQC  
Andrews S. (2018). FastQC: a quality control tool for high throughput sequence data. 
Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

* BOLD  
Ratnasingham, S. & Hebert, P. D. N. (2007). BOLD : The Barcode of Life Data System (www.barcodinglife.org).
Molecular Ecology Notes 7, 355–364. DOI: 10.1111/j.1471-8286.2006.01678.x