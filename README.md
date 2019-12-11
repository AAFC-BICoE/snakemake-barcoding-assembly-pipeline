# Snakemake Barcoding Assembly Pipeline

Pipeline for processing Illumina sequencing data consisting of AAFC Diptera COI PCR amplicons. Current release is not
compatible with other multi-fragment COI strategies. 

1) Trims adapters and bases below 10 quality score [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
2) Merges paired end reads into either Fragment A or Fragment B corresponding to PCR amplicons 
3) Removes degenerate primers from the fragments 
4) Dereplicates merged reads using VSEARCH  [VSearch](https://github.com/torognes/vsearch)
5) Merges top two merged reads using Emboss Merger [Emboss](http://emboss.open-bio.org/)
6) Evaluates results and generates a multi-fasta

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

* Create and enter a working directory:
* Clone repository
```bash
git clone https://github.com/AAFC-BICoE/snakemake-barcoding-assembly-pipeline.git .
```
* Create a folder named "fastq" and populate with COI Illumina reads in fastq.gz format 
* Initialize conda environment containing snakemake
```bash
source ~/miniconda3/bin/activate
```
* Invoke pipeline from within working directory. Adjust cores to suit computer 
```
snakemake --use-conda -k --cores 32
```
* Alternative pipeline to map reads to COI reference gene
```
snakemake -s barcoding_snakefile --use-conda -k --cores 32
```

## Methodology
Pipeline was designed to handle COI genes amplified in overlapping fragments from thousands of Diptera specimens. 

Reads are trimmed of adaptors and poor quality bases. Paired end reads are merged into single long fragments. 
Each fragment is examined for primers, specifically Fragment A forward primer, and Fragment B reverse primer. 
If Fragment A forward primer is detected, read is trimmed of all Fragment A primers. If Fragment B reverse
primer is detected, read is trimmed of all Fragment B primers. Reads are then de-replicated and abundance calculated. 
Ideal scenario is two very high abundance reads representing Fragment A and B, with minimal or no other reads. The top
two reads are merged together to form a single CO1 sequence. 

Curated CO1 sequences that contain <10% contamination/error reads, have a minimum of 40 combined reads
that pass all filtering steps, and are the correct 646 bp length are placed in Curated_Barcodes.fasta
Quality checks for all sequences are listed in Summary_Output.csv


## Built With

* [Python](https://www.python.org/doc/) - Programming language
* [Conda](https://conda.io/docs/index.html) - Package, dependency and environment management
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) - Workflow management system
* [BioPython](https://biopython.org/) - Tools for biological computation

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

* BBTools  
Brian-JGI (2018) BBTools is a suite of fast, multithreaded bioinformatics tools designed for analysis of DNA and RNA 
sequence data.https://jgi.doe.gov/data-and-tools/bbtools/ 

* BOLD  
Ratnasingham, S. & Hebert, P. D. N. (2007). BOLD : The Barcode of Life Data System (www.barcodinglife.org).
Molecular Ecology Notes 7, 355–364. DOI: 10.1111/j.1471-8286.2006.01678.x

* VSearch  
Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. 
PeerJ 4:e2584. doi: 10.7717/peerj.2584

* Emboss  
Rice P., Longden I. and Bleasby A. EMBOSS: The European Molecular Biology Open Software Suite. 
Trends in Genetics. 2000 16(6):276-277

## Author
Jackson Eyres \
Bioinformatics Programmer \
Agriculture & Agri-Food Canada \
jackson.eyres@canada.ca
