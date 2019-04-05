# Snakemake Barcoding Assembly Pipeline

Pipeline for processing Illumina sequencing data consisting of Diptera COI PCR amplicons. Current release is not
compatible with other multi-fragment COI strategies. 

1) Trims adapters and bases below <20 quality score [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
2) Merges paired end reads into either Fragment A or Fragment B corresponding to PCR amplicons 
3) Removes degenerate primers from the fragments 
4) Assembles reads using [SPAdes](http://cab.spbu.ru/software/spades/)
5) Detects and extracts target contigs
6) Aligns sequences to COI reference (Chrysomya putoria (NCBI accession number NC002697) to correct 5'-3' orientation [Mafft](https://mafft.cbrc.jp/alignment/software/) 
7) Submit sequences to BOLD API for bulk identification [Bold Retriever](https://bold-retriever.readthedocs.io/en/latest/)

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
* Using Bold Retriever after successful run of original pipeline
```bash
snakemake -s bold_retriever_snakefile --use-conda -k
```

## Methodology
Pipeline was designed to handle COI genes amplified in overlapping fragements from thousands of Diptera specimens. 

Reads are trimmed of adaptors and poor quality bases. Paired end reads are merged into single long fragments. 
Each fragment is examined for primers, specifically Fragment A forward primer, and Fragment B reverse primer. 
If Fragment A forward primer is detected, read is trimmed of the degenerate reverse primer. If Fragment B reverse
primer is detected, read is trimmed of degenerate forward primer. Reads are subsequently assembled with SPAdes 
which provides length and coverage information for all contigs. 
Each assembly is examined for a single high coverage contig above 600 basepairs, strongly indicating a successful COI
amplification without contamination. If the single high coverage contig is to short but still above 400 basepairs, 
it is put into a medium quality pool rather than good quality. 

Assemblies that failed, or contain a variety of high coverage contigs, or no contigs above 400 basepairs with 
high coverage are pooled into a folder called problem_fastas

In order to more fully understand the problems facing some of the assemblies, reads can be aligned to a CO1 reference, 
and the corresponding pileups examined using barcoding_snakefile. To examine the pileups of reads against the reference, 
[Tablet](https://ics.hutton.ac.uk/tablet/) can be used on the sorted.bam file, and sorted.bam.bai index file. 

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

## Author
Jackson Eyres \
Bioinformatics Programmer \
Agriculture & Agri-Food Canada \
jackson.eyres@canada.ca
