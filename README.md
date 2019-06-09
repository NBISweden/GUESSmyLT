# GUESSmyLT
Software to guess the RNA-Seq library type of paired and single end read files using mapping and gene annotation.  

## Table of contents

* [Background](#background)
* [Dependencies](#Dependencies)
* [Installation](#Installation)
  * [Installation with conda](#installation-with-conda)
  * [Installation with pip](#installation-with-pip)
  * [Installation with git](#installation-with-git)
  * [Check installation](#check-installation)
* [Result](#result)
* [Usage](#usage)
  * [File formats](#file-formats)
  * [Supported header formats](#supported-header-formats)
  * [Supported interleaved formats](#supported-interleaved-formats)
  * [Example commands](#example-commands)
    * [If you have no other information just the reads](#if-you-have-no-other-information-just-the-reads)
    * [If you only have a reference genome](#if-you-only-have-a-reference-genome)
    * [If you have reference genome and annotation](#if-you-have-reference-genome-and-annotation)
    * [If you have reference genome and mapped reads](#if-you-have-reference-genome-and-mapped-reads)
    * [If you have reference genome and annotation and mapped reads](#if-you-have-reference-genome-and-annotation-and-mapped-reads)
    * [If you only have transcript sequences](#if-you-only-have-transcript-sequences)
    * [If you have transcript sequences and annotation](#if-you-have-transcript-sequences-and-annotation)
    * [Other examples](#other-examples)
  * [Output](#output)
  * [Parameters](#parameters)
* [Overview of the pipeline](#overview-of-the-pipeline)
* [Overview of the different library types](#Overview-of-the-different-library-types)
* [Library prep methods](#library-prep-methods)
* [External resources](#external-resources)
* [Known issues](#known-issues)
* [TO DO](#to-do)
* [Citation](#citation)

## Background
The choice of RNA-Seq library type defines the read orientation of the sequencing and the order in which the strands of cDNA are sequenced, which means that RNA-Seq reads from different library types can differ significantly. The information regarding library type can be very useful for reads to be assembled into a transcriptome or mapped to a reference assembly. This is because the library type can help to discern where in the transcriptome shorter ambiguous reads belong by using the read’s relative orientation and from which strand it was sequenced. Unfortunately, this information regarding the library type used is not included in sequencing output files and is usually lost before the assembly of the data. Even when working with RNA-Seq data from public repositories there is no guarantee that the library type information is correct or that it exists at all. This is what GUESSmyLT aims to fix by looking at how reads map to a reference and together with gene annotation guess which library was used to generate the data.

## Dependencies:  
Developed for Unix systems. Depending installation approach more or less dependencies will be installed automatically. Check the installation paragraph.

Python and libraries:
 * Python >3
 * biopython (1.67)
 * bcbio-gff (0.6.4) - handling gff annotation
 * pysam (0.15.1) - handling mapped reads

Other programs:
 * Snakemake (5.4.0) - Workflow management
 * BUSCO (3.0.2) - Gene annotation
 * Bowtie2 (2.3.4.3) - Mapping
 * Trinity (2.8.4) - Reference assembly

Others:  
  - Prokaryote and eukaryote BUSCO datasets (from https://busco.ezlab.org) included in the package.

## Installation  

#### Installation with conda (in preparation): 

```bash

```

#### Installation with pip:  

Installation using **pip** will not install BUSCO, Bowtie2 and Trinity. These external programs can be installed using conda.

```bash
pip install GUESSmyLT
```

#### Installation with git:  

Installation using **git** will not install BUSCO, Bowtie2 and Trinity. These external programs can be installed using conda.

Clone the repository and move to the folder:

```bash
git clone https://github.com/NBISweden/GUESSmyLT.git
cd GUESSmyLT/
```

Launch the installation either:
```bash
python setup.py install
```

Or if you do not have administrative rights on your machine:

```bash
python setup.py install --user
```

#### Check installation

Executing:
```bash
GUESSmyLT
```

or

```bash
GUESSmyLT -h
```

to display help.

There is also an example run that takes roughly 5 mins. A folder called GUESSmyLT_example_out will be created in the working directory:
```bash
GUESSmyLT-example
```

## Result
The results are printed as stdout and to a result file. One example of a result would be:
```bash
Results of paired library inferring of reads 4_mapped_r1.sub.100000 on ref 4:

Library type    Reads     Percent     Vizualization according to firststrand

    fr_first     4019       47.2%     3' ----------<==1== 5'
                                      5' ==2==>---------- 3'


   fr_second     4454       52.3%     3' ----------<==2== 5'
                                      5' ==1==>---------- 3'


    rf_first       21        0.2%     3' ----------==1==> 5'
                                      5' <==2==---------- 3'


   rf_second       19        0.2%     3' ----------==2==> 5'
                                      5' <==1==---------- 3'


    ff_first        5        0.1%     3' ----------<==1== 5'
                                      5' <==2==---------- 3'


   ff_second        2        0.0%     3' ----------==2==> 5'
                                      5' ==1==>---------- 3'


   undecided        1        0.0%     3' -------??------- 5'
                                      5' -------??------- 3'


Roughly 50/50 split between the strands of the same library orientation should be interpreted as unstranded.

```

Based on the orientations of the reads we would assume that the library type is fr-unstranded as there is roughly a 50-50 split between fr-first and fr-second.

## Usage

### File formats
Read files: .fastq
Mapping:    .bam
Reference:  .fa


### Supported header formats
Tested for Old/New Illumina headers and downloads from SRA.
Should work, but not tested for all fastq header formats at: https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/

### Supported interleaved formats
If headers are in Old/New Illumina or if reads are alternating.  

Old Illumina: @HWUSI-EAS100R:6:73:941:1973#0/1  
New Illumina: @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG  
Alternating:  
&nbsp; &nbsp; &nbsp; @read1 (first mate)  
&nbsp; &nbsp; &nbsp; ..  
&nbsp; &nbsp; &nbsp; @read1 (second mate)  
&nbsp; &nbsp; &nbsp; ..  
&nbsp; &nbsp; &nbsp; @read2 (first mate)  
&nbsp; &nbsp; &nbsp; ..   
&nbsp; &nbsp; &nbsp; @read2 (second mate)  
&nbsp; &nbsp; &nbsp; ..  

### Example commands
**In top of your fastq RNA-Seq read file(s) (compressed or uncompressed):**  

#### If you have no other information just the reads
Example with paired reads in eukaryote.  
```bash
GUESSmyLT --reads read_1.fastq read_2.fastq
```
#### If you only have a reference genome  
Example with paired reads in eukaryote.  
```bash
GUESSmyLT --reads read_1.fastq read_2.fastq --reference ref.fa --mode genome --organism euk
```

#### If you have reference genome and annotation
Example with paired reads in eukaryote.  
```bash
GUESSmyLT --reads read_1.fastq read_2.fastq --reference ref.fa --mode genome --annotation annotation.gff --organism euk
```

#### If you have reference genome and mapped reads
Example with paired reads in eukaryote.
```bash
GUESSmyLT --reads read_1.fastq read_2.fastq --reference ref.fa --mode genome --mapped mapped.bam --organism euk
```

#### If you have reference genome and annotation and mapped reads
Example with paired reads in eukaryote.
```bash
GUESSmyLT --reads read_1.fastq read_2.fastq --reference ref.fa --mode genome --mapped mapped.bam --annotation annotation.gff --organism euk
```

#### If you only have transcript sequences
/!\ not yet implemented  (use genome mode instead it should work anyway)
Example with paired reads in eukaryote.  
```bash
GUESSmyLT --reads read_1.fastq read_2.fastq --reference ref.fa --mode transcriptome --organism euk
```

#### If you have transcript sequences and annotation
/!\ not yet implemented (use genome mode instead it should work anyway)
Example with paired reads in eukaryote.  (The annotation has to be the annotation within the trascriptome not the genome)
```bash
GUESSmyLT --reads read_1.fastq read_2.fastq --reference ref.fa --mode transcriptome --annotation annotation.gff --organism euk
```

#### Other examples

##### Paired end reads and reference with specified subsampled reads. Output directed to existing directory.
```bash
GUESSmyLT --reads read_1.fastq read_2.fastq --organism pro --reference ref.fa --subsample 100000 --output my_output/
```

##### Single end reads and prokaryote genome
```bash
GUESSmyLT --reads reads.fastq --reference ref.fa --organism pro
```

#### Interleaved paired reads
```bash
cd GUESSmyLT/
python3 GUESSmyLT.py --reads reads.fastq --reference ref.fa --organism euk
```
### Output
GUESSmyLT will print the result in the command line as well as write it to a file:
```bash
[output_dir]/result_[read_name]on_[refname].txt
```
Results from intermediate steps, such as the mapping from Bowtie2 or annotation from BUSCO are saved in
```bash
[output_dir]/intermediate_data/
```

## Parameters

### Mandatory

| Parameter | Input | Description |
| --- | --- | --- |
| --reads | .fastq file(s) | Full path(s) to RNA-Seq read file(s). Can be compressed or uncompressed. Order is not important. Can handle two paired end read files, one interleaved read file and single end read file. |
| --organism | euk or pro | Eukaryote or prokaryote (euk/pro) is an option needed for the BUSCO annotation. |



### Optional
| Parameter | Input | Description |
| --- | --- | --- |
| --subsample | Even integer | Number of reads that will be used for subsampling. |
| --reference | .fa file | Full paths to reference genome/transcriptome for mapping reads to (nucleotide fasta file). |
| --mode | genome or transcriptome | When no annotation is provided, tells the programm if the reference fasta file has to be considered as a genome or a transcriptome in order to use BUSCO properly. |
| --threads | Integer | Number of threads to use. |
| --memory | Number of GB ex: 10G | Maximum memory that can be used in GB. |
| --annotation | .gff file | Full path to annotation file for skipping BUSCO step. NOT DEVELOPED YET. |
| --mapped | Sorted .bam file | Full path to mapped read file for skipping Bowtie2 step. NOT DEVELOPED YET. |
| --output | File path | Full path to result file. If left out files will be written to working directory. |

## Overview of the pipeline
![alt text](https://github.com/NBISweden/GUESSmyLT/blob/master/Overview_of_GUESSmyLT.png "Pipeline of GUESSmyLT")  
GUESSmyLT uses Snakemake to build the pipeline it needs in order to predict the library type. Required arguments are organism (euk/pro) and reads (read file(s) in fastq format). Reference (genome or transcriptome in .fasta format) is optional, and if it is not provided, Trinity will first be executed to create a De novo assembly of the reads. Next, BUSCO is used for annotation. This is also a QC step because BUSCO looks for core genes, so called BUSCOs, in the reference. If they cannot be found, it indicates that the reference has bad quality and therefore the pipeline will terminate. If BUSCOs are found, the process continues with mapping the reads to the reference using Bowtie2. The mapping is done with unstranded option so that the reads can be mapped on both the strands and in both directions. Finally, the mapping and annotation is used for inference, which is done with a python script and the library type is returned.
On top of Snakemake, we have a python script, GUESSmyLT.py. Its purpose is to handle user arguments by:
1.	Checking that arguments are correct, files exists and are in correct format.
2.	Telling Snakemake what files exist by updating the config file.
3.	Executing snakemake.

The Snakefile subsample handles preparation of the readfiles:
1.	Subsamples reads into new read files that are used in the analysis. This makes GUESSmyLT faster and protects the original files from being modified.
2.	Modifying files:
a.	Changes read files that are in wrong format. Trinity and Pysam can only handle old Illumina format: @read_ID/pair#, where pair# is 1 or 2. They do not work with whitespaces, punctutations nor undescrores. Therefore, the script makes sure that the headers are converted into the correct format.
b.	Deinterleaves paired end read files if they are interleaved.

## Overview of the different library types:

![alt text](https://github.com/NBISweden/GUESSmyLT/blob/master/library_types.jpg)

## Library prep methods:

| kit | Description | Paired | Stranded | Strand according to mRNA | Strand according to `first strand`|
| --- | --- | --- | --- | --- | ---
| TruSeq RNA Sample Prep kit  | | yes | No | | fr-unstranded |
| SMARTer ultralow RNA protocol | | yes | No | | fr-unstranded |
| All dUTP methods, NSR, NNSR | | yes | Yes | RF | fr-firststrand
| TruSeq Stranded Total RNA Sample Prep Kit | | yes | Yes | RF | fr-firststrand
| TruSeq Stranded mRNA Sample Prep Kit | | yes | Yes | RF | fr-firststrand
| NEB Ultra Directional RNA Library Prep Kit | | yes | Yes | RF | fr-firststrand
| Agilent SureSelect Strand-Specific | | yes | Yes | RF | fr-firststrand
| Directional Illumina (Ligation) | | yes | Yes | FR | fr-secondstrand
| Standard SOLiD | | Yes | yes | FR | fr-secondstrand
| ScriptSeq v2 RNA-Seq Library Preparation Kit | | yes | Yes | FR | fr-secondstrand
| SMARTer Stranded Total RNA | | yes | Yes | FR | fr-secondstrand
| Encore Complete RNA-Seq Library Systems  | | yes | Yes | FR | fr-secondstrand
| NuGEN SoLo  | | yes | Yes | FR | fr-secondstrand
| Illumina ScriptSeq| |  yes | Yes | FR | fr-secondstrand
| SOLiD mate-pair protocol | | | | | ff

--rf orientation are produced using the Illumina mate-pair protocol?

## External resources:

[https://chipster.csc.fi/manual/library-type-summary.html](https://chipster.csc.fi/manual/library-type-summary.html)
[https://galaxyproject.org/tutorials/rb_rnaseq/](https://galaxyproject.org/tutorials/rb_rnaseq/)
[http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html](http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html)
[https://sailfish.readthedocs.io/en/master/library_type.html](https://sailfish.readthedocs.io/en/master/library_type.html)
[https://rnaseq.uoregon.edu](https://rnaseq.uoregon.edu)
[https://www.researchgate.net/post/What_is_the_difference_between_strand-specific_and_not_strand-specific_RNA-seq_data](https://www.researchgate.net/post/What_is_the_difference_between_strand-specific_and_not_strand-specific_RNA-seq_data)

## Known issues
1) Complains about gzip broken pipe when subsampling with compressed files (but works anyway).  
2) BUSCO sometimes looses the config path. Fix manually in terminal:
```bash
export AUGUSTUS_CONFIG_PATH=~/miniconda3/pkgs/augustus-3.2.3-boost1.60_0/config
```
3) BUSCO might not find any core genes. Fix by using more reads or by providing reference.  
4) Mapping, annotation, assembly or the entire pipeline is skipped.  This is most likely due to the fact that Snakemake checks which output files need to be generated and from there only performs the necessary steps of the pipeline. The result of this is that is you already have a .bam file, BUSCO/Trinity output folder or a result .txt file for the reads Snakemake will skip steps  
5) Installing Trinity for mac via Conda will give you a version from 2011 that doesn't work. Install using Homebrew instead.
## TO DO
	* Make BIOCONDA package for easy access. (Maybe snakemakes --use-conda)  
	* Add Travis using example data provided as reference.  
	* Look more into why some reads get undecided orientation. This is when a read's mate cannot be found and is probably due to a read is at the end of a gene and its mate is outside of the selected region.  

## Citation
If you use GUESSmyLT in your work, please cite us.
DOI will come soon

## Author
Berner Wik E.<sup>\*,1</sup>, Olin H.<sup>\*,1</sup>, Vigetun Haughey C.<sup>\*,1</sup>, Lisa Klasson<sup>1</sup>, Jacques Dainat<sup>2,3</sup> 

<sup>*</sup>These authors contributed equally to the work.</br>
<sup>1</sup>Molecular Evolution, Department of Cell and Molecular Biology, Uppsala University, 75124 Sweden.</br>
<sup>2</sup>National Bioinformatics Infrastructure Sweden (NBIS), SciLifeLab, Uppsala Biomedicinska Centrum (BMC), Husargatan 3, S-751 23 Uppsala, SWEDEN.</br>
<sup>3</sup>IMBIM - Department of Medical Biochemistry and Microbiology, Box 582, S-751 23 Uppsala, SWEDEN.</br>

