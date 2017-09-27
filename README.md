![alt text](http://empede.co.uk/imgrepos/OrthoFiller_head.png? "OrthoFiller logo")

What does OrthoFiller do?
==========
**OrthoFiller** simultaneously leverages data from multiple species to **mutually improve genome annotations**. It is designed specifically to address the problem of **“missing” genes** in sets of predicted genes: that is, to identify those genes that should be present in a genome’s annotation, but which have not been predicted, and whose existence can be verified through comparison with known gene families. OrthoFiller requires only a set of input genome files and a set of corresponding gene GTF files.

![alt text](http://empede.co.uk/imgrepos/Workflow-03.png "OrthoFiller workflow")


### Recent Updates

**2017-06-07** - Speed updates:
* Improved speed of early protein fasta file manipulation.
* Improved speed of nucleotide alignment acquisition

**2017-05-31** - OrthoFiller 1.1.0 released. Including:
* Introduction of target and reference species option
* Ability to search chromosomes individually to limit memory usage
* Bug fixes

For more details, see the OrthoFiller paper:

[Dunne, M.P. and Kelly, S. (2017) OrthoFiller: utilising data from multiple species to improve the completeness of genome annotations, BMC Genomics 18:390](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3771-x)

[<img src="http://empede.co.uk/imgrepos/OrthoFiller_paper.jpg?">](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3771-x)

Github link:

https://github.com/mdunne/OrthoFiller

Usage
=====

### Default usage

OrthoFiller runs as a single command that takes as input a tabulated text file containing locations of genomes and GTF files in the following format:

```
#gtf                            genome
/path/to/gtf_species1.gtf       /path/to/genome_species1.fasta
/path/to/gtf_species2.gtf       /path/to/genome_species2.fasta
/path/to/gtf_species3.gtf       /path/to/genome_species3.fasta
etc.
```

OrthoFiller is then run using:

`python OrthoFiller.py -i path/to/genome_locations_file.tdv -o output_folder -c num_cores`

If no output folder is specified then OrthoFiller will create one with a generic name. If the number of cores is not specified, OrthoFiller will run using only one core: this is not recommended as it will cause OrthoFiller to take a considerably long time. It is recommended that at least the same number of cores as number of species are used, and preferably at least double the number.

### Reference species

If you do not wish to search for new genes in every one of your selected species, you may divide them into *target* and *reference* species. Target species are specified in a file given by the `-i` option, and reference species are specified in the file given by the `-r` option.

`python OrthoFiller.py -i path/to/target_species.tdv -r path/to/reference_species.tdv -o output_folder -c num_cores`

The `target_species.tdv` file should *only* contain the species whose genomes you would like to search, and the `reference_species.tdv` file should contain *only* species which you *do not* wish to search.

### Pre-specified orthogroup and cds files

If OrthoFinder has already been run on a set of proteomes and the corresponding CDS nucleotide sequences are available, the `--prep` flag can be used, the input in this case being a genome FASTA file, GTF file, gene CDS sequences, and AA sequences for each species, along with the orthofinder results. This second method is intended to reduce processing time for proteomes that have already been analysed with OrthoFinder. All genomes and sequence files should be supplied in FASTA format. The locations of the genome, GTF, and sequence files should be put in a file in the following format:

```
#protein			gtf                            genome				cds
/path/to/aa_species1.fasta	/path/to/gtf_species1.gtf      /path/to/genome_species1.fasta	/path/to/cds_species1.fasta
/path/to/aa_species1.fasta	/path/to/gtf_species2.gtf      /path/to/genome_species2.fasta	/path/to/cds_species2.fasta
/path/to/aa_species1.fasta	/path/to/gtf_species3.gtf      /path/to/genome_species3.fasta	/path/to/cds_species3.fasta
etc.
```

Orthofiller is then run using:
`python OrthoFiller.py --prep -i path/to/genome_locations_file.tdv -o output_folder -c number_of_cores -g path/to/orthofinder_output_orthogroups.csv -s path/to/orthofinder_output_unassigned_genes.csv`

An example skeleton bash script for running OrthoFiller is included as runOrthoFil.sh. Paths for input files and for relevant packages must be added manually.

### Split by chromosome
For large genomes, HMMER sometimes runs into memory problems when applying HMM searches to the genome. If you run into such problems, you may wish to try the `--split` option for OrthoFiller, which splits the inputted genome files up and analyses each of their chromosomes individually. Although this limits memory usage, it does increase the runtime of OrthoFiller, especially for genomes with large numbers of chromosomes/contigs/scaffolds.

### OrthoFiller likes well-made GTF files.

In the event that a GTF file contains coordinates not present in the genome fasta file, OrthoFiller will throw an error and will fail to run. Ensure that all chromosome names in the GTF file match those in the fasta before running.

### OrthoFiller likes clean FASTA files.

The current iteration of OrthoFiller requires that chromosome names in the genome FASTA files consist of unspaced IDs only (i.e. no description lines or information other than the name). This will be fixed in future versions. In the meantime, you may wish simply to modify the names of your chromosome entries using:

```
sed -r "s/>([^ ]+*) .*/>\1/g" genome.fa > genome_clean.fa
```

Obtaining GTFs from GFF3 files
==============================
OrthoFiller uses GTF files for both its input and output, due to the superior uniformity of GTF naming and attribute conventions. To convert files from GFF3 to GTF format, we recommend using the simple tool fml_gff3togtf, from the Galaxy Tool Shed. This can be found at https://toolshed.g2.bx.psu.edu/repository?repository_id=afcb6456d8e300ed, and is implemented using python:

```
python gff_to_gtf.py infile.gff3 > outfile.gtf
```

Users should note that this tool truncates chromsome names to 15 characters. If this is going to be an issue, a wrapper for this script can be found in the utils directory in this repository (https://github.com/mpdunne/orthofiller/blob/master/utils/gff_to_gtf_safe.py). The above Galaxy tool should be downloaded first, and the path to its directory should be included in the appropriate place at the top of the `gff_to_gtf_safe.py` file. The full script can then be run as

```
python gff_to_gtf_safe.py infile.gff3 outfile.gtf
```

Note that, in order to function properly, the above conversion script requires that entries in the GFF3 input file are well-formed: that is they contina gene, mRNA, CDS, and exon entries for each gene. Ideally ensure that each GFF3 entry has each of these attributes before proceeding. Alternatively, if you simply wish to remove incomplete entries from your GFF3 file, you can use the `clean_gff.py` script, also included in the utils directory of this repository. The usage for this script is:

```
python clean_gff.py infile.gff3 infile_clean.gff3
```



Output File Format
==================
OrthoFiller output can be found in the `results` folder of the specified output directory. For each species, four files are produced:

1. `results/species.newGenes.gtf` is a GTF file containing new genes discovered by OrthoFiller;

2. `results/species.results.gtf` is a GTF file containing new genes discovered by OrthoFiller as well as all genes from the original annotation;

3. `results/species.newSequences.aa.fasta` is an amino acid FASTA file containing sequences of the new genes discovered by OrthoFiller;

4. `results/species.results.aa.fasta` is an amino acid FASTA file containing sequences of the new genes discovered by OrthoFiller as well as sequences of all genes from the original annotation.


Installing Dependencies
=======================
OrthoFiller is written to run on linux and requires the following to be installed and in the system path:

1. Python 2.7 together with the scipy and BioPython libraries 

2. BedTools Suite

3. OrthoFinder

4. R with Gamlss package

5. Augustus

6. Hmmer

### python and scipy

Up-to-date and clear instructions are provided here: http://www.scipy.org/install.html, be sure to chose a version using python 2.7. As websites can change, an alternative is to search online for "install scipy".

### BedTools

The BedTools suite can be downloaded from https://github.com/arq5x/bedtools2/releases

### Augustus

Augustus can be downloaded from http://augustus.gobics.de/. In addition to being in the system path, the Augustus config path must be explicitly set in a variable before running OrthoFiller. This can be done by executing, for example, `export AUGUSTUS_CONFIG_PATH=/path/to/augustus/augustus-x.x.x/config`.

### Hmmer

HMMER can be downloaded from http://hmmer.org/.

### OrthoFinder

OrthoFinder can be downloaded here: https://github.com/davidemms/OrthoFinder

Other useful software
=====================

### Alan

Alan is an in-terminal command-line tool for viewing alignments without the need for a GUI.

https://github.com/mpdunne/alan
