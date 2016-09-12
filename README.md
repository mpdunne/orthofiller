# OrthoFinder — identifying missing annotations for evolutionarily conserved genes.

What does OrthoFiller do?
==========
OrthoFiller simultaneously leverages data from multiple species to mutually improve genome annotations. It is designed specifically to address the problem of “missing” genes in sets of predicted genes: that is, to identify those genes that should be present in a genome’s annotation, but which have not been predicted, and whose existence can be verified through comparison with known gene families. OrthoFiller requires only a set of input genome files and a set of corresponding gene GTF files.

https://github.com/mdunne/OrthoFiller

Usage
=====
OrthoFiller runs as a single command that takes as input a tabulated text file containing locations of genomes and GTF files in the following format:

```
#gtf                            genome
/path/to/gtf_species1.gtf       /path/to/genome_species1.fasta
/path/to/gtf_species2.gtf       /path/to/genome_species2.fasta
/path/to/gtf_species3.gtf       /path/to/genome_species3.fasta
etc.
```

OrthoFiller is then run using:

`python OrthoFiller.py -i path/to/genome_locations_file.tdv -o output_folder -c number_of_cores`

If no output folder is specified then OrthoFiller will create one with a generic name. If the number of cores is not specified, OrthoFiller will run using only one core.

If OrthoFinder has already been run on a set of proteomes and the corresponding CDS nucleotide sequences are available, the `--prep` flag can be used, the input in this case being a genome FASTA file, GTF file, gene CDS sequences, and AA sequences for each species, along with the orthofinder results. This second method is intended to reduce processing time for proteomes that have already been analysed with OrthoFinder. All genomes and sequence files should be supplied in FASTA format. The locations of the genome, gtf, and sequence files should be put in a file in the following format:

```
#protein			gtf                            genome				cds
/path/to/aa_species1.fasta	/path/to/gtf_species1.gtf      /path/to/genome_species1.fasta	/path/to/cds_species1.fasta
/path/to/aa_species1.fasta	/path/to/gtf_species2.gtf      /path/to/genome_species2.fasta	/path/to/cds_species2.fasta
/path/to/aa_species1.fasta	/path/to/gtf_species3.gtf      /path/to/genome_species3.fasta	/path/to/cds_species3.fasta
etc.
```

Orthofinder is then run using:
`python OrthoFiller.py --prep -i path/to/genome_locations_file.tdv -o output_folder -c number_of_cores -g path/to/orthofinder_output_orthogroups.csv -s path/to/orthofinder_output_unassigned_genes.csv`

An example skeleton bash script for running OrthoFiller is included as runOrthoFil.sh. Paths for input files and for relevant packages must be added manually.

Output File Format
==================
OrthoFiller output can be found in the `results` folder of the specified output directory. For each species, four files are produced:

**1) results/species.newGenes.gff is a GTF file containing new genes discovered by OrthoFiller;
**2) results/species.results.gff is a GTF file containing new genes discovered by OrthoFiller as well as all genes from the original annotation;
**3) results/species.newGenes.gff is a FASTA file containing sequences of the new genes discovered by OrthoFiller;
**4) results/species.results.gff is a FASTA file containing sequences of the new genes discovered by OrthoFiller as well as sequences of all genes from the original annotation.


Installing Dependencies
=======================
OrthoFiller is written to run on linux and requires the following to be installed and in the system path:

1. Python 2.7 together with the scipy libraries stack (If you don't have python 2.7 you can still use the standalone binaries) 

2. BedTools Suite

3. OrthoFinder

4. MAFFT

5. R with Gamlss package

6. Augustus

7. Hmmer

python and scipy
----------------
Up-to-date and clear instructions are provided here: http://www.scipy.org/install.html, be sure to chose a version using python 2.7. As websites can change, an alternative is to search online for "install scipy".

BedTools
------
The BedTools suite can be downloaded from https://github.com/arq5x/bedtools2/releases

Augustus
--------
Augustus can be downloaded from http://augustus.gobics.de/. In addition to being in the system path, the Augustus config path must be explicitly set in a variable before running OrthoFiller. This can be done by executing, for example, `export AUGUSTUS_CONFIG_PATH=/path/to/augustus/augustus-x.x.x/config`.

Hmmer
-----
HMMER can be downloaded from http://hmmer.org/.
