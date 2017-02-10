# OrthoFiller — identifying missing annotations for evolutionarily conserved genes.

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

If OrthoFinder has already been run on a set of proteomes and the corresponding CDS nucleotide sequences are available, the `--prep` flag can be used, the input in this case being a genome FASTA file, GTF file, gene CDS sequences, and AA sequences for each species, along with the orthofinder results. This second method is intended to reduce processing time for proteomes that have already been analysed with OrthoFinder. All genomes and sequence files should be supplied in FASTA format. The locations of the genome, GTF, and sequence files should be put in a file in the following format:

In the event that a GTF file contains coordinates not present in the genome fasta file, OrthoFuller will throw an error and will fail to run. Enaure that all chromosome names in the GTF file matxh thoa in the fasta before running.

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

Obtaining GTFs from GFF3 files
==============================
OrthoFiller uses GTF files for both its input and output, due to the superior uniformity of GTF naming and attribute conventions. To convert files from GFF3 to GTF format, we recommend using the simple tool fml_gff3togtf, from the Galaxy Tool Shed. This can be found at https://toolshed.g2.bx.psu.edu/repository?repository_id=afcb6456d8e300ed, and is implemented using python:

```
python gff_to_gtf.py infile.gff3 > outfile.gtf
```

Users should note that this tool truncates chromsome names to 15 characters. If this is going to be an issue, you may wish to use placeholder names in the conversion step, for example, by using this python script:

```
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import csv
import sys
import os

def placehold(path_in, path_out):
        lookup={}
        with open(path_in, "r") as p:
                data = list(csv.reader(p, delimiter="\t"))
                chroms = list(set([x[0] for x in data]))
                for i, x in enumerate(chroms):
                        lookup["pl_" + str(i)] = x
                        lookup[x] = "pl_" + str(i)
                for l in data: l[0] = lookup[l[0]]
                write(data, path_out)
        return lookup

def unplacehold(path_in, path_out, lookup):
        with open(path_in, "r") as p:
                data = list(csv.reader(p, delimiter="\t"))
                for l in data: l[0] = lookup[l[0]]
                write(data, path_out)

def write(data, path_out):
        with open(path_out, "w") as o:
                datawriter = csv.writer(o, delimiter = '\t',quoting = csv.QUOTE_NONE, quotechar='')
                datawriter.writerows(data)

def convert(path_in, path_out):
        function="python PATH_TO_SCRIPT/gff_to_gtf.py " + path_in + " > " + path_out
        subprocess.call([function], shell = True)

if __name__ == '__main__':
        args = sys.argv[1:]
        infile = args[0]
        outfile = args[1]
        # Placehold and convert
        print("placeholding...")
        placeheld = infile + ".tmp"
        lookup = placehold(infile, placeheld)
        placeheldout = placeheld + ".plh"
        print("converting....")
        convert(placeheld, placeheldout)
        print("unplaceholding...")
        unplacehold(placeheldout, outfile, lookup)
        # Clean up
        os.remove(placeheld)
        os.remove(placeheldout)
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
