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

If OrthoFinder has already been run on a set of proteomes and the corresponding CDS nucleotide sequences are available, the `--prep` flag can be used, the input in this case being a genome FASTA file, GTF file, gene CDS sequences, and AA sequences for each species, along with the orthofinder results. This second method is intended to reduce processing time for proteomes that have already been analysed with OrthoFinder. All genomes and sequence files should be supplied in FASTA format. The locations of the genome, gtf, and sequence files should be put in a file in the following format:

In the event that a gtf file contains coordinates not present in the genome fasta file, OrthoFuller will throw an error and will fail to run. Enaure that all chromosome names in the gtf file matxh thoa in the fasta before running.

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
OrthoFiller uses gtf files for both its input and output, due to the superior uniformity of gtf naming and attribute conventions. To convert files from gff3 to gtf format, we recommend using the simple tool fml_gff3togtf, from the Galaxy Tool Shed. This can be found at https://toolshed.g2.bx.psu.edu/repository?repository_id=afcb6456d8e300ed, and is implemented using python:

```
python gff_to_gtf.py infile.gff3 > outfile.gtf
```

Users should note that this tool truncates chromsome names to 15 characters. If this is going to be an issue, you may wish to use placeholder names in the conversion step, for example:

```
file="file_to_convert.gff3"
cut -f1 $file | sort -u | grep -v "#" > $file.placeholder
awk '{print $1"\tTYRANT_"NR}' $file.placeholder > $file.placeholder.lookup
cp $file $file.placeholder.replaced
while read line; do
	echo "placeholding $line"; first=`echo "$line" | cut -f1`; second=`echo "$line" | cut -f2`;
	sed -ri "s/^$first\t/$second\t/g" $file.placeholder.replaced
done < $file.placeholder.lookup
python gff_to_gtf.py $file.placeholder.replaced > $file.placeholder.nearly;
cp $file.tyrant.nearly ${file%gff3}gtf
while read line; do
	echo "unplaceholding $line"; first=`echo "$line" | cut -f1`; second=`echo "$line" | cut -f2`;
	sed -ri "s/^$second\t/$first\t/g" ${file%gff3}gtf
done < $file.placeholder.lookup
rm $file*placeholder*
```

Output File Format
==================
OrthoFiller output can be found in the `results` folder of the specified output directory. For each species, four files are produced:

1. `results/species.newGenes.gff` is a GTF file containing new genes discovered by OrthoFiller;

2. `results/species.results.gff` is a GTF file containing new genes discovered by OrthoFiller as well as all genes from the original annotation;

3. `results/species.newGenes.gff` is a FASTA file containing sequences of the new genes discovered by OrthoFiller;

4. `results/species.results.gff` is a FASTA file containing sequences of the new genes discovered by OrthoFiller as well as sequences of all genes from the original annotation.


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
