#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import csv
import sys
import os

# Write the path to the Galaxy file here!
# e.g. path_to_gff_converter_dir = "/you/files/Code/fml_gff3togtf-5c6f33e20fcc/"
path_to_gff_converter_dir = ""

def placehold(path_in, path_out):
        lookup={}
        with open(path_in, "r") as p:
                data = list(csv.reader(p, delimiter="\t"))
                chroms = list(set([x[0] for x in data if not "#" in x[0]]))
                for i, x in enumerate(chroms):
                        lookup["pl_" + str(i)] = x
                        lookup[x] = "pl_" + str(i)
                for l in data:
                        if not "#" in l[0]: l[0] = lookup[l[0]]
                write(data, path_out)
        return lookup

def unplacehold(path_in, path_out, lookup):
        with open(path_in, "r") as p:
                data = list(csv.reader(p, delimiter="\t"))
                for l in data:
                        if not "#" in l[0]: l[0] = lookup[l[0]]
                write(data, path_out)

def write(data, path_out):
        with open(path_out, "w") as o:
                datawriter = csv.writer(o, delimiter = '\t',quoting = csv.QUOTE_NONE, quotechar='')
                datawriter.writerows(data)

def convert(path_in, path_out):
        function="python "+path_to_gff_converter_dir+"/gff_to_gtf.py " + path_in + " > " + path_out
        subprocess.call([function], shell = True)

def checkForScript():
	if not path_to_gff_converter or not os.path.isfile(path_to_gff_converter+"/gff_to_gtf.py"):
		print "ERROR: Please specify the path to the directory containing the Galaxy gtf_to_gff converter script.\n"
		print "       This can be downloaded from https://toolshed.g2.bx.psu.edu/repository?repository_id=afcb6456d8e300ed\n"
		print "       Or cloned from github from https://github.com/vipints/GFFtools-GX.git"
		print "       The path to the directory should then be placed in the appropriate place in the python file that you've just tried to run :)"
		sys.exit()

if __name__ == '__main__':
        args = sys.argv[1:]
        infile = args[0]
        outfile = args[1]
	checkForScript()
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
