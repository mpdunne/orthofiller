#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2017 Michael Dunne
#
# This program (OrthoFiller) is distributed under the terms of the GNU General Public License v3
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# For any enquiries send an email to Michael Dunne
# michael.dunne@worc.ox.ac.uk

__author__ = "Michael Dunne"
__credits__ = "Michael Dunne, Reed Roberts, David Emms, Steve Kelly"


########################################################
################ Safely import packages ################
########################################################
errors=[]
try:
	import csv
except ImportError as e:
	errors.append(e)
try:
	import re
except ImportError as e:
	errors.append(e)
try:
	import os
except ImportError as e:
	errors.append(e)
try:
	import sys
except ImportError as e:
	errors.append(e)
try:
	import itertools
except ImportError as e:
	errors.append(e)
try:
	import Bio
except ImportError as e:
	errors.append(e)
try:
	import subprocess
except ImportError as e:
	errors.append(e)
try:
	import multiprocessing
except ImportError as e:
	errors.append(e)
try:
	import datetime
except ImportError as e:
	errors.append(e)
try:
	import tempfile
except ImportError as e:
	errors.append(e)
try:
	import random
except ImportError as e:
	errors.append(e)
try:
	import string
except ImportError as e:
	errors.append(e)
try:
	import commands
except ImportError as e:
	errors.append(e)
try:
	import errno
except ImportError as e:
	errors.append(e)
try:
	from Bio import SeqIO
except ImportError as e:
	errors.append(e)
try:
	from Bio.Align.Applications import MafftCommandline
except ImportError as e:
	errors.append(e)
try:
	from Bio.SeqRecord import SeqRecord
except ImportError as e:
	errors.append(e)
try:
	import argparse
except ImportError as e:
	errors.append(e)
try:
	import time
except ImportError as e:
	errors.append(e)
try:
	import glob
except ImportError as e:
	errors.append(e)

if errors:
	print("Missing modules :(\nThe following module errors need to be resolved before running OrthoFiller:")
	for error in errors:
		print("-- " + str(error))
	sys.exit()

########################################################
################ Check Linux Packages ##################
########################################################

def callFunction(str_function):
	"""Call a function in the shell
	"""
	subprocess.call([str_function], shell = True)

def callFunctionQuiet(str_function):
	"""Call a function in the shell, but suppress output.
	"""
	#callFunction(str)
	with open(os.devnull, 'w') as FNULL:
		subprocess.call([str_function], shell = True, stdout=FNULL, stderr=subprocess.STDOUT)

def CanRunCommand(command, qAllowStderr = False):
	sys.stdout.write("Test can run \"%s\"" % command)
	capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout = [x for x in capture.stdout]
	stderr = [x for x in capture.stderr]
	if len(stdout) > 0 and (qAllowStderr or len(stderr) == 0):
		print(" - ok")
		return True
	else:
		print(" - failed")
		return False

def CanRunMinusH(package, packageFormatted):
	if CanRunCommand(package + " -h"):
		return True
	else:
		print("ERROR: Cannot run " + packageFormatted)
		print("    Please check "+ packageFormatted +" is installed and that the executables are in the system path\n")
		return False

def CanRunHelp(package, packageFormatted):
	if CanRunCommand(package + " -help"):
                return True
        else:
                print("ERROR: Cannot run " + packageFormatted)
                print("    Please check "+ packageFormatted +" is installed and that the executables are in the system path\n")
                return False


def CanRunMan(package, packageFormatted):
        if CanRunCommand("man " + package):
                return True
        else:
                print("ERROR: Cannot run " + packageFormatted)
                print("    Please check "+ packageFormatted +" is installed and that the executables are in the system path\n")
                return False

def CanRunBlank(package, packageFormatted):
        if CanRunCommand(package):
                return True
        else:
                print("ERROR: Cannot run " + packageFormatted)
                print("    Please check "+ packageFormatted +" is installed and that the executables are in the system path\n")
                return False

def CanRunOrthoFinder():
	if "ORTHOFINDER_DIR" in os.environ:
		return CanRunCommand("python " + os.environ["ORTHOFINDER_DIR"] + "/orthofinder.py -h")	
        elif os.path.isfile(os.path.dirname(os.path.abspath(__file__)) + "/orthofinder.py"):
		return CanRunCommand("python " + os.path.dirname(os.path.abspath(__file__)) + "/orthofinder.py -h ")
        elif os.path.isfile("orthofinder.py"):
                return CanRunCommand("python orthofinder.py -h")
        else:
                return CanRunCommand("orthofinder -h")	

def CanRunAwk():
	sys.stdout.write("Test can run \"awk '{print 4 + 5}'\"")
	path_tf = tempfile.mktemp()
	with open(path_tf, "w") as f:
		f.write("a")
	a=commands.getstatusoutput("awk '{print 4 + 5}' " + path_tf)[1]
	os.remove(path_tf)
	if a == "9":
		print(" - ok")
		return True
	else:
		print(" - failed")
		return False

def CanRunMafft():
	sys.stdout.write("Test can run \"mafft\"")
	path_tf = tempfile.mktemp()
	with open(path_tf, "w") as f:
		f.write(">1\nAAABBBCCCDDD\n>2\nAAACCCDDD\n")
	path_tf2 = tempfile.mktemp()
	callFunctionQuiet("mafft " + path_tf + " > " + path_tf2)
	with open(path_tf2, "r") as g:
		a=g.read()
	os.remove(path_tf)
	os.remove(path_tf2)
	if not a == ">1\nAAABBBCCCDDD\n>2\nAAA---CCCDDD\n":
                print(" - failed")
		return False
	else:
		print(" - ok")
		return True

def CanRunGeneric(package, packageFormatted):
	sys.stdout.write("Test can run \"" + package + "\"")
	runit = subprocess.call("type " + package, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
	if runit:
                print(" - ok")
                return True
        else:
                print(" - failed")
		print("ERROR: Cannot run " + packageFormatted)
                print("    Please check "+ packageFormatted +" is installed and that the executables are in the system path\n")
                return False
		
def CanRunAugTrain():
	sys.stdout.write("Test can run \"autoAugTrain.pl\"")
	runit = subprocess.call("type " + "autoAugTrain.pl", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
	if runit:
		print(" - ok")
		return True
	else:
		print(" - failed")
		print("    The path to the AUGUSTUS scripts folder must be exported into the path before running, i.e.\n        export PATH=$PATH:path_to_aug_scripts_folder")
		return False

def CanRunR():
	sys.stdout.write("Test can run \"R, Rscript\"")
	runit = subprocess.call("type " + "R", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
	rsunit = subprocess.call("type " + "Rscript", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
	# Check R can run
	path_r = tempfile.mktemp()
	with open(path_r, "w") as f:
		f.write("print(\"hello world\")")
	capture = subprocess.Popen("Rscript " + path_r, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = [x for x in capture.stdout]
        stderr = [x for x in capture.stderr]
	if not (len(stdout) > 0 and len(stderr) == 0):
		print(" - failed")
		print("    R was unable to run. The following error was outputted:")
		for line in stderr:
			print("    " + line)
		return False
	# Check gamlss is installed
	with open(path_r, "w") as f:
		f.write("library(\"gamlss\"); search(); print(\"blue whale\")")
	capture = subprocess.Popen("Rscript " + path_r, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	gamlFound = any(["blue whale" in i for i in capture.stdout])
	stderr = [x for x in capture.stderr]
	if not gamlFound:
		print(" - failed")
		print("    R was unable to run. Please make sure you have the \"gamlss\" package installed.")
		return False

	os.remove(path_r)
	a=commands.getstatusoutput("Rscript " + path_r)[1]
	if "WARNING: ignoring environment value of R_HOME" in a:
		print(" - failed")
		print("    R R_HOME WARNING encountered: please call \"unset R_HOME\" in bash before running OrthoFiller\n")
		return False
	else:
		print(" - ok")
		return True

def CanRunBedTools():
	""" New versions of bedtools often seem to break things. Make extra sure that everything works as it needs to.
	    This is not an exhaustive unit test.
	"""
	path_td = tempfile.mkdtemp()
	sys.stdout.write("Testing bedtools")
	try:
		path_mockGenome = path_td + "/genome.fasta"
		with open(path_mockGenome, "w") as f:
			f.write('>chr1\nCCCACGTAGCTAAGTGAATAAGTAGCCGCGCTCGACACACAGTGATGGATACGGCAGCTGGCACCAACAAGGAGCGTGGGATGTACGCTATTTTGGCGAT\n>chr2\nGGGCTCGCAACGGTTGGCCACGGCGTCTCTCATATTCGTTTAGTAGATTAACTTTCAAGATGAAAACGAGCTCCTACTAGAGAATCTCTCAGCGCACTAT\n')
		path_mockGtf1 = path_td + "/genes1.gtf"
		with open(path_mockGtf1, "w") as f:
			f.write('chr1\tfake\tCDS\t10\t39\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\nchr1\tfake\tCDS\t61\t90\t.\t+\t-\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\nchr1\tfake\tstart_codon\t10\t12\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\nchr1\tfake\tstop_codon\t88\t90\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\nchr2\tfake\tCDS\t15\t29\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\nchr2\tfake\tCDS\t46\t60\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\nchr2\tfake\tstart_codon\t58\t60\t.\t+\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\nchr2\tfake\tstop_codon\t15\t17\t.\t+\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\n')
		path_mockGtf2 = path_td + "/genes2.gtf"
		with open(path_mockGtf2, "w") as f:
			f.write('chr1\tfake\tCDS\t5\t34\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\nchr1\tfake\tCDS\t56\t85\t.\t+\t-\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\nchr1\tfake\tstart_codon\t5\t7\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\nchr1\tfake\tstop_codon\t83\t85\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\nchr2\tfake\tCDS\t10\t24\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\nchr2\tfake\tCDS\t41\t55\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\nchr2\tfake\tstart_codon\t53\t55\t.\t+\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\nchr2\tfake\tstop_codon\t10\t12\t.\t+\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\n')
		path_mockGtf3 = path_td + "/genes3.gtf"
		with open(path_mockGtf3, "w") as f:
			f.write('chr1\tfake\tCDS\t10\t30\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\nchr1\tfake\tCDS\t20\t50\t.\t+\t-\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\nchr1\tfake\tstart_codon\t10\t12\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\nchr1\tfake\tstop_codon\t48\t50\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\nchr2\tfake\tCDS\t50\t80\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\nchr2\tfake\tCDS\t70\t90\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\nchr2\tfake\tstart_codon\t88\t90\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\nchr2\tfake\tstop_codon\t50\t52\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\n')
		# Test merge
		path_merged = path_td + "/merged.gtf"
		callFunction("sort -k1,1V " + path_mockGtf3+ " | grep CDS | bedtools merge -i - -s > " + path_merged)
		with open(path_merged, "r") as f:
			if not f.read() == 'chr1\t9\t50\t+\nchr2\t49\t90\t-\n':
				raise ValueError()
		# Test intersect
		path_intersected = path_td + "/intersected.gtf"
		callFunction("bedtools intersect -wa -wb -a " + path_mockGtf1 + " -b "+ path_mockGtf2 + " > " + path_intersected)
		with open(path_intersected, "r") as f:
			if not f.read() == 'chr1\tfake\tCDS\t10\t39\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\tchr1\tfake\tCDS\t5\t34\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\nchr1\tfake\tCDS\t61\t90\t.\t+\t-\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\tchr1\tfake\tCDS\t56\t85\t.\t+\t-\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\nchr1\tfake\tCDS\t61\t90\t.\t+\t-\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\tchr1\tfake\tstop_codon\t83\t85\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\nchr1\tfake\tstart_codon\t10\t12\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\tchr1\tfake\tCDS\t5\t34\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\nchr2\tfake\tCDS\t15\t29\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\tchr2\tfake\tCDS\t10\t24\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\nchr2\tfake\tCDS\t46\t60\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\tchr2\tfake\tCDS\t41\t55\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\nchr2\tfake\tCDS\t46\t60\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\tchr2\tfake\tstart_codon\t53\t55\t.\t+\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\nchr2\tfake\tstop_codon\t15\t17\t.\t+\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\tchr2\tfake\tCDS\t10\t24\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\n':
				raise ValueError()
		# Test subtract
		path_subtracted = path_td + "/subtracted.gtf"
		callFunction("bedtools subtract -a " + path_mockGtf1 + " -b " + path_mockGtf2 + " > " + path_subtracted)
		with open(path_subtracted, "r") as f:
			if not f.read() == 'chr1\tfake\tCDS\t34\t39\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "1"; gene_name "posGene";\nchr1\tfake\tCDS\t85\t90\t.\t+\t-\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\nchr1\tfake\tstop_codon\t88\t90\t.\t+\t0\tgene_id "posGene"; transcript_id "posGene.t1"; exon_number "2"; gene_name "posGene";\nchr2\tfake\tCDS\t24\t29\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\nchr2\tfake\tCDS\t55\t60\t.\t-\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "1"; gene_name "negGene";\nchr2\tfake\tstart_codon\t58\t60\t.\t+\t0\tgene_id "negGene"; transcript_id "negGene.t1"; exon_number "2"; gene_name "negGene";\n':
				raise ValueError()
		# Test getfasta
		path_getFasta = path_td + "/getfasta.fa"
		callFunctionQuiet("infile=" +  path_mockGtf1 + "; outfile=" + path_getFasta + "; genome=" + path_mockGenome + """;
			tf=`mktemp -d`
			gffCds="$tf/gffCds"
			gffBed="$tf/gffBed"

			#Prepare the gff
			#echo "preparing gff..."
			grep -vP "^$" $infile | awk '$3=="CDS"' > $gffCds
			cut -f1-8 $gffCds > $gffBed.1
			sed -r "s/.*transcript_id[ =]\\"?([^\\";]*)\\"?;?.*/\\1/g" $gffCds > $gffBed.2
			paste $gffBed.1 $gffBed.2 | perl -ne 'chomp; @l=split; printf "$l[0]\\t%s\\t$l[4]\\t$l[8]\\t.\\t$l[6]\\n", $l[3]-1' | sort -u | sort -k1,1V -k2,2n > $gffBed
			#Negative strand
			#echo "negative strand..."
			awk '$6=="-"' $gffBed > $gffBed.neg
			bedtools getfasta -name -s -fullHeader -fi $genome -fo $gffBed.neg.tab -bed $gffBed.neg -tab
			tac $gffBed.neg.tab | awk '{a[$1]=a[$1]""$2} END {for (i in a) {print ">"i"\\n"a[i]}}' > $gffBed.neg.fa

			#Then positive strand
			#echo "positive strand..."
			awk '$6=="+"' $gffBed > $gffBed.pos
			bedtools getfasta -name -s -fullHeader -fi $genome -fo $gffBed.pos.tab -bed $gffBed.pos -tab
			cat $gffBed.pos.tab | awk '{a[$1]=a[$1]""$2} END {for (i in a) {print ">"i"\\n"a[i]}}' > $gffBed.pos.fa

			cat $gffBed.pos.fa $gffBed.neg.fa | sed -r "s/^>(.*)$/£££>\\1###/g" | sed -r \"s/$/###/g\" | tr '\\n' ' ' | sed -r "s/£££/\\n/g" | sed -r "s/### //g" | grep -v XXX | grep -v "\*[A-Z]" | grep -v "###$" | sed -r "s/###/\\n/g" | grep -vP "^$" > $outfile

			#echo $tf
			rm -r $tf
		""")
		with open(path_getFasta, "r") as f:
			if not f.read() == '>posGene.t1\nCTAAGTGAATAAGTAGCCGCGCTCGACACAGCACCAACAAGGAGCGTGGGATGTACGCTA\n>negGene.t1\nTCTTGAAAGTTAATCGAGACGCCGTGGCCA\n':
				raise ValueError()
		#######
		callFunction("rm -r " + path_td)
		print(" - ok")
		return True
	except ValueError:
		callFunction("rm -r " + path_td)
		print(" - fail")
		return False

def checkShell():
	""" Run quick checks to ensure the relevant packages are installed, and that they do
		the things we need them to do (applies in particular to bedtools and to R).
		- perl
		- awk
		- sed
		- echo
		- cut
		- sort
		- grep
		- tac/cat
		- R/Rscript
		- hmmbuild/nhmmer/makehmmerdb
		- augustus/autoAugTrain.pl
		- orthofinder
		- bedtools
		"""
	checks = []
	checks.append(CanRunMinusH("perl", "PERL"))
	checks.append(CanRunMan("sed", "sed"))
	checks.append(CanRunMan("echo", "echo"))
	checks.append(CanRunMan("cut", "cut"))
	checks.append(CanRunMan("sort", "sort"))
	checks.append(CanRunMan("grep", "grep"))
	checks.append(CanRunMan("uniq", "uniq"))
	checks.append(CanRunMan("tac", "tac"))
	checks.append(CanRunMan("cat", "cat"))
	checks.append(CanRunMan("Rscript", "Rscript"))
	checks.append(CanRunMinusH("nhmmer", "nhmmer"))
	checks.append(CanRunMinusH("hmmbuild", "hmmbuild"))
	checks.append(CanRunMinusH("makehmmerdb", "makehmmerdb"))
	checks.append(CanRunMinusH("makeblastdb", "blast+"))
	checks.append(CanRunMinusH("bedtools", "bedtools"))
	checks.append(CanRunBlank("augustus", "augustus"))
	checks.append(CanRunAugTrain())
	checks.append(CanRunMan("mktemp", "mktemp"))
	checks.append(CanRunAwk())
	checks.append(CanRunR())
	checks.append(CanRunMafft())
	checks.append(CanRunMinusH("mcl", "mcl"))
	#Check presence of orthofinder
	ortho=CanRunOrthoFinder()
	checks.append(ortho)
	if not ortho:
		print("    Cannot find the orthofinder.py file. Either\n        1) Set the orthofinder location as an environment variable using \"export ORTHOFINDER_DIR=dir\", where dir is the path to the directory containing orthofinder.py\n        2) include orthofinder.py in the same directory as OrthoFiller; or\n        3) install an orthofinder executable in your system PATH\n    Orthofinder can be downloaded from https://github.com/davidemms/OrthoFinder")
	if not all(checks):
		print("\nSome programs required to run orthoFiller are not installed or are not callable.\nPlease ensure all of the above programs are installed and in the system path.")
		sys.exit()
	if not CanRunBedTools():
		print("The version of BedTools you're using is causing problems with OrthoFiller.\nSometimes BedTools updates prevent BedTools itself from running properly.\nTry downloading bedtools v2.25.0 and run OrthoFiller again.")
		sys.exit()

########################################################
################ Function definitions ##################
########################################################

class SeqRef(object):
	def __init__(self, str_species, str_speciesNum, seqId):
		""" Basically a wrapper to hold information about each particular sequence.
		    Because dictionaries are rubbish. 
		"""
		self.species = str_species
		self.seqId = seqId
		self.uniqueId = str(str_speciesNum) + "_" +  seqId.replace("|", "-").replace(";", "-").replace(" ", "-")
	def __eq__(self, other):
		return self.uniqueId == other.uniqueId
	def __ne__(self, other):
		return not self.__eq__(other)
	def __repr__(self):
		return self.ToString()
	def ToString(self):
		return "UniqueId:%s; Species: %s; SeqId: %s" % (self.uniqueId, self.species, self.seqId)

def addSpecies(str_species, dict_speciesInfo):
	if not str_species in dict_speciesInfo:
		print(str_species)
		existingValues = [ x["number"] for x in dict_speciesInfo.values() ]
		highestVal = max(existingValues or [0])
		dict_speciesInfo[str_species] = {}
		dict_speciesInfo[str_species]["number"] = highestVal + 1


def readOrthoFinderOutput(path_orthoFinderOutputFile, path_orthoFinderSingletonsFile, dict_speciesInfo):
	"""Read CSV file into 3-tiered dictionary: orthogoup > species > sequences.
	   The last entry is an array of sequence strings.
	   The species names need a little extra attention, since some versions of orthofinder change the names
		of the input files to auto-generate a species name.
	   Output is a dictionary, keyed by unique id valued by a seqRef object, and a list of the orthogroups, and of the singletons.
	"""
	print("Reading orthogroups from " + path_orthoFinderOutputFile)
	dict_orthogroups, sequences_orthogroups = readOrthoFinderOutputIndividual(path_orthoFinderOutputFile, dict_speciesInfo)
	print("Reading singletons from " + path_orthoFinderSingletonsFile)
	dict_singletons, sequences_singletons = readOrthoFinderOutputIndividual(path_orthoFinderSingletonsFile, dict_speciesInfo)
	sequences_all = merge_dicts(sequences_orthogroups, sequences_singletons)
	return sequences_all, dict_orthogroups, dict_singletons
	

def readOrthoFinderOutputIndividual(path_orthoFinderOutputFile, dict_speciesInfo):
	dict_groups = {}
	sequences_local = {}
	with open(path_orthoFinderOutputFile) as csvfile:
		data = csv.reader(csvfile, delimiter="\t")
		# First line contains the names of the species, which orthoFinder automatically
		# generates from the fasta input names. Need to double check these names against the ones we have.
		speciesList_og		= data.next()
		speciesList_original	= list(dict_speciesInfo.keys())
		# The useful species list is just a copy of the dict_speciesInfo species, just making sure that 
		# it's in the same order as the species listed in the orthofinder output file header.
		speciesList_useful = convertOrthoFinderSpeciesNames(speciesList_og, speciesList_original)
		# Each subsequent line has orthogroup as first entry, and grouped sequence IDs
		# for each numbered column.
		for line in data:
			groupId = line[0]
			dict_groups[groupId] = []
			for i in range(1,len(speciesList_useful)):
				species = speciesList_useful[i]
				speciesNum = dict_speciesInfo[species]["number"]
				entry = re.split("[,]*", line[i])
				# Get rid of any empty entries.
				entryClean = itertools.ifilterfalse(lambda x: x=='', entry)
				for sequence in entryClean:
					seqRef = SeqRef(species, speciesNum, sequence.strip('"').strip(" "))
					dict_groups[groupId].append(seqRef.uniqueId)
					sequences_local[seqRef.uniqueId] = seqRef
	return dict_groups, sequences_local
				
def convertOrthoFinderSpeciesNames(speciesList_og, speciesList_original):
	speciesList_useful      = [""]
	for i in range(1,len(speciesList_og)):
		species_og = speciesList_og[i]
		catch = ""
		for species_or in speciesList_original:
			species_or_noxt = os.path.splitext(species_or)[0]
			species_or_fdot = species_or.split(".")[0]
			if species_og in [species_or, species_or_noxt, species_or_fdot]:
				catch = species_or
				break
		if not catch == "":
			speciesList_useful.append(catch)
			speciesList_original.remove(catch)
		else:
			print("OrthoFinder output contains one or more names that do not corrrespond to the inputted fasta names.\n Please ensure all protein fasta names are unique.\n If you ran OrthoFiller without using the --prep option, make sure you have the latest version of OrthoFinder and try again.")
			sys.exit(1)
	return speciesList_useful

def readInputLocations(path_speciesInfoFile):
	"""Read CSV file containing the locations for the sequence files, etc.
	   Each row is a species. The columns are [proteins, gffs, genome, cds].
	   The proteins string should be the same as the species element in the OrthoFinder output.
	   As such, we use the basename of this file as the species name.
	"""
	print("Loading and checking input data locations...")
	dict_speciesInfo = {}
	with open(path_speciesInfoFile) as path_locationsFile:
		# Ignore any commented lines, typically these are headers.
		data = csv.reader((row for row in path_locationsFile if not row.startswith('#')), delimiter="\t")
		for line in data:
			if not ''.join(line).strip():
				continue
			# The "species", i.e. the basename for the source protein file
			# will be the key in the dictionary.
			str_species = os.path.basename(line[0])
			addSpecies(str_species, dict_speciesInfo)
			# Then just build up the dictionary with human-readable names
			path_aa		= checkFileExists(line[0])
			path_gff 	= checkFileExists(line[1])
			path_genome	 = checkFileExists(line[2])
			path_cds	= checkFileExists(line[3])
			dict_speciesInfo[str_species]["protein"] = path_aa
			dict_speciesInfo[str_species]["gff"]	 = path_gff
			dict_speciesInfo[str_species]["genome"]  = path_genome
			dict_speciesInfo[str_species]["cds"]	 = path_cds
			checkChromosomes(path_gff, path_genome)
			checkSequences(path_gff, path_cds, path_aa)
	return dict_speciesInfo

def gffsForOrthoGroups(path_ogDir, path_orthogroups, path_singletons, dict_speciesInfo, int_cores):
	b=[]; bs=[]
	with open(path_orthogroups) as f:
		data = csv.reader(f, delimiter="\t")
		for line in data:
			b.append(line)
	with open(path_singletons) as f:
		data = csv.reader(f, delimiter="\t")
		for line in data:
			bs.append(line)
	speciesList_og		= b[0]
	speciesList_original    = list(dict_speciesInfo.keys())
	speciesList_useful	= convertOrthoFinderSpeciesNames(speciesList_og, speciesList_original)
	orthogroups = b[1:]
	singletons  = bs[1:]
	gffPool=multiprocessing.Pool(int_cores)
	for i in range(1,len(speciesList_useful)):
		str_species=speciesList_useful[i]
		a=[];
		with open(dict_speciesInfo[str_species]["gff"]) as f:
			data = csv.reader(f, delimiter="\t")
			for line in data:
				a.append(line)
		print("Extracting orthogroup and singleton gtf files for " + str_species)
		async(gffPool, gffsForGroups, args=(a, orthogroups, path_ogDir, str_species, "_orthoProtein.gtf", i))
		async(gffPool, gffsForGroups, args=(a, singletons, path_ogDir, str_species, "_singletonProtein.gtf", i))
	gffPool.close()
	gffPool.join()	
	
def gffsForGroups(list_gff, orthogroups, path_ogDir, str_species, str_outsuffix, int_speciesNum):
	#Gff entries by transcript name
	aa=[[re.sub(".*transcript_id[ =]\"([^\"]*)\".*", r'\1', x[8]), x] for x in list_gff]
	c={};
	for x in aa: c[x[0]]=[]
	for x in aa: c[x[0]].append(x[1])
	e = dict((x[0], [c[i] for i in itertools.ifilterfalse(lambda x: x=='', re.split("[ ,]*", x[int_speciesNum]))]) for x in orthogroups)
	for orthogroup in e:
		filename = path_ogDir + "/" + orthogroup+"." + str_species + str_outsuffix
		with open(filename, 'w') as mycsvfile:
			datawriter = csv.writer(mycsvfile, delimiter = '\t',quoting = csv.QUOTE_NONE, quotechar='')
			for row in list(itertools.chain.from_iterable(e[orthogroup])):
				datawriter.writerow(row + [str_species, orthogroup])

def trainAugustus(dict_speciesInfo, path_wDir, pool):
	"""trainAugustus - Trains augustus using the genomes of the input species
	"""
	path_augWDir = path_wDir + "/augustus"
	makeIfAbsent(path_augWDir)
	#Python doesn't like us to edit a dictionary while iterating over it.
	speciesList = [ x for x in dict_speciesInfo ]
	for str_species in speciesList:
		path_genome = dict_speciesInfo[str_species]["genome"]
		path_gffForTraining=dict_speciesInfo[str_species]["gffForTraining"]
		path_augustusSpecies=dict_speciesInfo[str_species]["augustusSpecies"]
		path_augSpeciesWDir = path_augWDir + "/" + str_species
		makeIfAbsent(path_augSpeciesWDir)
		if dict_speciesInfo[str_species]["needsTraining"]:
			print("training augustus on " + str_species)
			async(pool, trainAugustusIndividual, args=(path_augustusSpecies, path_genome, path_gffForTraining, path_augSpeciesWDir))

def makeGffTrainingFile(path_inputGff, path_outputGff):
	"""For AUGUSTUS to train its algorithms correctly, we need to format
	   the gff file in a certain way.
	"""
	#print("making training file " + path_outputGff + " from  " + path_inputGff + "...")
	#print("making training file " + path_outputGff + "...")	
	path_tmp=path_outputGff + ".tmp"
	path_bases=path_outputGff + ".bases"
	getBases(path_inputGff, path_bases)
	# Make sure there are no overlaps, by randomly choosing between overlapping entries, and sort the file.
	function="infile=\"" + path_inputGff + "\"; basefile=\"" + path_bases + "\"; outfile=\"" + path_tmp + "\"; " + """
		td=`mktemp -d`
		rm -f $outfile

		#echo "Assuring transcripts..."
		infile_td=`mktemp $td/infile_tid.XXXXXX`
		sed -r  '/transcript_id/! s/gene_id([ =])\\"([^\\"]*)\\";?( ?)/gene_id\\1\\"\\2\\"; transcript_id\\1\\"\\2.t999\\";\\3/g' $infile > $infile_td

		#echo "Grouping into regions.."
		sort -k1,1V -k4,4n $basefile | bedtools merge -s -i - > $td/gffmerged.bed.tmp

		cut -f1,2,3,4 $td/gffmerged.bed.tmp | sed -r "s/\\t([^\\t]*)$/\\t.\\t.\\t\\1/g" > $td/gffmerged.bed

		#echo "Intersecting..."
		bedtools intersect -a $td/gffmerged.bed -b $infile_td -wa -wb > $td/gffis.bed
		
		cat $td/gffis.bed | shuf | sed -r  "s/(.*transcript_id[ =]\\")([^\\"]*)(\\".*)/\\2\\t\\1\\2\\3\\t\\2/g" | awk 'BEGIN {FS="\\t"} {if (a[$2"."$3"."$4"."$7] == "") { a[$2"."$3"."$4"."$7]=$1 } ; if (a[$2"."$3"."$4"."$7]==$1) {v[$2"."$3"."$4"."$7]=v[$2"."$3"."$4"."$7]"\\n"$0"\\t"$2"."$3"."$4"."$7 } } END { for ( i in a ) {print v[i] } } ' | awk 'NF' | cut -f8- | sed -r "s/.\tgene_id/.\tgene_id/g" | sed -r "s/\.\-/\.neg/g" | sed -r "s/\.\+/\.pos/g" > $td/tmp1
		awk -F "\\t" '{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t.\\t"$7"\\t.\\tgene_id \\""$11".gene\\"; transcript_id \\""$11".gene.t1\\";\t"$11}' $td/tmp1 | sort -u > $outfile
		
		rm $basefile
		rm -r $td"""
	callFunction(function)
	print("check st(art|op) codon consistency")
	# Check each gene has a start codon and a stop codon and that they're in the right place
	checkCdsHealth(path_tmp, path_outputGff)
	os.remove(path_tmp)
	# Add exons as well as CDS
	function="infile=\"" + path_outputGff + "\"; tmpfile=`mktemp`; tmpfile2=`mktemp`; grep -P \"\\tCDS\\t\" $infile | sed -r \"s/\\tCDS\\t/\\texon\\t/g\" | sed -r \"s/\\t[^\\t]*\\tgene_id/\\t\\.\\tgene_id/g\" > $tmpfile2; cat $infile $tmpfile2 | sort -u | sort -k1,1V -k4,4n > $tmpfile; mv $tmpfile $infile; rm $tmpfile2"
	callFunction(function)

def getBases(path_gtf, path_gtfBases):
	"""Get the base for a gtf file
	"""
	with open(path_gtf) as c:
		gtf = [line for line in csv.reader(c, delimiter="\t") if ''.join(line).strip()]
	coords={}
	entries={}
	for line in gtf:
		if not line[2].lower() == "cds":
			continue
		transcript_id = re.sub(r'.*transcript_id \"([^\"]*)\".*', r'\1', line[8])
		if not transcript_id in coords:
			coords[transcript_id] = []
			entries[transcript_id] = line
		coords[transcript_id] += [int(line[3]), int(line[4])]
	for t_id in entries:
		entries[t_id][3] = min(coords[t_id])
		entries[t_id][4] = max(coords[t_id])
	with open(path_gtfBases, "w") as p:
		writer = csv.writer(p, delimiter="\t", quoting = csv.QUOTE_NONE, quotechar='')
		for entry in entries:
			writer.writerow(entries[entry])

def checkCdsHealth(path_inputGtf, path_outputGtf):
	with open(path_inputGtf) as c:
		gtf = [line for line in csv.reader(c, delimiter="\t") if ''.join(line).strip()]
	transcripts={}; cds={}; starts={}; stops={}; entries={}; strands={}
	for line in gtf:
		transcript_id = re.sub(r'.*transcript_id \"([^\"]*)\".*', r'\1', line[8])
		if not transcript_id in transcripts:
			transcripts[transcript_id] = []
			entries[transcript_id] = []
		if not transcript_id in strands:
			strands[transcript_id] = []
		strands[transcript_id] += [line[6]]
		if line[2].lower() == "start_codon":
			if not transcript_id in starts:
				starts[transcript_id] = []
			starts[transcript_id] += [[int(line[3]), int(line[4])+1]]
		elif line[2].lower() == "stop_codon":
			if not transcript_id in stops:
				stops[transcript_id] = []
			stops[transcript_id] += [[int(line[3]), int(line[4])+1]]
		elif line[2].lower() == "cds":
			if not transcript_id in cds:
				cds[transcript_id] = []
			cds[transcript_id] += range(int(line[3]), int(line[4])+1)
		entries[transcript_id].append(line)
	badGenes=[]
	for t_id in transcripts:
		if len(cds[t_id]) != len(set(cds[t_id])):
			badGenes.append(t_id); continue
		if len(set(strands[t_id])) != 1:
			badGenes.append(t_id); continue
		if not t_id in stops or not t_id in starts or len(stops[t_id]) !=1 or len(starts[t_id]) !=1:
			badGenes.append(t_id); continue
		if strands[t_id][0] == "+":
			startstart=starts[t_id][0][0]
			cdsstart=min(cds[t_id])
			if cdsstart != startstart:
				badGenes.append(t_id); continue
			endend=stops[t_id][0][1]
			cdsend=max(cds[t_id])
			if cdsend != endend -1 :
				badGenes.append(t_id); continue
		elif strands[t_id][0] == "-":
			startstart=starts[t_id][0][1]
			cdsstart=max(cds[t_id])
			if cdsstart != startstart -1:
				badGenes.append(t_id); continue
			endend=stops[t_id][0][0]
			cdsend=min(cds[t_id])
			if cdsend != endend:
				badGenes.append(t_id); continue
		else:
			badGenes.append(t_id); continue
	with open(path_outputGtf, "w") as p:
		writer = csv.writer(p, delimiter="\t", quoting = csv.QUOTE_NONE, quotechar='')
		for t_id in entries:
			if not t_id in badGenes:
				for entry in entries[t_id]:
					writer.writerow(entry[0:9])

def trainAugustusIndividual(str_augustusSpecies, path_genome, path_gff, path_augSpeciesWDir):
	callFunctionQuiet("autoAugTrain.pl --useexisting -v -v -v --workingdir=" + \
		path_augSpeciesWDir + " --species=" + str_augustusSpecies + \
		" --trainingset=" + path_gff + " --genome=" + path_genome)

def getProteinSequences(sequencesHolder, dict_speciesInfoDict):
	"""Get the protein sequences for each feature per species in an orthogroup
	   dict_speciesInfo should be the dictionary version of the file locations, indexed by species
	"""
	# To make protein sequence access faster, index.
	path_indexedProteinFiles = {}
	for str_species in dict_speciesInfoDict:
		path_indexedProteinFiles[str_species] = SeqIO.index(dict_speciesInfoDict[str_species]["protein"], "fasta")
	dict_proteinSequencesHolder = {}
	for sequence in sequencesHolder:
		seqRef = sequencesHolder[sequence]
		proteinSequence = path_indexedProteinFiles[seqRef.species][seqRef.seqId]
		# Replace the species-local id with the uniqueId
		proteinSequence.id = seqRef.uniqueId
		dict_proteinSequencesHolder[seqRef.uniqueId] = proteinSequence
	# Returns a dictionary of id / SeqIO sequence object.
	return dict_proteinSequencesHolder

def writeSequencesToFastaFile(dict_proteinSequencesHolder, path_outputFile):
	"""Write out sequences into fasta file. Overwrites file if necessary.
	"""
	actualSequences = dict_proteinSequencesHolder.values()
	SeqIO.write(actualSequences, path_outputFile, "fasta")

def makeProteinAlignment(path_proteinFastaFile, path_fastaOut):
	"""Makes an alignment and save it in fasta format to the fastaOut.
	"""
	# "auto" means l-ins-i is used when the protein set is small enough, FFT-NS2 otherwise.
	alignment = MafftCommandline(input=path_proteinFastaFile, auto="on")
	stdout, stderr = alignment()
	with open(path_fastaOut, "w") as outHandle:
		outHandle.write(stdout)

def getAlignmentStats(path_proteinAlignmentFastaFile):
	"""Calculates some statistics on the alignment, specifically gap quantities.
	"""
	protAl = SeqIO.parse(path_proteinAlignmentFastaFile, "fasta")
	# We're going to write the sequences into a matrix of characters.
	protAlSequencesCharArray = []
	for sequence in protAl:
		protAlSequencesCharArray.append(list(str(sequence.seq)))
	# The sequences should all be the same length
	# Count the number of gaps in each sequence and output it as a percentage of
	# The total length of the sequence.
	sequenceGapCounts = list(x.count('-') for x in protAlSequencesCharArray)
	sequenceLengths = list(len(x) for x in protAlSequencesCharArray)
	firstSequenceLength = len(protAlSequencesCharArray[0])
	if not sequenceLengths.count(firstSequenceLength) == len(protAlSequencesCharArray):
		raise ValueError('There is a problem with the protein alignment, lengths should be equal.')
	sequenceGapPercentages = []
	for value in sequenceGapCounts:
		sequenceGapPercentages.append(value / float(firstSequenceLength))
	# Now turn the array around and find the gap percentage per position.
	protAlSequencesCharArrayTransposed = zip(*protAlSequencesCharArray)
	positionGapCounts = list(x.count('-') for x in protAlSequencesCharArrayTransposed)
	positionLength = len(protAlSequencesCharArrayTransposed[0])
	positionGapPercentages = []
	for value in positionGapCounts:
		positionGapPercentages.append(value / float(positionLength))
	# Output
	return sequenceGapPercentages, positionGapPercentages
	
def threadGappedProteinSequenceThroughDNA(gappedProteinSequence, dnaSourceSequence):
	"""Protein must correspond exactly to DNA sequence.
	   Both sequences should be provided as strings.
	"""
	gappedDnaSequence = ""
	int_counter = 0
	for ch in gappedProteinSequence:
		if ch =="-":
			gappedDnaSequence = gappedDnaSequence + "---"
		else:
			codon=dnaSourceSequence[int_counter:int_counter+3]
			if not len(codon) == 3:	
				print "the length of the codon is " + str(len(codon))
				print "the codon is " + codon
				print "I got through the checks"
				if not (ch == "X" or ch == "x"):
					print "Slight mismatch, be careful"
				codon = codon + "-" * (-len(codon) % 3)
				print codon
			gappedDnaSequence = gappedDnaSequence + codon
			int_counter = int_counter + 3
	return gappedDnaSequence
	# CC - stop codons don't get included - is this important?

def getNucleotideAlignment(alignedProteinsFastaIn, fastaOut, sequencesHolder, dict_speciesInfoDict):
	"""Get a nucleotide alignment from the protein alignment.
	"""
	# CC - will want to have a global index at some point.
	path_indexedCdsFiles = {}
	for str_species in dict_speciesInfoDict:
		path_indexedCdsFiles[str_species] = SeqIO.index(dict_speciesInfoDict[str_species]["cds"], "fasta")
	with open(alignedProteinsFastaIn, "rU") as inputHandle:
		fastaRecords = SeqIO.parse(inputHandle, "fasta")
		nucleotideAlignments = []
		for record in fastaRecords:
			# The recordId is the absolute id made previously
			recordId = record.id
			# We need to use species name and local Id to look up in the cds file.
			str_species = sequencesHolder[recordId].species
			localId = sequencesHolder[recordId].seqId
			# Construct the new sequences
			cdsRecord = path_indexedCdsFiles[str_species][localId]
			# The ungapped dna sequence
			dnaSourceSequence = str(cdsRecord.seq)
			# The gapped protein sequence
			gappedProteinSequence = record.seq
			# Get the gapped dna sequence
			gDnaSeq = threadGappedProteinSequenceThroughDNA(gappedProteinSequence, dnaSourceSequence)
			# Construct the sequence object and put it out there.
			gDnaId = cdsRecord.id
			gDnaName = cdsRecord.name
			gDnaDesc = cdsRecord.description
			gDnaSeqOut = SeqRecord(Bio.Seq.Seq(gDnaSeq), id=gDnaId, name=gDnaName, description=gDnaDesc)
			nucleotideAlignments.append(gDnaSeqOut)
	# Write the sequences to the specified output file.
	with open(fastaOut, "w") as outputHandle:
		SeqIO.write(nucleotideAlignments, outputHandle, "fasta")

def buildHmm(nucAlignment, path_outputFile):
	"""Build an hmm based on a nucleotide alignment. Inputs are file names.
	"""
	callFunctionQuiet("hmmbuild --informat afa " + path_outputFile + " " + nucAlignment) #qgr

def makeHmmerDb(path_genomeFile, path_dbOutput):
	"""Makes a database per cds file for use with hmmer.
	"""
	callFunctionQuiet("makehmmerdb --block_size=10 " + path_genomeFile + " " + path_dbOutput) #qgr

def implementHmmSearch(path_hmmFile, path_db, path_hitsFileName):
	"""Runs across the genome and finds hmm hits
	"""
	callFunctionQuiet("nhmmer --tformat hmmerfm --dna --cpu 1 --tblout " + path_hitsFileName + " " + 	path_hmmFile + " " + path_db) #qgr

def processOg(orthogroup, list_orthogroupSequenceIds, orthogroupProteinSequences, dict_sequenceInfoById, dict_speciesInfo, path_wDir, str_ogFolder):
	"""Runs all alignments and markov models for a particular orthogroup.
	"""
	########################################################
	# Define output files
	########################################################
	path_protSeqFile = path_wDir + "/" +  str_ogFolder + "/" + orthogroup + "_ProteinSequences.fasta"
	path_proteinAlignmentFile = path_wDir + "/" + str_ogFolder + "/" + orthogroup + "_ProteinAlignment.fasta"
	path_nucAlignmentFile = path_wDir + "/" + str_ogFolder + "/" + orthogroup + "_NucAlignment.fasta"
	#######################################################
	# Write out the sequences and make the alignments
	#######################################################
	writeSequencesToFastaFile(orthogroupProteinSequences, path_protSeqFile)
	makeProteinAlignment(path_protSeqFile, path_proteinAlignmentFile)
	getNucleotideAlignment(path_proteinAlignmentFile, path_nucAlignmentFile, dict_sequenceInfoById, dict_speciesInfo)
	#######################################################
	# Debug - get some stats on the alignments
	#######################################################
	"""
	path_pAlnStatsFile = path_wDir + "/" +  str_ogFolder + "/" + orthogroup + "_ProteinAlignment.stats.txt"
	sequenceGapStats, positionGapStats = getAlignmentStats(path_proteinAlignmentFile)
	with open(path_pAlnStatsFile, "w") as p:
		p.write(re.sub('[\[\]]', '', str(sequenceGapStats,)))
		p.write("\n")
	p.write(re.sub('[\[\]]', '', str(positionGapStats,)))
	"""
	######################################################
	# Build the HMM.
	######################################################
	path_hmmFile = path_wDir + "/" + str_ogFolder + "/" + orthogroup + ".hmm"
	buildHmm(path_nucAlignmentFile, path_hmmFile)
	######################################################
	# Search the genome of each species in turn.
	######################################################
	for species in dict_speciesInfo:
		path_hitsFile = path_wDir + "/" + str_ogFolder + "/" + orthogroup + "." + species + ".hits"
		#print "Generating hits file: " + path_hitsFile
		implementHmmSearch(path_hmmFile, dict_speciesInfo[species]["hmmdb"], path_hitsFile)
		#Form a bed file from the resultant hits file
		path_hitsFileBed = path_hitsFile + ".bed"
		callFunction("grep -v \"#\" " + path_hitsFile + " | sed -r \"s/ +/\t/g\" | cut -f1,7,8,12,13,14,15 | perl -ne \
				'chomp;@l=split; printf \"%s\\t%s\\t%s\\t.\\t.\\t%s\\t" + species + "\\t" + orthogroup +
				"\\n\", $l[0], ($l[1] + $l[2] - abs($l[1] - $l[2])) / 2, ($l[1] + $l[2] + abs($l[2] - $l[1])) / 2,\
				 join(\"\\t\", @l[3..6])' > " + path_hitsFileBed)
	#print("Finished " + orthogroup)

def proposeNewGenes(path_hitsOgIntersectionFileNameAnnotated, path_allHitsOgIntersectionFileNameAnnotated, str_speciesName, path_candidatesFile, hitFilter):
	# Use R to find candidates.
	# Cupcakes, doesn't work if the total hit count is less than 1000.
	fitDistributions(path_hitsOgIntersectionFileNameAnnotated, \
				path_allHitsOgIntersectionFileNameAnnotated, \
				path_candidatesFile, \
				hitFilter)

def fitDistributions(path_hitsOgIntersectionFileNameAnnotated, path_allHitsOgIntersectionFileNameAnnotated, path_candidatesFile, hitFilter):
	"""Use R to fit distributions to the hit score data
	"""
	rfile=tempfile.mkstemp(suffix=".R")[1]
	if hitFilter:
		unpackFitDistributionScript(rfile)
	else:
		unpackFitDistributionScript_noFilter(rfile)
	callFunctionQuiet("rfile=\""+rfile+"\"; Rscript $rfile " + path_hitsOgIntersectionFileNameAnnotated + " " + \
				path_allHitsOgIntersectionFileNameAnnotated + " "  + path_candidatesFile + "; rm $rfile")

def prepareOutputFolder(path_outDir):
	if path_outDir == "":
		str_prefix = "OrthoFillerOut_" + datetime.datetime.now().strftime("%y%m%d") + "_RunId_"
		path_outDir = tempfile.mkdtemp(prefix=str_prefix, dir=".")
	path_wDir = path_outDir + "/working"
	path_ogDir = path_wDir + "/orthogroups"
	path_resultsDir = path_outDir + "/results"
	makeIfAbsent(path_outDir)
	makeIfAbsent(path_wDir)
	makeIfAbsent(path_ogDir)
	makeIfAbsent(path_resultsDir)
	return path_resultsDir, path_wDir

def prepareHmmDbs(dict_speciesInfo, path_wDir, int_cores):
	hmmdbpool=multiprocessing.Pool(int_cores)
	for str_speciesName in dict_speciesInfo:
		print("Preparing database for " + str_speciesName)
		if not "hmmdb" in dict_speciesInfo[str_speciesName]:
			path_genomeFile = dict_speciesInfo[str_speciesName]["genome"]
			path_db = path_wDir + "/" + str_speciesName + ".hmmdb"
			dict_speciesInfo[str_speciesName]["hmmdb"] = path_db
			async(hmmdbpool, makeHmmerDb,  args=(path_genomeFile, path_db))	
	hmmdbpool.close()
	hmmdbpool.join()

def run(dict_speciesInfo, dict_sequenceInfoById, orthogroups, singletons, path_resultsDir, path_wDir, path_orthoFinderOutputFile, path_singletonsFile, int_cores=16, firstPass=False, augOnly=False, hitFilter=True, hintFilter=True):
	"""Takes orthofinder output and a collection of genome info locations as input.
	   Processes orthogroups in parallel, filters hits, and generates gene models.
	"""
	######################################################
	# If we're running on first pass mode, we want to keep
	# all of the augustus-independent files for later, but
	# don't want the first pass stuff to interfere with
	# subsequent runs. So create a separate folder.
	######################################################
	if firstPass:
		path_wDirS=path_wDir + "/firstpass_working"
		makeIfAbsent(path_wDirS)
	else:
		path_wDirS=path_wDir
	#####################################################
	# Set off the Augustus training
	#####################################################
	trainingPool = multiprocessing.Pool(int_cores)
	print("Training AUGUSTUS")
	trainAugustus(dict_speciesInfo, path_wDir, trainingPool)
	#####################################################
	# If we're on a second pass, where the first pass was
	# used to create training files, we don't need to 
	# calculate all the hmms and hits.
	#####################################################
	if not augOnly:
		######################################################
		# Set up an hmm database for each species
		######################################################
		print("\n2.1. Preparing HMM databases")
	        print(  "============================")
		prepareHmmDbs(dict_speciesInfo, path_wDir, int_cores)#ql
		#####################################################
		# Produce gff files for each orthogroup/species pair
		#####################################################		
		str_ogDir = "orthogroups"	
		path_ogDir = path_wDir + "/" + str_ogDir
	        print("\n2.2. Extracting orthogroup gtf files")
                print(  "====================================")
		gffsForOrthoGroups(path_ogDir, path_orthoFinderOutputFile, path_singletonsFile, dict_speciesInfo, int_cores)#ql
		#####################################################
		# Process each individual orthogroup in parallel
		#####################################################
		proteinSequences = getProteinSequences(dict_sequenceInfoById, dict_speciesInfo)
		og_pool = multiprocessing.Pool(int_cores)
		jobs={}
		int_counter = 1
		str_total = str(len(orthogroups))
		print("\n2.2. Processing orthogroup HMMs")
                print(  "==============================================")
		for orthogroup in orthogroups:
			#print "Submitting " + orthogroup + "; " + str(int_counter) + " of " + str_total + " submitted."
			orthogroupProteinSequences = dict((x, proteinSequences[x]) for x in orthogroups[orthogroup] )
			#processOg(orthogroup, orthogroups[orthogroup], orthogroupProteinSequences, dict_sequenceInfoById, dict_speciesInfo, path_wDir, str_ogDir)
			jobs[orthogroup] = async(og_pool, processOg, args=(orthogroup, \
							orthogroups[orthogroup], \
							orthogroupProteinSequences, \
							dict_sequenceInfoById, \
							dict_speciesInfo, \
							path_wDir, \
							str_ogDir))#qr
			int_counter = int_counter + 1
		og_pool.close()
		# Keep track of how many orthogroups have been processed.
		while True:
			k=[jobs[j].ready() for j in jobs]
			sys.stdout.write("\r" + str(sum(k)) + " out of " + str_total +  " orthogroups processed")
			sys.stdout.flush()
			if all(k):
				print("\nFinished processing orthogroups")
				break
			#print([k for k in jobs if not jobs[k].ready()])
			time.sleep(2)
		og_pool.join()
		####################################################
		# Start a new pool for processing the hmm outfiles.
		####################################################
		pool = multiprocessing.Pool(int_cores)
		dict_ogIntersectionFileNamesAnnotated = {}
		print("\n3. Processing HMM output files")
		print(  "==============================")
		for str_speciesName in dict_speciesInfo:
			print("Submitting HMM output files for species " + str_speciesName + "...")
			# Prepare file names
			path_hitsBedFileName = path_wDir + "/" + str_speciesName + ".allHits.bed"
			path_ogBedFileName = path_wDir + "/" + str_speciesName + ".allOrthogroups.bed"
			path_hitsOgIntersectionFileName = path_wDir + "/" + str_speciesName + ".hitsIntersectOrthogroups.bed"
			# Get all hits into one file
			callFunction("find  " + path_ogDir + " -name \"OG*" + str_speciesName + "*hits.bed\" | xargs -n 32 cat | sed -r \"s/gene_id=*[^\\\"]*\\\"/gene_id=\\\"/g\" | sort -k1,1 -k2,2n | awk '$2 >0 && $3 > 0'  | sort -k1,1 -k2,2n > " + path_hitsBedFileName)#ql
			# Get all orthos and singletons into one file
			# Bedtools gets upset if we try to intersect with an empty file, so as a hack also provide a fake
			# entry in the same format. Hope that this never pops up in real life.
			callFunction("(find  " + path_ogDir + " -name \"OG*" + str_speciesName + "*Protein.gtf\" | xargs -n 32 cat | grep -v \"inary\" | awk '$3==\"CDS\"' | sed -r \"s/ \\\"/\\\"/g\" | sed -r \"s/; /;/g\" | cut -f1,4,5,6,7,9,10,11 | perl -ne 'chomp;@l=split; printf \"$l[0]\\t%s\\t%s\\t.\\t%s\\t\\n\", $l[1]-1, $l[2]-1, join(\"\\t\", @l[3..7])' > " + path_ogBedFileName + "; printf \"chr_FAKE_QKlWlKgGS4\\t0\\t1\\t.\\t.\\t-\\tgene_id=\\\"FAKE\\\"\\tfake.fasta\\tOG9999999\") |  sort -k1,1 -k2,2n >> " + path_ogBedFileName)
			#Now intersect
			callFunction("bedtools intersect -nonamecheck -a " + path_hitsBedFileName + " -b " + path_ogBedFileName + " -wa -wb > " + path_hitsOgIntersectionFileName)#ql
			callFunction("cat " + path_hitsOgIntersectionFileName + " " + path_hitsBedFileName + " | cut -f1-11 | sort | uniq -u | sed -r \"s/$/\\t.\\t.\\t.\\t.\\t.\\t.\\t.\\t.\\t./g\" > " + path_hitsOgIntersectionFileName + ".tmp; cat " + path_hitsOgIntersectionFileName + ".tmp " + path_hitsOgIntersectionFileName + " > " + path_hitsOgIntersectionFileName + ".tmp.tmp ; mv " + path_hitsOgIntersectionFileName + ".tmp.tmp " + path_hitsOgIntersectionFileName + "; rm " + path_hitsOgIntersectionFileName + ".tmp")#ql
			path_hitsOgIntersectionFileNameAnnotated = path_wDir + "/" + str_speciesName + ".hitsIntersectionOrthogroups.annotated.bed"
			#Annotate whether each line is a good match, a bad match, or a candidate match.
			#We don't need to distinguish singletons and orthos.
			async(pool, annotateIntersectedOutput, args=(path_hitsOgIntersectionFileName, path_hitsOgIntersectionFileNameAnnotated))#ql
			dict_ogIntersectionFileNamesAnnotated[str_speciesName] = path_hitsOgIntersectionFileNameAnnotated
		pool.close()
		pool.join()
		print("Done processing HMM output files")
		####################################################
		# Concatenate all the files in case we need to 
		# use the aggregate distribution.
		####################################################
		print("Generating concatenated version of HMM output")
		path_allHitsOgIntersectionFileNameAnnotated = path_wDir + "/allSpecies.hitsIntersectionOrthogroups.annotated.bed"
		deleteIfPresent(path_allHitsOgIntersectionFileNameAnnotated)#ql
		for str_speciesName in dict_ogIntersectionFileNamesAnnotated:
			path_annotatedFile = dict_ogIntersectionFileNamesAnnotated[str_speciesName]
			appendFileToFile(path_annotatedFile, path_allHitsOgIntersectionFileNameAnnotated)
		print("Checking quantity of reference hits...")
		with open(path_allHitsOgIntersectionFileNameAnnotated, "r") as f:
			data = list(csv.reader(f, delimiter="\t"))
			match_good = [a for a in data if a[20] == "match_good"]
			match_bad  = [a for a in data if a[20] == "match_bad"]
			if match_good < 1000 or match_bad < 1000:
				print("Error: hit quantity is not sufficient to classify potential new genes.")
				sys.exit()			
		####################################################
		# Fit a model for each individual species. If data
		# is insufficient, use aggregated data.
		####################################################
		pool = multiprocessing.Pool(int_cores)
		for str_speciesName in dict_speciesInfo:
			path_outFile = path_wDir + "/" + str_speciesName + ".proposedGenes"
			dict_speciesInfo[str_speciesName]["proposedgenes"] = path_outFile
			path_hitsOgIntersectionFileNameAnnotated = dict_ogIntersectionFileNamesAnnotated[str_speciesName]
			async(pool, proposeNewGenes, args=(path_hitsOgIntersectionFileNameAnnotated, path_allHitsOgIntersectionFileNameAnnotated, str_speciesName, path_outFile, hitFilter))#ql
		pool.close()
		pool.join()
#qx	sys.exit(1)
	####################################################
	# Run Augustus. We need the training pool to have 
	# finished by this point. Parse output
	####################################################
	print("Waiting for training to finish before continuing...")
	trainingPool.close()
	trainingPool.join()
	print("Done training.")
        print("\n4. Running Augustus")
        print(  "==============================")
	augustusPool = multiprocessing.Pool(int_cores)
	for str_speciesName in dict_speciesInfo:
		path_proposedGenes = dict_speciesInfo[str_speciesName]["proposedgenes"]
		path_genome = dict_speciesInfo[str_speciesName]["genome"]
		path_sourcegff=dict_speciesInfo[str_speciesName]["gff"]
		str_augustusOutNameStub= path_wDirS + "/" + str_speciesName + ".proposedGenes"
		path_augustusOut = str_augustusOutNameStub + ".AugustusModels.gff"
		path_fastaOut = str_augustusOutNameStub + ".AugustusParsed.sequences.fasta"
		path_augustusParsedOut = str_augustusOutNameStub + ".AugustusParsed.gff"	
		dict_speciesInfo[str_speciesName]["augustusoutput"] = path_augustusOut
		dict_speciesInfo[str_speciesName]["augustusparsed"] = path_augustusParsedOut
		dict_speciesInfo[str_speciesName]["augustussequences"] = path_fastaOut
		path_hintsFile = str_augustusOutNameStub + ".hints.gff"
		dict_speciesInfo[str_speciesName]["hints"] = path_hintsFile
		print("Running Augustus on " + str_speciesName)
		if not dict_speciesInfo[str_speciesName]["indirectAugustus"]:
			path_augustusSpeciesName = dict_speciesInfo[str_speciesName]["augustusSpecies"]
			async(augustusPool, runAndParseAugustus, args=(path_proposedGenes, path_genome, path_augustusOut, path_augustusParsedOut, path_fastaOut, path_augustusSpeciesName, path_hintsFile, path_sourcegff))
		else:
			path_otherSpeciesResults = path_wDirS + "/" + str_speciesName + ".augustus_otherSpecies"
			makeIfAbsent(path_otherSpeciesResults)
			dict_speciesInfo[str_speciesName]["indirectAugustusFolder"]=path_otherSpeciesResults
			otherSpecies = [ x for x in dict_speciesInfo if not dict_speciesInfo[x]["indirectAugustus"]]
			for str_otherSpecies in otherSpecies:
				print(otherSpecies)
				path_otherSpeciesAugustusOut = path_otherSpeciesResults + "/" + str_speciesName + ".proposedGenes." + str_otherSpecies + ".AugustusModels.gff"
				path_otherSpeciesAugustusParsedOut =  path_otherSpeciesResults + "/" + str_speciesName + ".proposedGenes." + str_otherSpecies + ".AugustusParsed.gff"
				path_otherSpeciesFastaOut =  path_otherSpeciesResults + "/" + str_speciesName + ".proposedGenes." + str_otherSpecies + ".AugustusParsed.sequences.fasta"
				otherSpeciesAugustusSpeciesName = dict_speciesInfo[str_otherSpecies]["augustusSpecies"]
				async(augustusPool, runAndParseAugustus, args=(path_proposedGenes, path_genome, path_otherSpeciesAugustusOut, path_otherSpeciesAugustusParsedOut, path_otherSpeciesFastaOut, otherSpeciesAugustusSpeciesName, path_hintsFile, path_sourcegff))
				print(otherSpeciesAugustusSpeciesName)
	augustusPool.close()
	augustusPool.join()
	####################################################
	# Combine data from the species on which augustus 
	# has been run indirectly
	####################################################
	combinePool = multiprocessing.Pool(int_cores)
	for str_speciesName in [ x for x in dict_speciesInfo if dict_speciesInfo[x]["indirectAugustus"]]:
		path_augustusParsedOut = dict_speciesInfo[str_speciesName]["augustusparsed"]
		path_fastaOut = dict_speciesInfo[str_speciesName]["augustussequences"]
		path_otherSpeciesResults = dict_speciesInfo[str_speciesName]["indirectAugustusFolder"]
		async(combinePool, combineIndirectAugustusResults, args=(path_otherSpeciesResults, path_augustusParsedOut, path_fastaOut))
	combinePool.close()
	combinePool.join()
	####################################################
	# Get a hint f score for the new genes and abandon
	# those genes whose score is not adequate.
	####################################################
	if hintFilter:
	        print("\n5.1 Filtering by hint F-score")
	        print(  "=============================")
	pool=multiprocessing.Pool(int_cores)
	for str_speciesName in dict_speciesInfo:
		path_augustusParsed=dict_speciesInfo[str_speciesName]["augustusparsed"]
		path_hintFile=dict_speciesInfo[str_speciesName]["hints"]
		path_augustusParsedHintFiltered=path_augustusParsed+".hintfiltered.gff"
		dict_speciesInfo[str_speciesName]["augustusparsed_hintfiltered"]=path_augustusParsedHintFiltered
		num_threshold=0.8
		path_augustusSequences=dict_speciesInfo[str_speciesName]["augustussequences"]	
		path_augustusSequencesHintFiltered=path_augustusSequences + ".hintfiltered.fasta"
		dict_speciesInfo[str_speciesName]["augustussequences_hintfiltered"]=path_augustusSequencesHintFiltered	
		if hintFilter:
			async(pool, hintFscoreFilter, args=(path_augustusParsed, path_hintFile, path_augustusParsedHintFiltered, num_threshold,  path_augustusSequences, path_augustusSequencesHintFiltered))
		else:
			dict_speciesInfo[str_speciesName]["augustusparsed_hintfiltered"]=path_augustusParsed
			dict_speciesInfo[str_speciesName]["augustussequences_hintfiltered"]=path_augustusSequences
	pool.close()
	pool.join()
	####################################################
	# Reinsert the sequences into the proteome and 
	# rerun OrthoFiller
	####################################################
	print("\n5.2 Re-running OrthoFinder")
        print(  "==========================")
	path_newProteomesDir = path_wDirS + "/newProteomes"
	deleteIfPresent(path_newProteomesDir)
	makeIfAbsent(path_newProteomesDir)
	dict_speciesInfo_modern = {}
	for str_speciesName in dict_speciesInfo:
		path_newProteome = path_newProteomesDir + "/" + str_speciesName + "newProteome.fasta"
		path_oldProteome = dict_speciesInfo[str_speciesName]["protein"]
		dict_speciesInfo[str_speciesName]["newProtein"] = path_newProteome
		dict_speciesInfo[str_speciesName]["newSpeciesName"] =  str_speciesName + "newProteome.fasta"
		path_predictedProteinSequences = dict_speciesInfo[str_speciesName]["augustussequences_hintfiltered"]
		makeNewProteome(path_oldProteome, path_predictedProteinSequences, path_newProteome)
		dict_speciesInfo_modern[str_speciesName + "newProteome.fasta"] = {}
		dict_speciesInfo_modern[str_speciesName + "newProteome.fasta"]["number"] = dict_speciesInfo[str_speciesName]["number"]
	runOrthoFinder(path_newProteomesDir, int_cores)
	####################################################
	# Check genes have ended up in the right orthogroup
	####################################################
	# Fetch the original species set before we add the new ones.
	print("\n5.3 Running orthogroup membership test")
        print(  "======================================")
	path_orthofinderOutputNew	= getOrthogroupsFile(path_newProteomesDir)
	path_orthofinderSingletonsNew	= getSingletonsFile(path_newProteomesDir)
	silNew, oNew, sNew = readOrthoFinderOutput(path_orthofinderOutputNew, \
									path_orthofinderSingletonsNew, dict_speciesInfo_modern)
	for str_speciesName in dict_speciesInfo:
		print("Double-checking membership for species " + str_speciesName)
		#Prepare file names.
		path_augustusParsed = dict_speciesInfo[str_speciesName]["augustusparsed_hintfiltered"]
		path_acceptedSequencesOut = path_resultsDir + "/" + str_speciesName + ".newSequences.fasta"
		path_augustusParsedUniq = path_augustusParsed + ".uniq"
		dict_speciesInfo[str_speciesName]["acceptedsequences"] = path_acceptedSequencesOut
		callFunction("grep transcript_id " + path_augustusParsed + " | sed -r 's/.*transcript_id/transcript_id/g' | sort -u > " + path_augustusParsedUniq)
		str_newSpeciesName=dict_speciesInfo[str_speciesName]["newSpeciesName"]
		print(str_newSpeciesName)
		acceptedSequences=[]
		potentialSequences={}
		with open(path_augustusParsedUniq, "r") as csvfile:
			data = csv.reader(csvfile, delimiter="\t")
			for entry in data:
				sequenceName=re.sub(";g", "; g", re.sub("_id=", "_id ", entry[0]))
				potentialSequences[sequenceName] = re.split("[, ]+", re.sub("possibleOrthos=", "", entry[1]))
		for seqId in potentialSequences:
			print(seqId)
			seqIdAlt=seqId.replace("_id ", "_id=").replace(" ", "")
			print(seqIdAlt)
			#Find out which orthogroup the sequence has been placed in.
			uniqueId=[x for x in silNew if (silNew[x].species==str_newSpeciesName and compareOutputSequences(silNew[x].seqId, seqId))][0]
			list_newOrthogroup = [x for x in oNew if uniqueId in oNew[x]]
			if not list_newOrthogroup:
				#If there is no orthogroup corresponding to this sequence, move on.
				print("There is no orthogroup containing this gene. This is an unsuccessful placement")
				continue
			newOrthogroup=list_newOrthogroup[0]
			oldOrthogroupsPotential=potentialSequences[seqId]
			success = 0
			for oldOrthogroup in oldOrthogroupsPotential:
				overlap = 0
				#Check each potential old orthogroup in turn to see
				#how much overlap it has with the new orthogroup
				for oldOrthogroupSequenceId in orthogroups[oldOrthogroup]:
					oldOrthSeq = dict_sequenceInfoById[oldOrthogroupSequenceId]
					correspondingSpecies = dict_speciesInfo[oldOrthSeq.species]["newSpeciesName"]
					for newOrthogroupSequenceId in oNew[newOrthogroup]:
						newOrthSeq = silNew[newOrthogroupSequenceId]
						if newOrthSeq.species == correspondingSpecies and compareOutputSequences(newOrthSeq.seqId, oldOrthSeq.seqId):
							overlap = overlap + 1
				if overlap > 0:
					success = 1
					break
			print("There are " + str(len(orthogroups[oldOrthogroup])) + " sequences in the trial old orthogroup. There are " + str(len(oNew[newOrthogroup])) 	+ " sequences in the new orthogroup. The overlap is " + str(overlap))
			if success == 1:
				print("This is a successful placement")
				acceptedSequences.append(seqId.replace(" ", "").replace("\"", ""))
			else:
				print("This is an unsuccessful placement")
		path_newProteome = path_newProteomesDir + "/" + str_speciesName + "newProteome.fasta"
		sequences=SeqIO.parse(path_newProteome, "fasta")
		protSequences=[]
		for s in sequences:
			protSequences.append(s)
		protSequencesAccepted = [x for x in protSequences if x.description.replace(" ", "").replace("\"", "") in acceptedSequences]
		dict_speciesInfo[str_speciesName]["newGenes"] = protSequencesAccepted
		###########################################################
		# Have to make sure gene names are not being duplicated,
		# which can be a problem with iterated runs, for example.
		###########################################################
		# get the list of existing gene names
		path_acceptedGff=path_resultsDir + "/" + str_speciesName + ".newGenes.gtf"
		path_geneNameConversionTable=path_augustusParsed + ".geneNamesConversion.txt"
		assignNames(str_speciesName, path_acceptedGff, path_geneNameConversionTable, protSequencesAccepted, dict_sequenceInfoById, path_augustusParsed, path_acceptedSequencesOut)
		# write everything out
		path_resultsFasta = path_resultsDir + "/" + str_speciesName + ".results.aa.fasta"
		path_resultsGff =path_resultsDir + "/" + str_speciesName + ".results.gtf"
		dict_speciesInfo[str_speciesName]["resultsgff"]=path_resultsGff
		# write out the gtf results file
		appendFileToFile(dict_speciesInfo[str_speciesName]["gff"], path_resultsGff)
		appendFileToFile(path_acceptedGff, path_resultsGff)
		# write out the fasta results file
		appendFileToFile(dict_speciesInfo[str_speciesName]["protein"], path_resultsFasta)
		appendFileToFile(path_acceptedSequencesOut, path_resultsFasta)
	print("\n6. Finishing up!")
	print(  "================")
	print("Results directory is " + path_resultsDir)
	for str_species in dict_speciesInfo:
		print("-------" + str(len(dict_speciesInfo[str_species]["newGenes"])) + " new genes found for " + str_species)

def compareOutputSequences(seq1, seq2):
	if seq1.replace(" ", "").replace("\"", "") == seq2.replace(" ", "").replace("\"", ""):
		return True
	else:
		return False

def makeNewProteome(path_oldProteome, path_predictedProteinSequences, path_newProteome):
	oldSeqs = [a for a in SeqIO.parse(path_oldProteome, "fasta") if a.seq != ""]
	newSeqs = [a for a in SeqIO.parse(path_predictedProteinSequences, "fasta") if a.seq != ""]
	SeqIO.write(oldSeqs + newSeqs, path_newProteome, "fasta")

def getOrthogroupsFile(path_aaDir):
	a = find("OrthologousGroups.csv", path_aaDir)
	if a:
		return a
	else:
		b = find("Orthogroups.csv", path_aaDir)
		if b: 
			return b
		else:
			print("Can't find Orthogroup output file... exiting...")
			sys.exit()
 
def getSingletonsFile(path_aaDir):
	a = find("OrthologousGroups_UnassignedGenes.csv", path_aaDir)
	if a:
		return a
	else:
		b = find("Orthogroups_UnassignedGenes.csv", path_aaDir)
		if b:
			return b
		else:
			print("Can't find Singletons output file... exiting...")
			sys.exit()

def assignNames(str_speciesName, path_acceptedGff, path_geneNameConversionTable, protSequencesAccepted, dict_sequenceInfoById, path_augustusParsed, path_acceptedSequencesOut):
	originalNames = [dict_sequenceInfoById[x].seqId for x in dict_sequenceInfoById if dict_sequenceInfoById[x].species == str_speciesName]
	originalNamesStubs = [x.split(".")[0] for x in originalNames ]
	allNames= originalNames + originalNamesStubs
	# for each new gene, give it a nice name and check that it hasn't been used before.
	counter=1
	print("updating names....")
	f = open(path_geneNameConversionTable, "w")
	writer = csv.writer(f, delimiter = '\t',quoting = csv.QUOTE_NONE, quotechar='')
	for s in protSequencesAccepted:
		newNameFound=False
		proposedGeneName=""
		while not newNameFound:
			proposedGeneName="orthofiller_g" + str(counter)	
			if not proposedGeneName in allNames:
				newNameFound=True
			counter = counter + 1
		writer.writerow([s.description , proposedGeneName])
		s.id=proposedGeneName
		s.description=proposedGeneName
		s.name=proposedGeneName
	# Write stuff out.
	f.close()
	SeqIO.write(protSequencesAccepted,  path_acceptedSequencesOut, "fasta")
	print("writing out results....")
	callFunction("rm -f "+path_acceptedGff + ";\
		while read line ; do \
			echo $line; \
			sourceId=`echo \"$line\" | cut -f1`; \
			replacementId=`echo \"$line\" | cut -f2 | sed -r \"s/ //g\"`;\
			tid=`echo $sourceId | sed -r \"s/.*transcript_id[= ]\\\"?([^\\\";]*)\\\"?;.*/\\1/g\"`;\
			gid=`echo $sourceId | sed -r \"s/.*gene_id[= ]\\\"?([^\\\";]*)\\\"?;.*/\\1/g\"`;\
			tnum=`echo $tid | sed -r \"s/^$gid\\.//g\"`;\
			grep -P \"(transcript_id[= ]\\\"$tid\\\")|(gene_id[= ]\\\"$gid\\\")|(transcript\\t.*\\t$tid$)|(gene\\t.*\\t$gid$)\" " + \
			path_augustusParsed + " | sed -r \"s/\\t$gid\\t/\\t$replacementId\\t/g\" \
						| sed -r \"s/\\t$tid\\t/\\t$replacementId\\.$tnum\\t/g\" \
						| sed -r \"s/transcript_id[= ]\\\"?$tid\\\"?; ?gene_id[= ]\\\"?$gid\\\"?;/transcript_id=\\\"$replacementId\\.$tnum\\\";gene_id=\\\"$replacementId\\\";/g\" | \
								sed -r 's/(.*)\\tpossibleOrthos.*/\\1/g' >> " + path_acceptedGff +";\
		done < " + path_geneNameConversionTable)

def annotateIntersectedOutput(path_hitsOgIntersectionFileName, path_hitsOgIntersectionFileNameAnnotated):
	callFunction("infile=\""+ path_hitsOgIntersectionFileName+"\"; outfile=\""+path_hitsOgIntersectionFileNameAnnotated+"\";\
		rm -f > $outfile; \
		awk -F \"\\t\" '$20 == \".\"' $infile | sed -r \"s/_id[= ]*[^\\\"]*\\\"/_id=\\\"/g\" | sed -r \"s/$/\\tmatch_none/g\" >> $outfile ;\
		awk -F \"\\t\" '$20 != \".\"' $infile | awk '$20 != $11' | sed -r \"s/_id[ =]*[^\\\"]*\\\"/_id=\\\"/g\" | sed -r \"s/$/\\tmatch_bad/g\" >> $outfile;\
		awk -F \"\\t\" '$20 != \".\"' $infile | awk '$20 == $11' | sed -r \"s/_id[ =]*[^\\\"]*\\\"/_id=\\\"/g\" | sed -r \"s/$/\\tmatch_good/g\" >> $outfile;")
	
def runAndParseAugustus(path_goodHits, path_genome, path_augustusOut, path_augustusParsedOut, path_fastaOut, path_augustusSpeciesName, path_hintsFile, path_sourcegff):
	print("augustus out is " + path_augustusOut)
	print("augustus parsed out is " + path_augustusParsedOut)
	runAugustus(path_goodHits, path_genome, path_augustusOut, path_augustusSpeciesName, path_hintsFile) #ql
	parseAugustusOutput(path_augustusOut, path_augustusParsedOut, path_fastaOut, path_sourcegff)


def combineIndirectAugustusResults(path_otherSpeciesResults, path_augustusParsedOut, path_fastaOut):
	print("Combining results for " + path_otherSpeciesResults +" into " + path_augustusParsedOut)
	function="predictions=\"" + path_otherSpeciesResults + "\"; gffout=\""+path_augustusParsedOut+"\"; fastaout=\""+path_fastaOut+"\"; " + """
		working=`mktemp -d`
		concat=`mktemp $working/concat.XXXXXXX`
		doubles=`mktemp $working/doubles.XXXXXXX`
		interim=`mktemp $working/interim.XXXXXXX`

		### Get genes as merged consensus regions ###
		cat  $predictions/*AugustusParsed.gff | grep -P "\\tgene\\t" | sort -k1,1V -k4,4n | bedtools merge -s -i - | cut -f1,2,3,5 | sed -r "s/\\t([^\\t]*)$/\\t.\\t.([^\\t]*)/g" |  perl -ne 'chomp; @l=split; printf "%s\\t%s\\t%s\\t.\\t.\\t%s\\tg%s\\n", $l[0], $l[1]-1, $l[2], $l[3], $., ' > $concat

		### Get only those gene regions which are predicted more than one time  ###
		for file in `find $predictions -type "f" -name "*AugustusParsed.gff"`; do base=`basename $file` ; bfile=`mktemp`; grep -P "\\tgene\\t" $file | sed -r "s/$/\\t$base/g" > $bfile; bedtools intersect -a $concat -b $bfile -wa -wb ; done | sort -k1,1V -k2,2n | rev | uniq -D -f11 | rev  > $doubles

		### Pick one gene from each region ###
		mkdir $working/split
		awk -v dir="$working/split" '{ f = dir "/" $1 "-" $2 "-" $3 ".tsv"; print  > f }' $doubles
		for file in `find $working/split -type "f"`; do shuf $file | sort -s -k13,13nr | head -n1 >> $interim; done
		almost=`mktemp $working/almost.XXXXXX`

		echo "concatenating genes from all sources..."
		while read line ; do file=`echo "$line" | cut -f18`; gid=`echo "$line" | cut -f16`; newgid=`echo "$line" | cut -f7`; grep -P "\\tgene\\t.*\\t$gid(\\t|\\z)|\\ttranscript\\t.*\\t$gid\\.t[^\\t]*(\\t|\\z)|\\t[^\\t]*gene_id[= ]\\"$gid\\";" $predictions/$file | sed -r "s/$/\\t$file/g" >> $almost; done < $interim

		### Give new names to all the genes and reflect this in the sequences file ###
		counter=1
		mkdir $working/bysource
		awk -v dir="$working/bysource" 'BEGIN { FS = "\\t" } { f = dir "/" $11 ".almost.tsv"; print  > f }' $almost

		fastatmp=`mktemp $working/fasta.XXXXXX`

		for almostsplit in `find $working/bysource -type "f"`; do
			echo "working on $almostsplit"
			gffsourcebase=`cut -f 11 $almostsplit | sort | uniq | head -n1`
			gffsource="$predictions/$gffsourcebase"
			fastasource=`echo "${gffsource%AugustusParsed.gff}AugustusParsed.sequences.fasta"`
			echo "the gff source is $gffsource"
			echo "the fasta source is $fastasource"
			for tid in `grep -P "\\ttranscript\\t" $almostsplit | cut -f9`; do
				###Edit names in the gff file
				gid=`echo $tid | sed -r "s/(g[^.]*)\\.t.*/\\1/g"`
				newgid="g$counter"
				newtid=`echo $tid | sed -r "s/g[^.]*(\.t.*)/$newgid\\1/g"`
				echo "old gene is called $gid, transcript $tid"
				echo "writing out new gene $newgid, transcript $newtid."
				sed -ri "s/(\\tgene\\t.*\\t)$gid(\\tpos.*)/\\1$newgid\\2/g" $almostsplit
				sed -ri "s/(\\ttranscript\\t.*\\t)$tid(\\tpos.*)/\\1$newtid\\2/g" $almostsplit
				sed -ri "s/transcript_id[= ]\\"$tid\\"; ?gene_id[= ]\\"$gid\\";/transcript_id \\"$newtid\\"; gene_id \\"$newgid\\";/g" $almostsplit

				#Get the relevant fasta sequence and change the names
				awk -v patt="gene_id[= ]\\"$gid\\";\\n" 'BEGIN {RS=">"} $0 ~ patt {print ">"$0}' $fastasource | grep -v "^$" | sed -r "s/transcript_id[= ]\\"$tid\\"; ?gene_id[= ]\\"$gid\\";/transcript_id \\"$newtid\\"; gene_id \\"$newgid\\";/g" >> $fastatmp
				counter_tmp=`echo $[counter + 1]`
				counter=$counter_tmp
			done
		done
		mv $fastatmp $fastaout
		cut -f1-10 $working/bysource/* > $gffout
		"""
	callFunction(function)

def runAugustus(path_goodHits, path_genome, path_augustusOut, path_augustusSpeciesName, path_hitsHintsGff):
	print("making hints file....")	
	makeHintsFile(path_goodHits, path_hitsHintsGff)
	print("running augustus, hints file: " + path_hitsHintsGff + "; augustus species: " + path_augustusSpeciesName + "; genome: " + path_genome)
	callFunctionQuiet("augustus --genemodel=complete --hintsfile=" + path_hitsHintsGff + \
			" --species=" + path_augustusSpeciesName + " " + path_genome + " > " + path_augustusOut)

def makeHintsFile(path_goodHits, path_hitsHintsGff):
	"""Converts our hits output file into a gff file which can be used to 
	   give hints to AUGUSTUS.
	"""
	callFunction("grep -v \"#\" " + path_goodHits + " | sed -r \"s/ +/\\t/g\" | perl -ne 'chomp;@l=split; printf \"%s\\tOrthoFiller\\texonpart\\t%s\\t%s\\t%s\\t%s\\t.\\torthogroup=%s;source=M\\n\", $l[0], $l[1], $l[2], $l[6], $l[5], $l[10]' | sed -r \"s/ +/\\t/g\"  > " + path_hitsHintsGff)

def parseAugustusOutput(path_augustusOutput, path_outputGff, path_outputFasta, path_sourcegff):
	function = "infile=\"" + path_augustusOutput + "\"; outfile=\"" +\
			path_outputGff + "\"; fastaout=\"" + path_outputFasta + "\"; sourcegff=\"" + path_sourcegff + "\";"+ \
			"""ot=`mktemp -d`; mkdir $ot/augsplit;
			#echo "parsing in $ot";
			awk -v RS="# start gene" -v ot="$ot" '{print "#"$0 > ot"/augsplit/augSplit."NR }' $infile
			mkdir $ot/success
			rm -f $fastaout
			rm -f $outfile
			find $ot/augsplit -type "f" | xargs grep -P "transcript supported by hints \(any source\): [^0]" | cut -f 1 -d ":" | xargs -I '{}' mv '{}' $ot/success/
			for file in `find $ot/success -type "f"`; do
				flatstring=`grep "#" $file | sed -r "s/# //g" | sed ':a;N;$!ba;s/\\n/ /g'`
				possibleOrthos=`echo "$flatstring" | sed -r "s/.*hint groups fully obeyed://g" | grep -oP "OG[0-9]{7}" | paste -sd, | sed -r "s/^,//g" | sed -r "s/,$//g"`
				grep -v "#" $file | grep -P "\\tAUGUSTUS\\t" | grep -v "^$" | sed -r "s/$/\\tpossibleOrthos=$possibleOrthos/g" >> $outfile
				id=`grep "transcript_id" $file | sed -r "s/.*(transcript_id \\"[^\\"]*\\"; gene_id \\"[^\\"]*\\";).*/\\1/g" | sort | head -n1 `
				#echo "printing sequences for $id"
				sequence=`echo "$flatstring" | sed -r "s/.*protein sequence = \[([A-Z ]*)\].*/\\1/g" | sed -r "s/[ \\t]//g"`
				echo ">$id###$sequence" >> $fastaout.tmp
			done
			sort -u $outfile > $outfile.tmp;

			grep CDS $outfile.tmp |  awk 'BEGIN {FS="\\t"}{if (b[$9]==""){b[$9]=$4; e[$9]=$5; c[$9]=$0}; if (b[$9] > $4){b[$9]=$4}; if(e[$9] < $5){e[$9]=$5}} END {for (i in b) {print c[i]"\\t"b[i]"\\t"e[i]}}' | awk  'BEGIN{OFS="\\t"; FS="\\t"} $4=$11; $5=$12' | sort -u | sed -r "s/CDS/gene/g" | cut -f1-10 > $outfile.tmp.genes

			bedtools intersect -a $outfile.tmp.genes -b $sourcegff -wa -wb | cut -f1-9 > $outfile.tmp.genes.remove
			bedtools subtract -A -a $outfile.tmp -b $outfile.tmp.genes.remove > $outfile.tmp.out

			while read line; do 
				tag=`echo "$line" | cut -f9`
				grep -v "$tag" $fastaout.tmp > $fastaout.tmptmp
				mv $fastaout.tmptmp $fastaout.tmp
			done < $outfile.tmp.genes.remove
			
			sed -r "s/###/\\n/g" $fastaout.tmp > $fastaout

			mv $outfile.tmp.out $outfile

			rm -f $outfile.tmp.genes $outfile.tmp.genes.remove $outfile.tmp $fastaout.tmp
			rm -r $ot
			"""
	callFunction(function)

def hintFscoreFilter(path_augustusParsed, path_hintFile, path_augustusParsedHintFiltered, num_threshold, path_augustusSequences, path_augustusSequencesHintFiltered):
	implementHintFscoreFilter(path_augustusParsed, path_hintFile, path_augustusParsedHintFiltered, num_threshold)
	extractFromFastaByName(path_augustusParsedHintFiltered, path_augustusSequences, path_augustusSequencesHintFiltered)

def implementHintFscoreFilter(path_augustusParsed, path_hintFile, path_outFile, num_threshold):
	function="augParsed=\"" + path_augustusParsed + "\"; hintFile=\"" + path_hintFile + "\"; of=\"" + path_outFile + "\"; threshold=\"" + str(num_threshold) + "\"" + """
	augParsedBed="$augParsed.bed"
	rm -f $of
	grep -P "\\tCDS\\t" $augParsed | sed -r "s/transcript_id \\"([^\\"]*)\\"; gene_id \\"([^\\"]*)\\";/transcript_id=\\"\\1\\";gene_id=\\"\\2\\";/g" | sort -u | perl -ne 'chomp; @l=split; printf "%s\\t%d\\t%d\\t.\\t.\\t%s\\t%s\\n", $l[0], $l[3]-1, $l[4], $l[6], $l[8]' > $augParsedBed

	hintsFileBed="$hintFile.bed"

	perl -ne 'chomp; @l=split; printf "%s\\t%d\\t%d\\t.\\t.\\t%s\\t%s\\n", $l[0], $l[3]-1, $l[4], $l[6], $l[8]' $hintFile > $hintsFileBed

	gids=`mktemp`
	sed -r "s/.*gene_id=\\"([^\\"]*)\\";.*/\\1/g" $augParsedBed | sort -u > $gids

	IFS='\n'
	nsuccess=0
	nfailure=0
	for gid in `cat $gids`; do
		#echo "checking hint scores for $gid from $augParsed"
		entrytmp=`mktemp`
		grep -P "gene_id[ =]\\"$gid\\"" $augParsedBed > $entrytmp
		compatibleHints=`bedtools intersect -wa -s -b $entrytmp -a $hintsFileBed | sort -u`
		gL=`awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $entrytmp`
		successes=""
		for hint in `echo "$compatibleHints" | cut -f7`; do
			hEntry=`mktemp`
			echo "$compatibleHints" | awk -v a="$hint" '$7 == a' > $hEntry
			hL=`awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $hEntry`
			iL=`bedtools intersect -s -a $hEntry -b $entrytmp | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'`
			successes=`awk -v gL="$gL" -v hL="$hL" -v iL="$iL" -v thresh="$threshold" 'BEGIN{hR=iL/hL; hP=iL/gL; hF=2*hR*hP/(hP+hR); if (hF >= thresh) {print "success"}}'`
		done
		if [ "$successes" != "" ]; then
			nsuccess=`echo $[$nsuccess+1]`
			grep -P "(gene_id[= ]\\"$gid\\";|\\t$gid\\t|\\t$gid\\.t.*\\t)" $augParsed | sort -u >> $of
		else
			nfailure=`echo $[$nfailure+1]`
		fi
	done
	echo "Finished implementing hint score filter for $augParsed. $nsuccess survivors and $nfailure non-survivors."
	sort -u $of > $of.tmp; mv $of.tmp $of

	rm $gids
	"""
	callFunction(function)


def extractFromFastaByName(path_gffFile, path_fastaFile, path_fastaOut):
	function="gff=\""+path_gffFile+"\"; fasta=\""+path_fastaFile+"\"; fastaout=\"" + path_fastaOut + "\"; " + """
	rm -f $fastaout; tids=`grep -P "\\ttranscript_id[ =]\\"[^\\"]*\\"; ?gene_id[ =]\\"[^\\"]*\\";" $gff | sed -r "s/.*\\ttranscript_id[ =]\\"([^\\"]*)\\"; ?gene_id[= ]\\"[^\\"]*\\";.*/\\1/g" | sort -u `
	IFS='\n'
	for tid in `echo "$tids"`; do
		#echo "fetching sequence for $tid"
		awk -v patt=".*transcript_id[= ][^a-zA-Z0-9_.]$tid[^a-zA-Z0-9._];.*;" 'BEGIN {RS=">"} $0 ~ patt {print ">"$0}' $fasta | grep -v "^$" >> $fastaout
	done"""
	callFunction(function)

def fetchSequences(path_gffIn, path_genome, path_cdsFastaOut, path_aaFastaOut, int_translationTable):
	path_nucleotideSequences=path_cdsFastaOut
	print("fetching nucleotide sequences for " + path_gffIn)
	callFunction("infile=" +  path_gffIn+ "; outfile=" + path_nucleotideSequences + "; genome=" + path_genome + """;
		tf=`mktemp -d`
		gffCds="$tf/gffCds"
		gffBed="$tf/gffBed"
		
		#Prepare the gff
		echo "preparing gff..."
		grep -vP "^$" $infile | awk '$3=="CDS"' > $gffCds
		cut -f1-8 $gffCds > $gffBed.1
		sed -r "s/.*transcript_id[ =]\\"?([^\\";]*)\\"?;?.*/\\1/g" $gffCds > $gffBed.2
		paste $gffBed.1 $gffBed.2 | perl -ne 'chomp; @l=split; printf "$l[0]\\t%s\\t$l[4]\\t$l[8]\\t.\\t$l[6]\\n", $l[3]-1' | sort -u | sort -k1,1V -k2,2n > $gffBed

		#Negative strand
		#echo "negative strand..."
		awk '$6=="-"' $gffBed > $gffBed.neg
		bedtools getfasta -name -s -fullHeader -fi $genome -fo $gffBed.neg.tab -bed $gffBed.neg -tab
		tac $gffBed.neg.tab | awk '{a[$1]=a[$1]""$2} END {for (i in a) {print ">"i"\\n"a[i]}}' > $gffBed.neg.fa

		#Then positive strand
		#echo "positive strand..."
		awk '$6=="+"' $gffBed > $gffBed.pos
		bedtools getfasta -name -s -fullHeader -fi $genome -fo $gffBed.pos.tab -bed $gffBed.pos -tab
		cat $gffBed.pos.tab | awk '{a[$1]=a[$1]""$2} END {for (i in a) {print ">"i"\\n"a[i]}}' > $gffBed.pos.fa

		cat $gffBed.pos.fa $gffBed.neg.fa | sed -r "s/^>(.*)$/£££>\\1###/g" | sed -r \"s/$/###/g\" | tr '\\n' ' ' | sed -r "s/£££/\\n/g" | sed -r "s/### //g" | grep -v XXX | grep -v "\*[A-Z]" | grep -v "###$" | sed -r "s/###/\\n/g" | grep -vP "^$" > $outfile

		#echo $tf
		rm -r $tf
		""")
	print("translating to protein...")
	sequences=SeqIO.parse(path_nucleotideSequences, "fasta")
	protSequences=[]
	for s in sequences:
		s_p = s
		s_p.seq = s.seq.translate(table=int_translationTable)
		protSequences.append(s_p)
	SeqIO.write(protSequences, path_aaFastaOut, "fasta")

def start(path_speciesInfoFile, path_orthoFinderOutputFile, path_singletonsFile, path_outDir, path_resultsDir, path_wDir, hitFilter, hintFilter, int_cores):
	######################################################
	# Read in the locations of the input files and the
	# orthofinder output.
	# MD-CC: will need to have a consistency check for this.:
	######################################################
        print("\n1.1. Reading and checking orthofinder output and sequence info files")
        print(  "===================================================================")
	dict_speciesInfo = readInputLocations(path_speciesInfoFile)
	dict_sequenceInfoById, orthogroups, singletons = readOrthoFinderOutput(path_orthoFinderOutputFile, path_singletonsFile, dict_speciesInfo)
	#####################################################
	# How many genes are there for each species? If any
	# species has less than 100, it needs special training.
	#####################################################
	firstPassMode=False
	print("\n1.2. Preparing gtf files for Augustus training")
        print(  "==============================================")
	pool = multiprocessing.Pool(int_cores)
	for str_species in dict_speciesInfo:
		sequences = [ dict_sequenceInfoById[x].seqId for x in dict_sequenceInfoById if dict_sequenceInfoById[x].species == str_species ]
		dict_speciesInfo[str_species]["needsTraining"] = False
		dict_speciesInfo[str_species]["indirectAugustus"] = False
		if len(sequences) < 100:
			firstPassMode=True
			dict_speciesInfo[str_species]["indirectAugustus"] = True
			dict_speciesInfo[str_species]["augustusSpecies"] = ""
			dict_speciesInfo[str_species]["gffForTraining"] = ""
		else:
			path_gff=dict_speciesInfo[str_species]["gff"]
			dict_speciesInfo[str_species]["needsTraining"] = True#ql
			path_gffForTraining = path_wDir + "/" + str_species + ".training.gff"
#qr			dict_speciesInfo[str_species]["augustusSpecies"]=commands.getstatusoutput("a=`find " + path_wDir+ "/Augustus/"+str_species+"/autoAugTrain -name \"tmp_opt*\" -exec stat {} --printf=\"%y\\t%n\\n\" \\;  | sort -t\"-\" -k1,1n -k2,2n -k3,3n | head -n1  | cut -f2`; echo ${a##*/} | sed -r \"s/tmp_opt_//g\"")[1]
			dict_speciesInfo[str_species]["augustusSpecies"]=str_species+ ".orthofiller." + datetime.datetime.now().strftime("%y%m%d") + "." + ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(9))
			dict_speciesInfo[str_species]["gffForTraining"] = path_gffForTraining
			async(pool, makeGffTrainingFile, args=(path_gff, path_gffForTraining))#ql
	pool.close()
	pool.join()
	if firstPassMode:
		path_firstPassOutDir=path_outDir + "/firstPass"
		makeIfAbsent(path_firstPassOutDir)
		run(dict_speciesInfo, dict_sequenceInfoById, orthogroups, singletons, path_firstPassOutDir, path_wDir, path_orthoFinderOutputFile, path_singletonsFile, int_cores, True, False, hitFilter, hintFilter)
		pool = multiprocessing.Pool(int_cores)
		for str_species in dict_speciesInfo:
			sequences = [ dict_sequenceInfoById[x].seqId for x in dict_sequenceInfoById if dict_sequenceInfoById[x].species == str_species ]
			dict_speciesInfo[str_species]["indirectAugustus"]=False
			if len(sequences) < 100:
				# Make a unique (statistically..!) name for this iteration of the training. This loses
				# some efficiency, but AUGUSTUS seems to have problems when we try to train to the same
				# set too many times. (and it complains quietly).
				dict_speciesInfo[str_species]["augustusSpecies"]=str_species+ ".orthofiller." + datetime.datetime.now().strftime("%y%m%d") + "." + ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(9))
				path_gff=dict_speciesInfo[str_species]["resultsgff"]
				#path_gff="/cellar/michael/OrthoFiller/testing/testSet_jgi_sgd/OrthoFiller_20160407_nomt_removed_100_generic_firstPass_new/firstPass/Sacce_S288C_genes_nomt.aa.fasta.removed_100.fasta.results.gtf"
				path_gffForTraining = path_wDir + "/" + str_species + ".training.gff"
				dict_speciesInfo[str_species]["gffForTraining"] = path_gffForTraining
				dict_speciesInfo[str_species]["needsTraining"] = True
				async(pool, makeGffTrainingFile, args=(path_gff, path_gffForTraining))
			else:
				dict_speciesInfo[str_species]["needsTraining"] = False
		pool.close()
		pool.join()
		######################################################
		# Run it
		######################################################
		run(dict_speciesInfo, dict_sequenceInfoById, orthogroups, singletons, path_resultsDir, path_wDir, path_orthoFinderOutputFile, path_singletonsFile, int_cores, False, True, hitFilter, hintFilter)
	else:
		run(dict_speciesInfo, dict_sequenceInfoById, orthogroups, singletons, path_resultsDir, path_wDir, path_orthoFinderOutputFile, path_singletonsFile, int_cores, False, False, hitFilter, hintFilter)

def unpackFitDistributionScript(path_scriptDestination):
	str_script='library("gamlss")\n\nargs <- commandArgs(TRUE)\n\n\nsourceF=args[1]\naltSourceF= args[2]\noutF=args[3]\n\nprint(paste("reading in source table: ", args[1], sep=""))\na <- read.table(sourceF, sep="\\t", header=FALSE)\n\nnames(a) <- c("hitChr", "hitStart", "hitEnd", "mystery1", "mystery2", "hitStrand", "eVal", "score", "bias", "hitSpecies", "hitOg", "targetChr", "targetStart", "targetEnd", "mystery5", "mystery6", "targetStrand", "geneLabel", "targetSpecies", "targetOg", "match")\n\na <- cbind(a, score_adj=a$score/(a$hitEnd - a$hitStart))\n\nh <- hist(a$score_adj, breaks=50, plot=FALSE)\nb <- h$breaks\n\na_none <- a[a$match=="match_none",]\na_good <- a[a$match=="match_good",]\n\n#Declare variables\ng_good = ""\ng_bad = ""\n\ng_prob_good = ""\ng_prob_bad = "" \n\n# Sampling 1000 data points makes the curve-fitting quicker and hardly affects the fit.\ngetGamlss <- function(theData) { print(theData$t); theData_s <- as.data.frame(sample(theData$t, 1000)); colnames(theData_s) <- c("t"); gamlss(t ~ 1, data=theData_s, family="ST1", method=RS(), gd.tol=10000000, c.cyc=0.001, control=gamlss.control(n.cyc=200)) } \n\nif(nrow(a_good) > 1000) {\n\tprint("source table is good, going ahead...")\n\ta_bad_og <- a[a$match=="match_bad",]\n\ta_bad_singleton <- a[a$match=="match_singleton",]\n\ta_bad <- rbind(a_bad_og, a_bad_singleton)\n\n\tprint("fitting good hits")\n\tgoodScores=as.data.frame(a_good$score_adj); colnames(goodScores) <- c("t")\n\t\n\tg_good <- getGamlss(goodScores)\n\tprint("fitting bad hits")\
\n\tbadScores=as.data.frame(a_bad$score_adj); colnames(badScores) <- c("t")\n\tg_bad <- getGamlss(badScores)\n\t\n\tg_prob_good =  nrow(a_good) / (nrow(a_bad) + nrow(a_good))\n\tg_prob_bad =  1 - g_prob_good\n\n} else {\n\tprint("source table too sparse, using aggregate distribution")\n\tz = read.table(altSourceF, sep="\\t", header=FALSE)\n\n\tnames(z) <- c("hitChr", "hitStart", "hitEnd", "mystery1", "mystery2", "hitStrand", "eVal", "score", "bias", "hitSpecies", "hitOg", "targetChr", "targetStart", "targetEnd", "mystery5", "mystery6", "targetStrand", "geneLabel", "targetSpecies", "targetOg", "match")\n\n\tz <- cbind(z, score_adj=z$score/(z$hitEnd - z$hitStart))\n\t\n\tz_good <- z[z$match=="match_good",]\n\tz_bad_og <- z[z$match=="match_bad",]\n	z_bad_singleton <- z[z$match=="match_singleton",]\n\tz_bad <- rbind(z_bad_og, z_bad_singleton)\n\n\tprint("fitting good hits")\t\n\tgoodScores=as.data.frame(z_good$score_adj); colnames(goodScores) <- c("t")\n	g_good <- getGamlss(goodScores)\n	\tprint("fitting bad hits")\n\tbadScores=as.data.frame(z_bad$score_adj); colnames(badScores) <- c("t")\n	g_bad <- getGamlss(badScores)\n\n	g_prob_good =  nrow(z_good) / (nrow(z_bad) + nrow(z_good))\n	g_prob_bad =  1 - g_prob_good\n\n}\n\nxg_fun <- function(x) { k=dST1(x, mu=g_good$mu.coefficients ,sigma=exp(g_good$sigma.coefficients), nu=g_good$nu.coefficients, tau=exp(g_good$tau.coefficients)) }\nxb_fun <- function(x) { k=dST1(x, mu=g_bad$mu.coefficients ,sigma=exp(g_bad$sigma.coefficients), nu=g_bad$nu.coefficients, tau=exp(g_bad$tau.coefficients)) }\n\ng_val <- function(x) { (xg_fun(x)*g_prob_good - xb_fun(x)*g_prob_bad) / (g_prob_good*xg_fun(x) + g_prob_bad*xb_fun(x)) }\n\nnone_g_scores = cbind(a_none, g_val=g_val(a_none$score_adj))\nnone_good = none_g_scores[none_g_scores$g_val >= 0,]\nnone_bad = none_g_scores[none_g_scores$g_val < 0,]\n\nprint(paste("writing to ", outF, sep=""))\n\nwrite.table(none_good, outF, quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\\t")\n'
	f=open(path_scriptDestination, "w")
	f.write(str_script)

def unpackFitDistributionScript_noFilter(path_scriptDestination):
	str_script='library("gamlss")\n\nargs <- commandArgs(TRUE)\n#args=c("all_coverage_annotated.orthogroupintersection.Ashgo1_1_GeneCatalog_proteins_20140830.aa.fasta.bed.hitTypes.new.tmp")\n\nsourceF=args[1]\naltSourceF= args[2]\noutF=args[3]\n\nprint(paste("reading in source table: ", args[1], sep=""))\na <- read.table(sourceF, sep="\\t", header=FALSE)\n\nnames(a) <- c("hitChr", "hitStart", "hitEnd", "mystery1", "mystery2", "hitStrand", "eVal", "score", "bias", "hitSpecies", "hitOg", "targetChr", "targetStart", "targetEnd", "mystery5", "mystery6", "targetStrand", "geneLabel", "targetSpecies", "targetOg", "match")\n\na <- cbind(a, score_adj=a$score/(a$hitEnd - a$hitStart))\n\nh <- hist(a$score_adj, breaks=50, plot=FALSE)\nb <- h$breaks\n\na_none <- a[a$match=="match_none",]\na_good <- a[a$match=="match_good",]\n\n#Declare variables\ng_good = ""\ng_bad = ""\n\ng_prob_good = ""\ng_prob_bad = "" \n\n# Sampling 1000 data points makes the curve-fitting quicker and hardly affects the fit.\n#getGamlss <- function(theData) { print(theData$t); theData_s <- as.data.frame(sample(theData$t, 1000)); colnames(theData_s) <- c("t"); gamlss(t ~ 1, data=theData_s, family="ST1", method=RS(), gd.tol=10000000, c.cyc=0.001, control=gamlss.control(n.cyc=200)) } \ngetGamlss <- function(theData) {}\n\nif(nrow(a_good) > 1000) {\n\tprint("source table is good, going ahead...")\n\ta_bad_og <- a[a$match=="match_bad",]\n\ta_bad_singleton <- a[a$match=="match_singleton",]\n\ta_bad <- rbind(a_bad_og, a_bad_singleton)\n\n\tprint("fitting good hits")\n\tgoodScores=as.data.frame(a_good$score_adj); colnames(goodScores) <- c("t")\n\t\n\tg_good <- getGamlss(goodScores)\n\tprint("fitting bad hits")\
		\n\tbadScores=as.data.frame(a_bad$score_adj); colnames(badScores) <- c("t")\n\tg_bad <- getGamlss(badScores)\n\t\n\tg_prob_good =  nrow(a_good) / (nrow(a_bad) + nrow(a_good))\n\tg_prob_bad =  1 - g_prob_good\n\n} else {\n\tprint("source table too sparse, using aggregate distribution")\n\tz = read.table(altSourceF, sep="\\t", header=FALSE)\n\n\tnames(z) <- c("hitChr", "hitStart", "hitEnd", "mystery1", "mystery2", "hitStrand", "eVal", "score", "bias", "hitSpecies", "hitOg", "targetChr", "targetStart", "targetEnd", "mystery5", "mystery6", "targetStrand", "geneLabel", "targetSpecies", "targetOg", "match")\n\n\tz <- cbind(z, score_adj=z$score/(z$hitEnd - z$hitStart))\n\t\n\tz_good <- z[z$match=="match_good",]\n\tz_bad_og <- z[z$match=="match_bad",]\n	z_bad_singleton <- z[z$match=="match_singleton",]\n\tz_bad <- rbind(z_bad_og, z_bad_singleton)\n\n\tprint("fitting good hits")\t\n\tgoodScores=as.data.frame(z_good$score_adj); colnames(goodScores) <- c("t")\n	g_good <- getGamlss(goodScores)\n	\tprint("fitting bad hits")\n\tbadScores=as.data.frame(z_bad$score_adj); colnames(badScores) <- c("t")\n	g_bad <- getGamlss(badScores)\n\n	g_prob_good =  nrow(z_good) / (nrow(z_bad) + nrow(z_good))\n	g_prob_bad =  1 - g_prob_good\n\n}\n\nnone_good = a_none\n\nprint(paste("writing to ", outF, sep=""))\n\nwrite.table(none_good, outF, quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\\t")\n'
	f=open(path_scriptDestination, "w")
	f.write(str_script)

def checkChromosomes(path_gff, path_genome):
	print("checking chromosomes...")
	res=commands.getstatusoutput("gff=\"" + path_gff + "\"; genome=\"" + path_genome + "\"; " + """
		a=`mktemp`; b=`mktemp`;
		grep ">" $genome | sed -r "s/>//g" > $a;
		cut -f1 $gff | sort -u > $b
		
		genchrdup=`sort $a | uniq -d | sort -u | wc -l`
		errMsg=""
		if [ "$genchrdup" != 0 ]; then
			errMsgTmp=`echo "Genome file $genome has duplicate chromosomes. Please adjust and try again."`
			errMsg=$errMsgTmp
		fi
		
		gffonly=`cat $a $a $b | sort | uniq -u | wc -l`
		if [ "$gffonly" != 0 ]; then
			errMsgTmp=`echo "$errMsg\\nGff file  $gff contains coordinates that do not exist in genome file $genome. Please adjust and try again."`
			errMsg=$errMsgTmp
		fi
		echo "$errMsg" """)[1]
	if res != "":
		sys.exit(res)

def checkSequences(path_gff, path_cds, path_aa):
	print("checking gff, fasta, and cds files for consistency")
	res=commands.getstatusoutput("gff=\"" + path_gff + "\"; cds=\"" + path_cds + "\"; aa=\"" + path_aa + "\"; " + """
		a=`mktemp -d`
		grep ">" $cds | sed -r "s/>//g" | sort -u > $a/cds
		grep ">" $aa | sed -r "s/>//g" | sort -u > $a/aa
		cut -f9 $gff | grep "transcript_id" | sed -r "s/.*transcript_id[= ]\\"([^\\"]*)\\".*/\\1/g" | sort -u  > $a/gff
		comm -23 $a/aa $a/cds > $a/noCdsEntries.err
		comm -23 $a/aa $a/gff > $a/noGffEntries.err
		mCds=`cat $a/noCdsEntries.err | wc -l `
		mGff=`cat $a/noGffEntries.err | wc -l `
		errMsg=""
		if [ "$mCds" != 0 ]; then
			errMsgTmp=`echo "$mCds entries from aa fasta file are missing in cds fasta file. Entries must be identically named."`
			errMsg=$errMsgTmp
		fi			
		if [ "$mGff" != 0 ]; then
			errMsgTmp=`echo "$errorMsg$mGff entries from aa fasta file are missing in gtf file $gff. Gtf entries must have attribute 'transcript_id \\"[sequence name]\\";'.\\n"`	
			errMsg=$errMsgTmp
		fi
		echo "$errMsg" """)[1]
	if res != "":
		sys.exit(res)

def prepareFromScratch(path_infile, path_outDir, int_cores):
	#Pull out Gtf files, extract dna and then rna
	#Then run orthofinder
	path_seqDir=path_outDir + "/sequences"
	path_cdsDir=path_seqDir + "/cds"
	path_aaDir=path_seqDir + "/aa"
	makeIfAbsent(path_seqDir)
	makeIfAbsent(path_cdsDir)
	makeIfAbsent(path_aaDir)
	path_speciesInfoFile=path_seqDir + "/inputs.csv"
	dict_basicInfo={}
	with open(path_infile) as p:
		# Ignore any commented lines, typically these are headers.
		data = csv.reader((row for row in p if not row.startswith('#')), delimiter="\t")
		for line in data:
			# Use genome name as key in dictionary
			key = os.path.basename(line[1])
			dict_basicInfo[key]={}
			dict_basicInfo[key]["gff"]	= checkFileExists(line[0])
			dict_basicInfo[key]["genome"]	= checkFileExists(line[1])
	with open(path_speciesInfoFile, "w") as f:
		writer = csv.writer(f, delimiter = '\t',quoting = csv.QUOTE_NONE, quotechar='')
		writer.writerow(["#protein", "gff", "genome", "cds"])
		for key in dict_basicInfo:
			path_gffIn=dict_basicInfo[key]["gff"]
			path_genome=dict_basicInfo[key]["genome"]
			checkChromosomes(path_gffIn, path_genome)
			path_cdsFastaOut=path_cdsDir+"/"+key+".cds.fasta"
			path_aaFastaOut=path_aaDir+"/"+key+".aa.fasta"
			fetchSequences(path_gffIn, path_genome, path_cdsFastaOut, path_aaFastaOut, 1)
			writer.writerow([path_aaFastaOut, path_gffIn, path_genome, path_cdsFastaOut])
	for path_file in glob.glob(path_aaDir + "/Results*"):
		deleteIfPresent(path_file)
	runOrthoFinder(path_aaDir, int_cores)
	path_orthoFinderOutputFile	= getOrthogroupsFile(path_aaDir)
	path_singletonsFile		= getSingletonsFile(path_aaDir)
	return path_speciesInfoFile, path_orthoFinderOutputFile, path_singletonsFile

def runOrthoFinder(path_aaDir, int_cores=16):
	if "ORTHOFINDER_DIR" in os.environ:
		finderString="python " + os.environ["ORTHOFINDER_DIR"] + "/orthofinder.py"		
	elif os.path.isfile(os.path.dirname(os.path.abspath(__file__)) + "/orthofinder.py"):
		finderString="python " + os.path.dirname(os.path.abspath(__file__)) + "/orthofinder.py"
	elif os.path.isfile("orthofinder.py"):
		finderString="python orthofinder.py"
	else:
		try:
			callFunctionQuiet("orthofinder -h")
			finderString="orthofinder"
		except OSError:
			sys.stderr.write("Error: Can't find orthofinder. Looked for orthofinder in the following order: OrthoFiller.py directory, execution directory, system PATH. Please ensure orthofinder is either installed and included in your PATH or that the orthofinder.py file is included in the same directory as the OrthoFiller.py file. Orthofinder can be downloaded from https://github.com/davidemms/OrthoFinder")
			sys.exit()
	# Call orthofinder
	version=tuple([int(i) for i in commands.getstatusoutput(finderString + " -h | grep -P \"version +[0-9\.]+\" | sed -r \"s/.*version ([0-9.]+).*/\\1/g\"")[1].split(".")])
	if version < (1,0,2):
		callFunction(finderString + " -f " + path_aaDir + " -t " + str(int_cores)) 
	elif version < (1,1,2):
		callFunction(finderString + " -f " + path_aaDir + " -t " + str(int_cores) + " -g")
	else:
		callFunction(finderString + " -f " + path_aaDir + " -t " + str(int_cores) + " -og")

####################################
############ Utilities #############
####################################

def merge_dicts(*dict_args):
	"""
	Given any number of dicts, shallow copy and merge into a new dict,
	precedence goes to key value pairs in latter dicts.
	"""
	result = {}
	for dictionary in dict_args:
		result.update(dictionary)
	return result

def async(pool, function, args):
	"""Run asynchronously
	"""
	return pool.apply_async(function, args=args) 

def suppressStdOut():
	"""http://thesmithfam.org/blog/2012/10/25/temporarily-suppress-console-output-in-python/#
	   The various programs spit out a lot of stuff. Sometimes it's helpful to suppress it.
	"""
	with open(os.devnull, "w") as devnull:
		old_stdout = sys.stdout
		sys.stdout = devnull
		try:
			yield
		finally:
			sys.stdout = old_stdout
def find(name, path):
	"""Find the relative path of a named file in a folder (returns the first one it finds)
	"""
	for root, dirs, files in os.walk(path):
		if name in files:
			return os.path.join(root, name)

def makeIfAbsent(path_dir):
	try:
		os.makedirs(path_dir)
	except OSError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path_dir):
			pass
		else:
			raise

def appendFileToFile(path_source, path_dest):
	with open(path_dest, "a") as f:
		with open(path_source, "r") as g:
			f.write(g.read())

def deleteIfPresent(path):
	"""Delete a folder which may or may not exist
	"""
	try:
		if os.path.isdir(path):
			callFunction("rm -r " + path)
		else:
			os.remove(path)
	except OSError:
		pass

def checkFileExists(path_file):
	if not os.path.exists(path_file):
		sys.exit("File does not exist: " + path_file)
	else:
		return path_file

####################################
########### Entry code #############
####################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Run OrthoFiller")
	parser.add_argument("--noHintFilter", action="store_true", dest="noHintFilter", default=False)
	parser.add_argument("--noHitFilter", action="store_true", dest="noHitFilter", default=False)
	parser.add_argument("-c", "--cores", metavar="cores", help="The maximum number of cores you wish to use", dest="CO", default=1)
	parser.add_argument("-o", "--outdir", metavar="outdir", help="The output directory", dest="OD", default="")
	parser.add_argument("-i", "--infoFiles", metavar="info", dest="IN")
	parser.add_argument("--prep", help="Input data in pre-prepared form", dest="prep", action="store_true")
	parser.add_argument("-g", "--orthogroups", metavar="orthogroups", help="An orthofinder output file (orthogroups)", dest="OG")
	parser.add_argument("-s", "--singletons", metavar="singletons", help="An orthofinder output file (singles)", dest="SN")
	parser.add_argument("-t", "--translationtable", metavar="transtable", help="Which translation table to use", dest="TT")
	parser.add_argument("--checkPrograms", action="store_true", dest="checkOnly", default=False)

	args = parser.parse_args()
	prep=args.prep
	
	print("\n0.1. Checking installed programs")
	print(  "================================")
	checkShell()
	if args.checkOnly:
		print("Finished checking installed programs. Everything looks fine.")
		sys.exit()

	#Check existence and non-confliction of arguments
	if args.IN == None:
		sys.exit("Input file list -i required.")
	if prep:
		if (args.IN == None) | (args.OG == None) | (args.SN == None):
			sys.exit("Option --prep requires options -i [input file list], -g [orthogroups file], and -s [singletons file].")
	else:
		if (args.OG != None) | (args.SN != None):
			sys.exit("Options -g and -s can only be used with option --prep for pre-prepared data.")

	path_outDir = args.OD
	int_cores = args.CO

	print("\n0.2. Checking and unpacking input data")
	print(  "======================================")

	path_resultsDir, path_wDir = prepareOutputFolder(path_outDir)
	
	#If the data isn't pre-prepared, we must prepare it.
	#Else simply check each file exists. Later we will make sure every entry in the info file exists.
	if not prep:
		path_speciesInfoFile, path_orthoFinderOutputFile, path_singletonsFile = prepareFromScratch(args.IN, path_outDir, int_cores)
	else:
		path_orthoFinderOutputFile = checkFileExists(args.OG)
		path_singletonsFile = checkFileExists(args.SN)
		path_speciesInfoFile = checkFileExists(args.IN)

	start(path_speciesInfoFile, path_orthoFinderOutputFile, path_singletonsFile, path_outDir, path_resultsDir, path_wDir, not args.noHitFilter, not args.noHintFilter, int(int_cores))

