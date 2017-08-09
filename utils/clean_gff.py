#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import csv
import sys
import os
import re

##########################
# 
# Usage instructions:
# 
# python clean_gff.py in_file.gff3 out_file.gff3
#
##########################

def readCsv(path_csv):
	with open(path_csv, "r") as p:
		data = list(csv.reader((row for row in p if not row.startswith('#')), delimiter="\t"))
	return data

def write(data, path_out):
	with open(path_out, "w") as o:
		datawriter = csv.writer(o, delimiter = '\t',quoting = csv.QUOTE_NONE, quotechar='')
		datawriter.writerows(data)

def getParents(lines):
	res_parents ={}
	res_lines = {}
	res_children = {}
	for l in lines:
		l_id = re.sub(r".*ID=([^;]*);.*",r"\1", l[8])
		l_pa = re.sub(r".*Parent=([^;]*);.*",r"\1", l[8])
		res_parents[l_id] = l_pa
		res_children[l_pa] = res_children.get(l_pa, []) + [l_id]
		res_lines[l_id]  = res_lines.get(l_id, []) + [l]
	return res_parents, res_lines, res_children

def deDup(dlist):
	nlist = []
	for i in dlist:
		if i not in nlist:
			nlist.append(i)
	return nlist


def go(infile, outfile):
	# Grab all the bits of gff
	print "Reading gff..."
	gff   = [a for a in readCsv(infile) if not a == []]
	genes = [a for a in gff if a[2].lower() == "gene"]
	mrna  = [a for a in gff if a[2].lower() == "mrna"]
	exon  = [a for a in gff if a[2].lower() == "exon"]
	cds   = [a for a in gff if a[2].lower() == "cds"]
	# We only want well-formed genes, that is genes with all of gene, mrna, cds, and exon attributes.
	# ----gene
	print "Grabbing genes..."
	gene_lines = {}
	for g in genes:
		g_id = re.sub(r".*ID=([^;]*);.*",r"\1", g[8])
		# There shouldn't be more than one gene line.
		gene_lines[g_id] = [g]
	# ----mRNA
	print "Grabbing mrna..."
	mrna_parents, mrna_lines, gene_children = getParents(mrna)
	# ----CDS
	print "Grabbing cds..."
	cds_parents, cds_lines, mrna_children_cds = getParents(cds)
	# ----exon
	print "Grabbing exons..."
	exon_parents, exon_lines, mrna_children_exon = getParents(exon)
	# Now cross check to ensure only genes with mrna, exon, and cds children are obtained. Genes may have multiple mrna children.
	# Note that we don't check that the cds and exons actually match the mrna and gene in their contents.
	res = []
	genecount = len(gene_lines.keys())
	survivors = 0
	for i, g_id in enumerate(gene_lines.keys()):
		sys.stdout.write("\rChecking gene " + str(i+1) + " of " + str(genecount))
		sys.stdout.flush()
		if not g_id in gene_children: continue
		g_mrna = gene_children[g_id]
		added = False
		for m in g_mrna:
			# Check that each mrna has cds's and exons.
			if (not m in mrna_children_exon) or (not m in mrna_children_cds): continue
			m_exon = deDup(mrna_children_exon[m])
			m_cds  = deDup(mrna_children_cds[m])
			if not added:
				res += gene_lines[g_id]
				survivors += 1
				added = True
			res += mrna_lines[m]
			for e in m_exon:
				res += exon_lines[e]
			for c in m_cds:
				res += cds_lines[c]
	print "\nThere were " + str(survivors) + " survivors and " + str(genecount - survivors) + " non-survivors."
	write(res, outfile)

if __name__ == '__main__':
	args = sys.argv[1:]
	infile = args[0]
	outfile = args[1]
	go(infile, outfile)
