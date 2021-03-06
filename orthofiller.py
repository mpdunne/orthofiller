#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2017 Michael Dunne
#
# OrthoFiller version 1.1.1
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
libnames = ["csv", "re", "os", "sys", "signal", "shutil", "itertools", "Bio", "subprocess", \
            "multiprocessing", "datetime", "tempfile", "random", "string", "commands", \
            "errno", "argparse", "time", "glob"]

for libname in libnames:
    try:
        lib = __import__(libname)
    except ImportError as e:
		errors.append(e)
    else:
        globals()[libname] = lib

try:
	from Bio.Align.Applications import MafftCommandline
except ImportError as e:
	errors.append(e)
try:
	from Bio.SeqRecord import SeqRecord
except ImportError as e:
	errors.append(e)
try:
	from Bio import SeqIO
except ImportError as e:
	errors.append(e)

#sublibs = [["Bio", "SeqIO"], ["Bio.Align.Applications", "MafftCommandline"], \
#           ["Bio.SeqRecord", "SeqRecord"]]

#for libpair in sublibs:
#    try:
#        lib = __import__(libpair[0], fromlist=libpair[1])
#    except ImportError as e:
#		errors.append(e)
#    else:
#        globals()[libpair[1]] = lib

if errors:
	exitstring = "Missing modules :(\nThe following module errors need to be resolved before running OrthoFiller:"
	for error in errors: exitstring += "\n-- " + str(error)
	exit(exitstring)
	
########################################################
################ Check Linux Packages ##################
########################################################

def check(boolVal, trueMsg=" - ok", falseMsg=" - failed", trueExtra="", falseExtra=""):
	return outputBool(True, trueMsg, trueExtra) if boolVal else outputBool(False, falseMsg, falseExtra)

def outputBool(boolVal, msg, msgExtra):
	print msg
	if msgExtra: print msgExtra
	return boolVal

def canRunCommand(command, qAllowStderr = False):
	sys.stdout.write("Test can run \"%s\"" % command)
	capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout = [x for x in capture.stdout]
	stderr = [x for x in capture.stderr]
	return check(len(stdout) > 0 and (qAllowStderr or len(stderr) == 0))

def canRunSpecific(line, lineFormatted):
	return True if canRunCommand(line) else outputBool(False, "ERROR: Cannot run " + lineFormatted, \
		"    Please check "+ lineFormatted +" is installed and that the executables are in the system path\n")

def returnsExpected(message, command, contentsIn, expectedOut):
	sys.stdout.write("Test can run \""+message+"\"\n")
	path_tf1 = tempfile.mktemp(); path_tf2 = tempfile.mktemp()
	write(contentsIn, path_tf1)
	callFunction(command + " " + path_tf1 + " > " + path_tf2)
	res = read(path_tf2)
	return res == expectedOut	
		
def canRunMinusH(package, packageFormatted):
	return canRunSpecific(package + " -h", packageFormatted)

def canRunHelp(package, packageFormatted):
	return canRunSpecific(package + " -help", packageFormatted)

def canRunMan(package, packageFormatted):
        return canRunSpecific("man " + package, packageFormatted)

def canRunBlank(package, packageFormatted):
	return canRunSpecific(package, packageFormatted)

def canRunAwk():
	return returnsExpected("awk '{print 4 + 5}'", "awk '{print 4 + 5}'", "a", "9\n")

def canRunMafft():
	return returnsExpected("awk 'mafft", "mafft --quiet", ">1\nAAABBBCCCDDD\n>2\nAAACCCDDD\n", ">1\nAAABBBCCCDDD\n>2\nAAA---CCCDDD\n")
	
def canRunOrthoFinder():
	success = False
	if "ORTHOFINDER_DIR" in os.environ:
		success = canRunCommand("python " + os.environ["ORTHOFINDER_DIR"] + "/orthofinder.py -h")	
        elif os.path.isfile(os.path.dirname(os.path.abspath(__file__)) + "/orthofinder.py"):
		success = canRunCommand("python " + os.path.dirname(os.path.abspath(__file__)) + "/orthofinder.py -h ")
        elif os.path.isfile("orthofinder.py"):
                success = canRunCommand("python orthofinder.py -h")
        else:
                success = canRunCommand("orthofinder - h")
	if not success:
		print("    Cannot find the orthofinder.py file. Either\n        1) Set the orthofinder location as an environment variable using \"export ORTHOFINDER_DIR=dir\", where dir is the path to the directory containing orthofinder.py\n        2) include orthofinder.py in the same directory as OrthoFiller; or\n        3) install an orthofinder executable in your system PATH\n    Orthofinder can be downloaded from https://github.com/davidemms/OrthoFinder")
	return success

def canRunGeneric(package, packageFormatted):
	sys.stdout.write("Test can run \"" + package + "\"")
	runit = subprocess.call("type " + package, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
	return check(runit, falseExtra="ERROR: Cannot run " + packageFormatted+"\n    Please check "+ packageFormatted +" is installed and that the executables are in the system path\n")
		
def canRunAugTrain():
	sys.stdout.write("Test can run \"autoAugTrain.pl\"")
	runit = subprocess.call("type " + "autoAugTrain.pl", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
	return check(runit, falseExtra="    The path to the AUGUSTUS scripts folder must be exported into the path before running, i.e.\n        export PATH=$PATH:path_to_aug_scripts_folder")

def canRunR():
	sys.stdout.write("Test can run \"R, Rscript\"")
	runit  = subprocess.call("type " + "R", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
	rsunit = subprocess.call("type " + "Rscript", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
	# Check R can run
	path_r = tempfile.mktemp()
	write("print(\"hello world\")", path_r)
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
	write("library(\"gamlss\"); search(); print(\"blue whale\")", path_r)
	capture = subprocess.Popen("Rscript " + path_r, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	gamlFound = any(["blue whale" in i for i in capture.stdout])
	stderr = [x for x in capture.stderr]
	if not gamlFound:
		print(" - failed")
		print("    R was unable to run. Please make sure you have the \"gamlss\" package installed.")
		return False
	os.remove(path_r)
	a=commands.getstatusoutput("Rscript " + path_r)[1]
	return check(not "WARNING: ignoring environment value of R_HOME" in a, \
			falseExtra="    R R_HOME WARNING encountered: please call \"unset R_HOME\" in bash before running OrthoFiller\n")

def mockGtf(chrom, coords, types, strand, gid, tid):
	sline = [chrom, "fake", "#", "#", "#", ".", strand, ".", "gene_id \""+gid+"\"; transcript_id \""+tid+"\";"]
	if len(coords) != len(types): raise ValueError()
	res = []
	for i in range(len(coords)):
		nsline = sline[:]
		nsline[2] = types[i]
		nsline[3] = str(coords[i][0])
		nsline[4] = str(coords[i][1])
		res.append(nsline)
	return res

def canRunBedTools():
	""" New versions of bedtools often seem to break things. Make extra sure that everything works as it needs to.
	    This is not an exhaustive unit test.
	"""
	path_td = tempfile.mkdtemp()
	sys.stdout.write("Testing bedtools")
	try:
		# Mock genome
		path_mockGenome = path_td + "/genome.fasta"
		seq1 = SeqRecord(Bio.Seq.Seq("CCCACGTAGCTAAGTGAATAAGTAGCCGCGCTCGACACACAGTGATGGATACGGCAGCTGGCACCAACAAGGAGCGTGGGATGTACGCTATTTTGGCGAT"), id="chr1", description="")
		seq2 = SeqRecord(Bio.Seq.Seq("GGGCTCGCAACGGTTGGCCACGGCGTCTCTCATATTCGTTTAGTAGATTAACTTTCAAGATGAAAACGAGCTCCTACTAGAGAATCTCTCAGCGCACTAT"), id="chr2", description="")
		SeqIO.write([seq1, seq2], path_mockGenome, "fasta")
		mockGtf1 = mockGtf("chr1", [[10,39],[61,90],[10,12],[88,90]], ["CDS", "CDS", "start_codon", "stop_codon"], "+", "posGene", "posGene.t1")
		mockGtf2 = mockGtf("chr2", [[15,29],[46,60],[58,60],[15,17]], ["CDS", "CDS", "start_codon", "stop_codon"], "-", "negGene", "negGene.t1")
		mockGtf3 = mockGtf("chr1", [[5,34], [56,85],[5,7],  [83,85]], ["CDS", "CDS", "start_codon", "stop_codon"], "+", "posGene", "posGene.t1")
                mockGtf4 = mockGtf("chr2", [[10,24],[41,55],[53,55],[10,12]], ["CDS", "CDS", "start_codon", "stop_codon"], "-", "negGene", "negGene.t1")
		mockGtf5 = mockGtf("chr1", [[10,30],[20,50],[10,12],[48,50]], ["CDS", "CDS", "start_codon", "stop_codon"], "+", "posGene", "posGene.t1")
                mockGtf6 = mockGtf("chr2", [[50,80],[70,90],[88,90],[50,52]], ["CDS", "CDS", "start_codon", "stop_codon"], "-", "negGene", "negGene.t1")
		path_mockGtf1 = writeCsv(mockGtf1 + mockGtf2, path_td + "/genes1.gtf")
		path_mockGtf2 = writeCsv(mockGtf3 + mockGtf4, path_td + "/genes2.gtf")
		path_mockGtf3 = writeCsv(mockGtf5 + mockGtf6, path_td + "/genes3.gtf")
		# Test merge
		path_merged = path_td + "/merged.gtf"
		callFunction("sort -k1,1V " + path_mockGtf3+ " | grep CDS | bedtools merge -i - -s > " + path_merged)
		with open(path_merged, "r") as f:
			if not f.read() == 'chr1\t9\t50\t+\nchr2\t49\t90\t-\n':
				raise ValueError()
		# Test intersect
		path_intersected = path_td + "/intersected.gtf"
		callFunction("bedtools intersect -wa -wb -a " + path_mockGtf1 + " -b "+ path_mockGtf2 + " > " + path_intersected)
		isec = readCsv(path_intersected)
		expected =  [mockGtf1[0] + mockGtf3[0], mockGtf1[1] + mockGtf3[1], mockGtf1[1] + mockGtf3[3], mockGtf1[2] + mockGtf3[0]]
		expected += [mockGtf2[0] + mockGtf4[0], mockGtf2[1] + mockGtf4[1], mockGtf2[1] + mockGtf4[2], mockGtf2[3] + mockGtf4[0]]
		if not isec == expected: raise ValueError()
		# Test subtract
		path_subtracted = path_td + "/subtracted.gtf"
		callFunction("bedtools subtract -a " + path_mockGtf1 + " -b " + path_mockGtf2 + " > " + path_subtracted)
		subGtf1 = mockGtf("chr1", [[34,39],[85,90],[88,90]], ["CDS", "CDS", "stop_codon"], "+", "posGene", "posGene.t1")
		subGtf2 = mockGtf("chr2", [[24,29],[55,60],[58,60]], ["CDS", "CDS", "start_codon"], "-", "negGene", "negGene.t1")
		if not readCsv(path_subtracted) == subGtf1 + subGtf2: raise ValueError()
		# Test getfasta
		path_getFasta = path_td + "/getfasta.fa"
		path_getFastaAa = path_td + "/getfasta.aa.fa"
		fetchSequences(path_mockGtf1, path_mockGenome, path_getFasta, path_getFastaAa, 1)
		with open(path_getFasta, "r") as f:
			if not f.read() == '>posGene.t1\nCTAAGTGAATAAGTAGCCGCGCTCGACACAGCACCAACAAGGAGCGTGGGATGTACGCTA\n>negGene.t1\nTCTTGAAAGTTAATCGAGACGCCGTGGCCA\n':
				raise ValueError()
		#######
		deleteIfPresent(path_td)
		print(" - ok")
		return True
	except ValueError:
		callFunction("rm -r " + path_td)
		print(" - fail")
		print("The version of BedTools you're using is causing problems with OrthoFiller.\nSometimes BedTools updates prevent BedTools itself from running properly.\nTry downloading bedtools v2.25.0 and run OrthoFiller again.")
		return False

def checkShell():
	""" Run quick checks to ensure the relevant packages are installed, and that they do
		the things we need them to do (applies in particular to bedtools and to R).
		- perl          - awk          -orthofinder
		- echo          - sed          - bedtools
		- cut           - grep         - R/Rscript
		- sort          -tac/cat       - hmmer/hmmbuild
		- augustus/autoAugTrain.pl
		"""
	checks = []
	for program in ["sed", "echo", "cut", "sort", "grep", "uniq", "tac", "cat", "Rscript", "mktemp"]:
		checks.append(canRunMan(program, program))
	for program in ["nhmmer", "hmmbuild", "makehmmerdb", "bedtools", "mcl"]:
		checks.append(canRunMinusH(program, program))
	checks.append(canRunMinusH("perl", "PERL"))
	checks.append(canRunMinusH("makeblastdb", "blast+"))
	checks.append(canRunBlank("augustus", "augustus"))
	#Check presence of orthofinder
	checks += [canRunAwk(), canRunR(), canRunMafft(), canRunAugTrain(), canRunOrthoFinder(), canRunBedTools()]
	if not all(checks):
		print("\nSome programs required to run orthoFiller are not installed or are not callable.\nPlease ensure all of the above programs are installed and in the system path.")
		sys.exit()

####################################
############ Utilities #############
####################################

def stage(str_step):
	print("\n" + str_step)
	print("="*len(str_step))

def callFunction(str_function):
	"""Call a function in the shell
	"""
	subprocess.call([str_function], shell = True)

def callFunctionQuiet(str_function):
	"""Call a function in the shell, but suppress output.
	"""
	with open(os.devnull, 'w') as FNULL: subprocess.call([str_function], shell = True, stdout=FNULL, stderr=subprocess.STDOUT)

def callFunctionMezzoPiano(str_function):
	"""Call function in the shell, and suppress output but not stderr.
	"""
	with open(os.devnull, 'w') as FNULL: subprocess.call([str_function], shell = True, stdout=FNULL)

def grabLines(command):
	return commands.getstatusoutput(command)[1].split("\n")

def readCsv(path_csv, ignoreHash=False):
	with open(path_csv, "r") as p: data = list(csv.reader((row for row in p if (ignoreHash or not row.startswith('#'))), delimiter="\t"))
	return data

def writeCsv(data, path_csv):
	with open(path_csv, "w") as f:
		writer = csv.writer(f, delimiter = '\t',quoting = csv.QUOTE_NONE, quotechar='')
		writer.writerows(data)
	return path_csv
		
def write(strObj, path_file, tag="w"):
	with open(path_file, tag) as f: f.write(strObj)

def read(path_file, tag="r"):
	with open(path_file, tag) as f: return f.read()

def flattenLol(lol):
        res = []
        for l in lol: res += l
        return res

def move(path_orig, path_target):
	callFunction("mv " + path_orig + " " + path_target)

def bashIt(bashLine, path_inFile, path_outFile=""):
	if path_outFile:
		callFunction(bashLine + " " + path_inFile + " > " + path_outFile)
	else:
		path_tmp = tempfile.mktemp()
		callFunction(bashLine + " " + path_inFile + " > " + path_tmp)
		move(path_tmp, path_inFile)

def init_worker():
	signal.signal(signal.SIGINT, signal.SIG_IGN)

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

def concatFiles(list_files, path_outfile, removeOrig=False):
	with open(path_outfile,'wb') as o:
		for path_indi in list_files:
			with open(path_indi,'rb') as p:
				shutil.copyfileobj(p, o)
			if removeOrig: os.remove(path_indi)

def concatGlob(str_patt, path_outfile, removeOrig=False):
	concatFiles(glob.glob(str_patt), path_outfile, removeOrig)

def exit(msg):
	sys.exit(msg)

def merge_dicts(*d_args):
	"""
	Given any number of dicts, shallow copy and merge into a new dict,
	precedence goes to key value pairs in latter dicts.
	"""
	result = {}
	for dictionary in d_args:
		result.update(dictionary)
	return result

def runAndTrackJobs(jobs, int_cores, message="", outMsg=True, ret = False, silent=False):
	pool = multiprocessing.Pool(int_cores, init_worker)
	reporters = {}
	for j in jobs:
		reporters[j] = async(pool, jobs[j][0], args = jobs[j][1])
	pool.close()
	trackProgress(reporters, message, outMsg, silent)
	pool.join()
	if ret:
		return reporters
#	try:
 #               time.sleep(1)
  #      except KeyboardInterrupt:
   #             pool.terminate()
    #            pool.join()
     #   else:
      #          pool.close()
       #         pool.join()


def trackProgress(jobs, prefix="", outMsg=True, silent=False):
	str_total = str(len(jobs))
	while True:
		k=[jobs[j].ready() for j in jobs]
                if not silent: sys.stdout.write("\r    " + prefix + str(sum(k)) + " out of " + str_total +  " orthogroups processed")
                sys.stdout.flush()
                if all(k):
			if outMsg and not silent: print("\n    Finished processing orthogroups")
			break
		time.sleep(2)

def runJobs(jobs, int_cores, ret = False):
	pool = multiprocessing.Pool(int_cores, init_worker)
	reporters = []
	for j in jobs:
		reporters.append(async(pool, j[0], args = j[1]))
	pool.close()
	pool.join()
	if ret:
		return reporters
#	try:
 #		time.sleep(5)
  #	except KeyboardInterrupt:
   #		pool.terminate()
    #		pool.join()
     #	else:
      #		pool.close()
       #	pool.join()

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

def makeIfAbsent(path_dir, fresh=False):
	try:
		if fresh: deleteIfPresent(path_dir)
		os.makedirs(path_dir)
		return path_dir
	except OSError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path_dir):
			return path_dir
		else:
			raise

def blankFile(path_file):
	open(path_file,"w").close()
	return path_file

def appendFileToFile(path_source, path_dest):
	with open(path_dest, "a") as f:
		with open(path_source, "r") as g: f.write(g.read())

def deleteIfPresent(path):
	"""Delete a folder which may or may not exist
	"""
	try:
		if os.path.isdir(path):
			callFunction("rm -r " + path)
		else:
			os.remove(path)
		return path
	except OSError:
		pass

def checkFileExists(path_file):
	if not os.path.exists(path_file):
		exit("File does not exist: " + path_file)
	else:
		return path_file

######################################
############ Preparation #############
######################################

def prepareOutputFolder(path_outDir):
	if path_outDir == "":
		str_prefix = "OrthoFillerOut_" + datetime.datetime.now().strftime("%y%m%d") + "_RunId_"
		path_outDir = tempfile.mkdtemp(prefix=str_prefix, dir=".")
	path_wDir       = makeIfAbsent(path_outDir + "/working")
	path_ogDir      = makeIfAbsent(path_wDir + "/orthogroups")
	path_resultsDir = makeIfAbsent(path_outDir + "/results")
	return path_resultsDir, path_wDir


def checkChromosomes(path_gtf, path_genome):
	genomeChr = grabLines("grep \">\" " + path_genome + " | sed -r \"s/>//g\"")
	gtfChr    = grabLines("cut -f1 " + path_gtf + " | sort -u")
	errMsg    = ""
	if len(genomeChr) != len(set(genomeChr)):
		errMsg += "Genome file $genome has duplicate chromosomes. Please adjust and try again."
	if not all([g in genomeChr for g in gtfChr]):
		errMsg += "\nGtf file " + path_gtf + " contains coordinates that do not exist in genome file $genome. Please adjust and try again."
	if errMsg: exit(errMsg)

def checkSequences(path_gtf, path_cds, path_aa):
	cdsEntries = grabLines("grep \">\" " + path_cds + " | sed -r \"s/>//g\" | sort -u")
	aaEntries  = grabLines("grep \">\" " + path_aa + " | sed -r \"s/>//g\" | sort -u")
	gtfEntries = grabLines("cut -f9 " + path_gtf + " | grep \"transcript_id\" | sed -r \"s:.*transcript_id[= ]\\\"([^\\\"]*)\\\".*:\\1:g\" | sort -u")
	errMsg = ""
	missingCds = [a for a in aaEntries if not a in cdsEntries]
	if missingCds:
		errMsg += str(len(missingCds)) + " entries from aa fasta file are missing in cds fasta file. Entries must be identically named."
	missingAa = [a for a in aaEntries if not a in gtfEntries]
	if missingAa:
		errMsg += "\n"+str(len(missingCds)) + " entries from aa fasta file are missing in gtf file $gtf. gtf entries must have attribute 'transcript_id \\\"[sequence name]\\\";'."
	if errMsg: exit(errMsg)

def grabBasicInfo(line, sptype):
	return {"gtf": checkFileExists(line[0]), "genome": checkFileExists(line[1]), "type": sptype}

def prepareFromScratch(path_infile, path_outDir, int_cores, path_refFile = ""):
	#Pull out gtf files, extract dna and then rna
	#Then run orthofinder
	useReference = not path_refFile == ""
	path_seqDir = makeIfAbsent(path_outDir + "/sequences")
	path_cdsDir = makeIfAbsent(path_seqDir + "/cds")
	path_aaDir  = makeIfAbsent(path_seqDir + "/aa")
	path_speciesInfoFile = path_seqDir + "/inputs.csv"
	path_refInfoFile     = path_seqDir + "/inputs.ref.csv"
	d_basicInfo={}
	data = readCsv(path_infile)
	if useReference: rel = readCsv(path_refFile)
	for line in data:
		if not ''.join(line).strip(): continue
		# Use genome name as key in dictionary
		key = os.path.basename(line[1])
		d_basicInfo[key] = grabBasicInfo(line, "target")
	if useReference:
		for line in rel:
			if not ''.join(line).strip(): continue
			key = os.path.basename(line[1])
			d_basicInfo[key] = grabBasicInfo(line, "reference")
	spInfo  = [["#protein", "gtf", "genome", "cds"]]
	refInfo = [["#protein", "gtf", "genome", "cds"]]
	jobs=[]
	for key in d_basicInfo:
		path_gtfIn       = d_basicInfo[key]["gtf"]
		path_genome      = d_basicInfo[key]["genome"]
		path_cdsFastaOut = path_cdsDir+"/"+key+".cds.fasta"
		path_aaFastaOut  = path_aaDir+"/"+key+".aa.fasta"
		if d_basicInfo[key]["type"] == "target":
			spInfo += [[path_aaFastaOut, path_gtfIn, path_genome, path_cdsFastaOut]]
		else:
			refInfo += [[path_aaFastaOut, path_gtfIn, path_genome, path_cdsFastaOut]]
		jobs.append([prepareSpecies, (path_gtfIn, path_genome, path_cdsFastaOut, path_aaFastaOut)])
	runJobs(jobs, int_cores)
	writeCsv(spInfo, path_speciesInfoFile)
	writeCsv(refInfo, path_refInfoFile)
	# Run orthofinder on the newly extracted proteomes
	for path_file in glob.glob(path_aaDir + "/Results*"):
		deleteIfPresent(path_file)
	runOrthoFinder(path_aaDir, int_cores)
	# Grab the output files from orthofinder
	path_orthoFinderOutputFile	= getOrthogroupsFile(path_aaDir)
	path_singletonsFile		= getSingletonsFile(path_aaDir)
	if not useReference: path_refInfoFile == ""
	return path_speciesInfoFile, path_orthoFinderOutputFile, path_singletonsFile, path_refInfoFile

def prepareSpecies(path_gtfIn, path_genome, path_cdsFastaOut, path_aaFastaOut):
	print("Extracting sequences from " + path_gtfIn)
	checkChromosomes(path_gtfIn, path_genome)
	fetchSequences(path_gtfIn, path_genome, path_cdsFastaOut, path_aaFastaOut, 1)
	print("Finished extracting sequences from " + path_gtfIn)

def checkCdsHealth(path_inputGtf, path_outputGtf):
	gtf = readCsv(path_inputGtf)
	transcripts=[]; cds={}; starts={}; stops={}; entries={}; strands={}
	for line in gtf:
		tid = re.sub(r'.*transcript_id \"([^\"]*)\".*', r'\1', line[8])
		transcripts.append(tid)
		if line[2].lower() == "start_codon":
			starts[tid] = starts.get(tid, []) + [[int(line[3]), int(line[4])+1]]
		elif line[2].lower() == "stop_codon":
			stops[tid] = stops.get(tid, []) + [[int(line[3]), int(line[4])+1]]
		elif line[2].lower() == "cds":
			cds[tid] = cds.get(tid, []) + range(int(line[3]), int(line[4])+1)
		strands[tid] = strands.get(tid, []) + [line[6]]
		entries[tid] = entries.get(tid, []) + [line]
	badGenes=[]
	for t_id in set(transcripts):
		if not t_id in cds or len(cds[t_id]) != len(set(cds[t_id])):
			badGenes.append(t_id); continue
		if not t_id in strands or len(set(strands[t_id])) != 1:
			badGenes.append(t_id); continue
		if not t_id in stops or not t_id in starts or len(stops[t_id]) !=1 or len(starts[t_id]) !=1:
			badGenes.append(t_id); continue
		if strands[t_id][0] == "+":
			if starts[t_id][0][0] != min(cds[t_id]):
				badGenes.append(t_id); continue
			if max(cds[t_id]) != stops[t_id][0][1] -1:
				badGenes.append(t_id); continue
		elif strands[t_id][0] == "-":
			if max(cds[t_id]) != starts[t_id][0][1] -1:
				badGenes.append(t_id); continue
			if min(cds[t_id]) != stops[t_id][0][0]:
				badGenes.append(t_id); continue
		else:
			badGenes.append(t_id); continue
	writeCsv(flattenLol([entries[t_id] for t_id in entries if not t_id in badGenes]), path_outputGtf)

def addSpecies(str_species, d_spInfo):
	if not str_species in d_spInfo:
		highestVal = max([ x["number"] for x in d_spInfo.values() ] or [0])
		d_spInfo[str_species] = {"number": highestVal + 1}
	return d_spInfo

def readInputSet(dataset, d_spInfo, tag):
	for line in dataset:
		if not ''.join(line).strip(): continue
		d_spInfo, str_species = readInputIndividual(line, d_spInfo, tag)
	return d_spInfo	

def readInputLocations(path_speciesInfoFile, path_referenceFile):
	"""Read CSV file containing the locations for the sequence files, etc.
	   Each row is a species. The columns are [proteins, gtfs, genome, cds].
	   The proteins string should be the same as the species element in the OrthoFinder output.
	   As such, we use the basename of this file as the species name.
	"""
	print("Loading and checking input data locations...")
	d_spInfo = readInputSet(readCsv(path_speciesInfoFile), {}, "target")
	if not path_referenceFile == "":
		d_spInfo = readInputSet(readCsv(path_referenceFile), d_spInfo, "reference")
	return d_spInfo

def readInputIndividual(line, d_spInfo, tag):
	# The "species", i.e. the basename for the source protein file
	# will be the key in the dictionary.
	str_species = os.path.basename(line[0])
	d_spInfo = addSpecies(str_species, d_spInfo)
	# Then just build up the dictionary with human-readable names
	d_spInfo[str_species]["protein"] = path_aa      = checkFileExists(line[0])
	d_spInfo[str_species]["gtf"]     = path_gtf     = checkFileExists(line[1])
	d_spInfo[str_species]["genome"]  = path_genome  = checkFileExists(line[2])
	d_spInfo[str_species]["cds"]     = path_cds     = checkFileExists(line[3])
	d_spInfo[str_species]["type"]    = tag
	print("Checking chromosomes...")
	checkChromosomes(path_gtf, path_genome)
	checkSequences(path_gtf, path_cds, path_aa)
	return d_spInfo, str_species

###########################################
######## Orthogroup GTF extraction ########
###########################################
def gtfsForOrthoGroups(path_ogGtfDir, path_orthogroups, path_singletons, d_spInfo, int_cores):
	b = readCsv(path_orthogroups); 	bs = readCsv(path_singletons)
	speciesList_og        = b[0]
	speciesList_original  = list(d_spInfo.keys())
	speciesList_useful    = convertOrthoFinderSpeciesNames(speciesList_og, speciesList_original)
	orthogroups = b[1:]; singletons  = bs[1:]
	jobs=[]
	for i in range(1,len(speciesList_useful)):
		str_species=speciesList_useful[i]
		a = readCsv(d_spInfo[str_species]["gtf"])
		print("Extracting orthogroup and singleton gtf files for " + str_species)
		jobs.append([gtfsForGroups, (a, orthogroups, path_ogGtfDir, str_species, "_orthoProtein.gtf", i)])
		jobs.append([gtfsForGroups, (a, singletons, path_ogGtfDir, str_species, "_singletonProtein.gtf", i)])
	runJobs(jobs, int_cores)
	
def gtfsForGroups(list_gtf, orthogroups, path_ogGtfDir, str_species, str_outsuffix, int_speciesNum):
	#gtf entries by transcript name
	aa = [[re.sub(".*transcript_id[ =]\"([^\"]*)\".*", r'\1', x[8]), x] for x in list_gtf]
	c  = dict((x[0], []) for x in aa)
	for x in aa: c[x[0]].append(x[1])
	e = dict((x[0], [c[i] for i in itertools.ifilterfalse(lambda x: x=='', re.split("[ ,]*", x[int_speciesNum]))]) for x in orthogroups)
	for orthogroup in e:
		filename = path_ogGtfDir + "/" + orthogroup+"." + str_species + str_outsuffix
		with open(filename, 'w') as mycsvfile:
			datawriter = csv.writer(mycsvfile, delimiter = '\t',quoting = csv.QUOTE_NONE, quotechar='')
			for row in list(itertools.chain.from_iterable(e[orthogroup])):
				datawriter.writerow(row + [str_species, orthogroup])

######################################
############ Alignments ##############
######################################

def implementGetProteinFastaFiles(orthogroup, orthogroupProteinSequences, path_ogAlDir):
        # Define output files
	path_protSeqFile = path_ogAlDir + "/" + orthogroup + "_ProteinSequences.fasta"
	# Write out the sequences
	writeSequencesToFastaFile(orthogroupProteinSequences, path_protSeqFile)

def implementGetProteinAlignments(path_proteinFastaFile, path_fastaOut):
	alignment = MafftCommandline(input=path_proteinFastaFile, auto="on")
        stdout, stderr = alignment()
	write(stdout, path_fastaOut)

def getProteinFastaFiles(orthogroups, proteinSequences, d_sequenceInfoById, d_spInfo, path_ogAlDir, int_cores):
	jobs={}
	for orthogroup in orthogroups:
		orthogroupProteinSequences = [proteinSequences[x] for x in orthogroups[orthogroup]]
		jobs[orthogroup] = [implementGetProteinFastaFiles, (orthogroup, orthogroupProteinSequences, path_ogAlDir)]
	runAndTrackJobs(jobs, int_cores)

def getNucleotideAlignment(alignedProteinsFastaIn, fastaOut, d_cds):
	"""Get a nucleotide alignment from the protein alignment.
        """
	aa = SeqIO.parse(alignedProteinsFastaIn, "fasta")
	nuc = []
	for sequence in aa:
		# Grab the cds sequence
		cds = d_cds[sequence.id]; dnaSourceSequence = str(cds.seq);
		# Grab the gapped aa sequence
		gappedProteinSequence = str(sequence.seq)
		# Run one through the other
		gDnaSeq = threadGappedProteinSequenceThroughDNA(gappedProteinSequence, dnaSourceSequence)
		# Construct the new alignment
		gDnaId = cds.id; gDnaName = cds.name; gDnaDesc = cds.description
		nuc.append(SeqRecord(Bio.Seq.Seq(gDnaSeq), id=gDnaId, name=gDnaName, description=gDnaDesc))
	SeqIO.write(nuc, fastaOut, "fasta")

def getProteinAlignments(orthogroups, path_ogAlDir, int_cores):
    jobs={};
    for orthogroup in orthogroups:
        path_prots      = path_ogAlDir + "/" + orthogroup + "_ProteinSequences.fasta"
        path_alignment	= path_ogAlDir + "/" + orthogroup + "_ProteinAlignment.fasta"
        jobs[orthogroup]= [implementGetProteinAlignments, (path_prots, path_alignment)]
    runAndTrackJobs(jobs, int_cores)

def getNucleotideAlignments(orthogroups, path_ogAlDir, d_sequenceInfoById, d_spInfo, int_cores):
        jobs={}; d_indexedCds = {}
	for str_species in d_spInfo:
                d_indexedCds[str_species] = SeqIO.index(d_spInfo[str_species]["cds"], "fasta")
	l = len(orthogroups)
	for i, orthogroup in enumerate(orthogroups):
		path_proteinAlignmentFile = path_ogAlDir + "/" + orthogroup + "_ProteinAlignment.fasta"
		path_nucAlignmentFile	  = path_ogAlDir + "/" + orthogroup + "_NucAlignment.fasta"
		orthogroupProteinSequences = [x for x in orthogroups[orthogroup]]
		d_cds = {}
		for o in orthogroupProteinSequences:
			seqid    = d_sequenceInfoById[o].seqId
			species  = d_sequenceInfoById[o].species
			cds      = d_indexedCds[species][seqid]
			d_cds[o] = cds
		jobs[orthogroup] = [getNucleotideAlignment, (path_proteinAlignmentFile, path_nucAlignmentFile, d_cds)]
		sys.stdout.write("\r    Preparing "+str(i)+" of "+str(l)+" orthogroups...")
		sys.stdout.flush()
	sys.stdout.write("\r    Finished preparing orthogroups            ")
	sys.stdout.flush()
	runAndTrackJobs(jobs, int_cores)

def threadGappedProteinSequenceThroughDNA(gappedProteinSequence, dnaSourceSequence):
	"""Protein must correspond exactly to DNA sequence.
	   Both sequences should be provided as strings.
	"""
	gappedDnaSequence = ""; int_counter = 0
	for ch in gappedProteinSequence:
		if ch =="-":
			gappedDnaSequence = gappedDnaSequence + "---"
		else:
			codon=dnaSourceSequence[int_counter:int_counter+3]
			if not len(codon) == 3:	codon += "-"*(-len(codon) % 3)
			gappedDnaSequence = gappedDnaSequence + codon
			int_counter = int_counter + 3
	return gappedDnaSequence

def fetchSequences(path_gtfIn, path_genome, path_cdsFastaOut, path_aaFastaOut, int_translationTable):
	path_nucleotideSequences = path_cdsFastaOut
	#print("fetching nucleotide sequences for " + path_gtfIn)
	callFunction("infile=" +  path_gtfIn+ "; outfile=" + path_nucleotideSequences + "; genome=" + path_genome + """;
		tf=`mktemp -d`
		gtfCds="$tf/gtfCds"
		gtfBed="$tf/gtfBed"
		
		#Prepare the gtf
		#echo "preparing gtf..."
		grep -vP "^$" $infile | awk '$3=="CDS"' > $gtfCds
		cut -f1-8 $gtfCds > $gtfBed.1
		sed -r "s/.*transcript_id[ =]\\"?([^\\";]*)\\"?;?.*/\\1/g" $gtfCds > $gtfBed.2
		paste $gtfBed.1 $gtfBed.2 | perl -ne 'chomp; @l=split; printf "$l[0]\\t%s\\t$l[4]\\t$l[8]\\t.\\t$l[6]\\n", $l[3]-1' | sort -u | sort -k1,1V -k2,2n > $gtfBed

		#Negative strand
		#echo "negative strand..."
		awk '$6=="-"' $gtfBed > $gtfBed.neg
		bedtools getfasta -name -s -fullHeader -fi $genome -fo $gtfBed.neg.tab -bed $gtfBed.neg -tab
		tac $gtfBed.neg.tab | awk '{a[$1]=a[$1]""$2} END {for (i in a) {print ">"i"\\n"a[i]}}' > $gtfBed.neg.fa

		#Then positive strand
		#echo "positive strand..."
		awk '$6=="+"' $gtfBed > $gtfBed.pos
		bedtools getfasta -name -s -fullHeader -fi $genome -fo $gtfBed.pos.tab -bed $gtfBed.pos -tab
		cat $gtfBed.pos.tab | awk '{a[$1]=a[$1]""$2} END {for (i in a) {print ">"i"\\n"a[i]}}' > $gtfBed.pos.fa

		cat $gtfBed.pos.fa $gtfBed.neg.fa | sed -r "s/^>(.*)$/£££>\\1###/g" | sed -r \"s/$/###/g\" | tr '\\n' ' ' | sed -r "s/£££/\\n/g" | sed -r "s/### //g" | grep -v XXX | grep -v "\*[A-Z]" | grep -v "###$" | sed -r "s/###/\\n/g" | grep -vP "^$" > $outfile

		#echo $tf
		rm -r $tf
		""")
	#print("translating to protein...")
	sequences=SeqIO.parse(path_nucleotideSequences, "fasta")
	protSequences=[]
	for s in sequences:
		s_p = s
		s_p.seq = s.seq.translate(table=int_translationTable)
		protSequences.append(s_p)
	SeqIO.write(protSequences, path_aaFastaOut, "fasta")

def getProteinSequences(sequencesHolder, d_spInfoDict):
	"""Get the protein sequences for each feature per species in an orthogroup
	   d_spInfo should be the dictionary version of the file locations, indexed by species
	"""
	# To make protein sequence access faster, index.
	d_indexedProteinFiles = {}
	for str_species in d_spInfoDict:
		d_indexedProteinFiles[str_species] = SeqIO.index(d_spInfoDict[str_species]["protein"], "fasta")
	d_proteinSequencesHolder = {}
	for sequence in sequencesHolder:
		seqRef = sequencesHolder[sequence]
		proteinSequence = d_indexedProteinFiles[seqRef.species][seqRef.seqId]
		# Replace the species-local id with the uniqueId
		proteinSequence.id = seqRef.uniqueId
		d_proteinSequencesHolder[seqRef.uniqueId] = proteinSequence
	# Returns a dictionary of id / SeqIO sequence object.
	return d_proteinSequencesHolder

def writeSequencesToFastaFile(proteinSequences, path_outputFile):
	"""Write out sequences into fasta file. Overwrites file if necessary.
	"""
	SeqIO.write(proteinSequences, path_outputFile, "fasta")

######################################
############ OrthFinder ##############
######################################

def rerunOrthoFinder(path_wDir, d_spInfo):
	path_newProteomesDir = makeIfAbsent(path_wDir + "/newProteomes", fresh = True)
	d_spInfo_modern = {}
	for str_species in d_spInfo:
		d_spInfo[str_species]["newProtein"]     = path_newProteome = path_newProteomesDir + "/" + str_species + "newProteome.fasta"
		d_spInfo[str_species]["newSpeciesName"] = str_species + "newProteome.fasta"
		path_oldProteome = d_spInfo[str_species]["protein"]
		if d_spInfo[str_species]["type"] == "target":
			path_predictedProteinSequences = d_spInfo[str_species]["augustussequences_hintfiltered"]
			makeNewProteome(path_oldProteome, path_predictedProteinSequences, path_newProteome)
		else:
			callFunction("cp " + path_oldProteome + " " + path_newProteome)
		d_spInfo_modern[str_species + "newProteome.fasta"] = {"number": d_spInfo[str_species]["number"]}
	runOrthoFinder(path_newProteomesDir, int_cores)
	return d_spInfo_modern, path_newProteomesDir

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
			exit("Error: Can't find orthofinder. Looked for orthofinder in the following order: OrthoFiller.py directory, execution directory, system PATH. Please ensure orthofinder is either installed and included in your PATH or that the orthofinder.py file is included in the same directory as the OrthoFiller.py file. Orthofinder can be downloaded from https://github.com/davidemms/OrthoFinder")
	# Call orthofinder
	version=tuple([int(i) for i in commands.getstatusoutput(finderString + " -h | grep -P \"version +[0-9\.]+\" | sed -r \"s/.*version ([0-9.]+).*/\\1/g\"")[1].split(".")])
	basicCall = finderString + " -f " + path_aaDir + " -t " + str(int_cores)
	ending = "" if version < (1,0,2) else (" -g" if version < (1,1,2) else " -og")
	callFunction(basicCall + ending)

def readOrthoFinderOutput(path_orthoFinderOutputFile, path_orthoFinderSingletonsFile, d_spInfo):
	"""Read CSV file into 3-tiered dictionary: orthogoup > species > sequences.
	   The last entry is an array of sequence strings.
	   The species names need a little extra attention, since some versions of orthofinder change the names
		of the input files to auto-generate a species name.
	   Output is a dictionary, keyed by unique id valued by a seqRef object, and a list of the orthogroups, and of the singletons.
	"""
	print("Reading orthogroups from " + path_orthoFinderOutputFile)
	d_orthogroups, sequences_orthogroups = readOrthoFinderOutputIndividual(path_orthoFinderOutputFile, d_spInfo)
	print("Reading singletons from " + path_orthoFinderSingletonsFile)
	d_singletons, sequences_singletons = readOrthoFinderOutputIndividual(path_orthoFinderSingletonsFile, d_spInfo)
	sequences_all = merge_dicts(sequences_orthogroups, sequences_singletons)
	return sequences_all, d_orthogroups, d_singletons
	
def readOrthoFinderOutputIndividual(path_orthoFinderOutputFile, d_spInfo):
	d_groups = {}
	sequences_local = {}
	with open(path_orthoFinderOutputFile) as csvfile:
		data = csv.reader(csvfile, delimiter="\t")
		# First line contains the names of the species, which orthoFinder automatically
		# generates from the fasta input names. Need to double check these names against the ones we have.
		# The useful species list is just a copy of the d_spInfo species, just making sure that 
		# it's in the same order as the species listed in the orthofinder output file header.
		speciesList_original = list(d_spInfo.keys())
		speciesList_useful   = convertOrthoFinderSpeciesNames(data.next(), speciesList_original)
		# Each subsequent line has orthogroup as first entry, and grouped sequence IDs
		# for each numbered column.
		for line in data:
			og = line[0]
			for i in range(1,len(speciesList_useful)):
				species    = speciesList_useful[i]
				speciesNum = d_spInfo[species]["number"]
				# Get rid of any empty entries.			
				for sequence in itertools.ifilterfalse(lambda x: x=='', re.split("[,]*", line[i])):
					seqRef       = SeqRef(species, speciesNum, sequence.strip('"').strip(" "))
					d_groups[og] = d_groups.get(og, []) + [seqRef.uniqueId]
					sequences_local[seqRef.uniqueId] = seqRef
	return d_groups, sequences_local
				
def convertOrthoFinderSpeciesNames(speciesList_og, speciesList_original):
	speciesList_useful = [""]
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
			exit("OrthoFinder output contains one or more names that do not corrrespond to the inputted fasta names.\n Please ensure all protein fasta names are unique.\n If you ran OrthoFiller without using the --prep option, make sure you have the latest version of OrthoFinder and try again.")
	return speciesList_useful

def findFile(options, directory, errMsg):
	for f in options:
		a = find(f, directory)
		if a: return a
	exit(errMsg)

def getOrthogroupsFile(path_aaDir):
	return findFile(["OrthologousGroups.csv", "Orthogroups.csv"], path_aaDir, "Can't find Orthogroup output file... exiting...")

def getSingletonsFile(path_aaDir):
	return findFile(["OrthologousGroups_UnassignedGenes.csv", "Orthogroups_UnassignedGenes.csv"], path_aaDir, "Can't find Singletons output file... exiting...")

######################################
############### HMMs #################
######################################

def buildHmm(nucAlignment, path_outputFile):
	"""Build an hmm based on a nucleotide alignment. Inputs are file names.
	"""
	callFunctionMezzoPiano("hmmbuild --informat afa " + path_outputFile + " " + nucAlignment) #qgr

def makeHmmerDb(path_fa, path_dbOutput):
	"""Makes a database per cds file for use with hmmer.
	"""
	callFunctionQuiet("makehmmerdb --block_size=10 " + path_fa + " " + path_dbOutput) #qgr

def implementHmmSearch(path_hmmFile, path_db, path_hitsFile, species, orthogroup):
	"""Runs across the genome and finds hmm hits
	"""
	callFunctionMezzoPiano("nhmmer --tformat hmmerfm --dna --cpu 1 --tblout " + path_hitsFile + " " + 	path_hmmFile + " " + path_db) #qgr
	path_hitsFileBed = path_hitsFile + ".bed"
	makeBed(path_hitsFile, species, orthogroup, path_hitsFileBed, remove=True)

def prepareHmmDbInfo(d_spInfo, path_hmmDbDir, splitByChr):
	"""Separate this out to make it easier to turn off hmmdb creation
	"""
	for str_species in d_spInfo:
		if not d_spInfo[str_species]["type"] == "target": continue
		print("Preparing database directories for " + str_species)
		path_genome = d_spInfo[str_species]["genome"]
		#Working by chromosome reduces memory consumption.
		speciesTag=path_hmmDbDir+"/"+str_species
		if splitByChr:
			makeIfAbsent(path_hmmDbDir+"/"+str_species)
			chromosomes = list(SeqIO.parse(path_genome, "fasta"))
			d_spInfo[str_species]["hmmdb"]={}
			for chrInfo in chromosomes:
				path_chr = speciesTag + "/"+ str_species + "."+chrInfo.id+".fa"
				path_db  = path_chr + ".hmmdb"
				d_spInfo[str_species]["hmmdb"][chrInfo.id] = [path_chr, path_db]
				SeqIO.write(chrInfo, path_chr, "fasta")
		else:
			d_spInfo[str_species]["hmmdb"] = {"full": [path_genome, speciesTag+".hmmdb"]}
		
def prepareHmmDbs(d_spInfo, path_hmmDbDir, int_cores):
	jobs = []
	for str_species in d_spInfo:
		if not d_spInfo[str_species]["type"] == "target": continue
		print("Preparing databases for " + str_species)
		#Working by chromosome reduces memory consumption.
		for i in d_spInfo[str_species]["hmmdb"]:
			path_fa, path_db = d_spInfo[str_species]["hmmdb"][i]
			jobs.append([makeHmmerDb, (path_fa, path_db)])
	runJobs(jobs, int_cores)

def buildHmms(orthogroups, path_ogAlDir, path_ogHmmDir, int_cores):
	jobs={};
	for orthogroup in orthogroups:
		path_nucAlignmentFile  = path_ogAlDir + "/" + orthogroup + "_NucAlignment.fasta"
		path_hmmFile           = path_ogHmmDir + "/" + orthogroup + ".hmm"
		jobs[orthogroup]       = [buildHmm, (path_nucAlignmentFile, path_hmmFile)]
	runAndTrackJobs(jobs, int_cores)

def prepareHitDirs(orthogroups, d_spInfo, path_ogHmmDir, path_ogHitsDir, splitByChr):
	print("Preparing HMMs and HMM directories...")
	# Split the genomes up by chromosome if necessary, and run the hmms.
	for species in d_spInfo:
		if not d_spInfo[species]["type"] == "target": continue
		path_spHitsDir = makeIfAbsent(path_ogHitsDir + "/" + species)
		d_spInfo[species]["hitDirs"] = []
		if splitByChr:
			tChr=len(d_spInfo[species]["hmmdb"])
			for nChr, chromosome in enumerate(d_spInfo[species]["hmmdb"]):
				path_hmmdb =d_spInfo[species]["hmmdb"][chromosome][1]
				str_nChr   = str(nChr+1)
				path_chrHitDir = makeIfAbsent(path_spHitsDir + "/chrscf_" + str_nChr)
				d_spInfo[species]["hitDirs"] += [{"id": chromosome, "num": str_nChr, "hitDir": path_chrHitDir, "hmmdb": path_hmmdb}]
		else:
			path_hmmdb=d_spInfo[species]["hmmdb"]["full"][1]
			d_spInfo[species]["hitDirs"] = [{"id": "full", "num": 1, "hitDir": path_spHitsDir, "hmmdb": path_hmmdb}]

def runHmms(orthogroups, d_spInfo, path_ogHmmDir, path_ogHitsDir, int_cores, splitByChr):
	#Sort hmms by size so that we clear the biggest ones first.
	weights=dict((o, os.path.getsize(path_ogHmmDir + "/" + o + ".hmm")) for o in orthogroups)
	orthogroups_sorted = sorted(weights.keys(), key=lambda x: int(weights[x]), reverse=True)
	# Split the genomes up by chromosome if necessary, and run the hmms.
	for species in d_spInfo:
		if not d_spInfo[species]["type"] == "target": continue
		print("Running HMMs on species " + species)
		jobs={}
		if splitByChr:
			tChr = len(d_spInfo[species]["hmmdb"])
			for h in d_spInfo[species]["hitDirs"]:
				path_chrHitDir = h["hitDir"]
				path_hmmdb     = h["hmmdb"]
				str_nChr       = h["num"]
				chromosome     = h["id"]
				for orthogroup in orthogroups_sorted:
					path_hmmFile  = path_ogHmmDir + "/" + orthogroup + ".hmm"
					path_hitsFile = path_chrHitDir + "/" + orthogroup + "." + species + "." + chromosome + ".hits"
					jobs[orthogroup] = [implementHmmSearch, (path_hmmFile, path_hmmdb, path_hitsFile, species, orthogroup)]
				runAndTrackJobs(jobs, int_cores, "Chromosome/scaffold " + str_nChr + " of " + str(tChr) + ": ", False)
			print("\n")
		else:
			path_hitDir = d_spInfo[species]["hitDirs"][0]["hitDir"]
			path_hmmdb = d_spInfo[species]["hitDirs"][0]["hmmdb"]
			for orthogroup in orthogroups_sorted:
				path_hmmFile = path_ogHmmDir + "/" + orthogroup + ".hmm"
				path_hitsFile = path_hitDir + "/" + orthogroup + "." + species + ".hits"
				jobs[orthogroup] = [implementHmmSearch, (path_hmmFile, path_hmmdb, path_hitsFile, species, orthogroup)]
			runAndTrackJobs(jobs, int_cores, "", True)
			print("")

def processHmmOutput(d_spInfo, path_candidates, path_ogHitsDir, path_ogGtfDir, splitByChr):
	d_ogIntersectionFileNamesAnnotated = {}	
	jobs = []
	for str_species in d_spInfo:
		if not d_spInfo[str_species]["type"] == "target": continue
		print("Submitting HMM output files for species " + str_species + "...")
		d_ogIntersectionFileNamesAnnotated[str_species] = path_hitsOgIntersectionFileNameAnnotated \
						= path_candidates + "/" + str_species + ".hitsIntersectionOrthogroups.annotated.bed"
		jobs.append([processHmmOutputIndividual, (str_species, \
						path_candidates, \
						path_ogHitsDir, \
						path_ogGtfDir, \
						path_hitsOgIntersectionFileNameAnnotated, \
						d_spInfo, \
						splitByChr)])
	runJobs(jobs, int_cores)
	return d_ogIntersectionFileNamesAnnotated

def processHmmOutputIndividual(str_species, path_wDir, path_ogHitsDir, path_ogGtfDir, path_hitsOgIntersectionFileNameAnnotated, d_spInfo, splitByChr):
	path_hitsBedFileName    = path_wDir + "/" + str_species + ".allHits.bed"
	path_ogBedFileName      = path_wDir + "/" + str_species + ".allOrthogroups.bed"
	path_hitsOgIntersectionFileName	= path_wDir + "/" + str_species + ".hitsIntersectOrthogroups.bed"
	# Get all hits together, and double check that they are good files.
	# The hits directories should be neatly listed, regardless of whether they were calculated splitwise.
	for i, hDir in enumerate(d_spInfo[str_species]["hitDirs"]):
		searchString = hDir["hitDir"] + "/*bed"
		concatGlob(searchString, path_hitsBedFileName + ".chr." + str(i)+".tmp")
	# Get all hits together, and all protein annotations.
	concatGlob(path_hitsBedFileName + ".chr.*.tmp", path_hitsBedFileName, removeOrig = True)
	concatGlob(path_ogGtfDir + "/*" + str_species + "*Protein.gtf", path_ogBedFileName)
	# Clean out the rubbish
	bashIt("awk '$3 > 0 && $2 > 0'", path_hitsBedFileName)
	bashIt("awk '$3==\"CDS\"' ", path_ogBedFileName)
	# Intersect
	callFunction("bedtools intersect -nonamecheck -a " + path_hitsBedFileName + " -b " + path_ogBedFileName + " -wa -wb > " + path_hitsOgIntersectionFileName)
	# Begin annotation
	callFunction("cat " + path_hitsOgIntersectionFileName + " " + path_hitsBedFileName + " | cut -f1-11 | sort | uniq -u | sed -r \"s/$/\\t.\\t.\\t.\\t.\\t.\\t.\\t.\\t.\\t.\\t.\\t./g\" > " + path_hitsOgIntersectionFileName + ".tmp; cat " + path_hitsOgIntersectionFileName + ".tmp " + path_hitsOgIntersectionFileName + " > " + path_hitsOgIntersectionFileName + ".tmp.tmp ; mv " + path_hitsOgIntersectionFileName + ".tmp.tmp " + path_hitsOgIntersectionFileName + "; rm " + path_hitsOgIntersectionFileName + ".tmp")
	annotateIntersectedOutput(path_hitsOgIntersectionFileName, path_hitsOgIntersectionFileNameAnnotated)

def annotateIntersectedOutput(path_hitsOgIntersectionFileName, path_hitsOgIntersectionFileNameAnnotated):
	callFunction("infile=\""+ path_hitsOgIntersectionFileName+"\"; outfile=\""+path_hitsOgIntersectionFileNameAnnotated+"\";\
		awk -F \"\\t\" '$22 == \".\"' $infile | sed -r \"s/_id[= ]*[^\\\"]*\\\"/_id=\\\"/g\" | sed -r \"s/$/\\tmatch_none/g\" > $outfile ;\
		awk -F \"\\t\" '$22 != \".\"' $infile | awk -F \"\\t\" '$22 != $11' | sed -r \"s/_id[ =]*[^\\\"]*\\\"/_id=\\\"/g\" | sed -r \"s/$/\\tmatch_bad/g\" >> $outfile;\
		awk -F \"\\t\" '$22 != \".\"' $infile | awk -F \"\\t\" '$22 == $11' | sed -r \"s/_id[ =]*[^\\\"]*\\\"/_id=\\\"/g\" | sed -r \"s/$/\\tmatch_good/g\" >> $outfile;")

def getBases(path_gtf, path_gtfBases):
	"""Get the base for a gtf file
	"""	
	gtf = readCsv(path_gtf); coords  = {}; 	entries = {}
	for line in gtf:
		if not line[2].lower() == "cds": continue
		transcript_id = re.sub(r'.*transcript_id \"([^\"]*)\".*', r'\1', line[8])
		if not transcript_id in coords:
			coords[transcript_id] = []
			entries[transcript_id] = line
		coords[transcript_id] += [int(line[3]), int(line[4])]
	for t_id in entries:
		entries[t_id][3] = min(coords[t_id])
		entries[t_id][4] = max(coords[t_id])
	writeCsv(entries.values(), path_gtfBases)

def makeBed(path_hitsFile, species, orthogroup, path_hitsFileBed, remove = False):
	callFunction("grep -v \"#\" " + path_hitsFile + " | sed -r \"s/ +/\t/g\" | cut -f1,7,8,12,13,14,15 | perl -ne \
                                'chomp;@l=split; printf \"%s\\t%s\\t%s\\t.\\t.\\t%s\\t" + species + "\\t" + orthogroup +
                                "\\n\", $l[0], ($l[1] + $l[2] - abs($l[1] - $l[2])) / 2, ($l[1] + $l[2] + abs($l[2] - $l[1])) / 2,\
                                 join(\"\\t\", @l[3..6])' > " + path_hitsFileBed)
	if remove: os.remove(path_hitsFile)

def extractFromFastaByName(path_gtfFile, path_fastaFile, path_fastaOut):
	with open(path_gtfFile, "r") as f:
		data = list(csv.reader((row for row in f if not row.startswith('#')), delimiter="\t"))
		#tids = [re.sub(r".*transcript_id[ =]\"([^\"]*)\".*", r"\1", i[8]) for i in data]
		tids = [i[8] for i in data]
	sequences = [a for a in SeqIO.parse(path_fastaFile, "fasta") if a.description in tids]
	SeqIO.write(sequences, path_fastaOut, "fasta")

########################################################
###################### Augustus ########################
########################################################

def trainAugustus(d_spInfo, path_wDir, pool):
	"""trainAugustus - Trains augustus using the genomes of the input species
	"""
	path_augWDir = makeIfAbsent(path_wDir + "/augustus")
	#Python doesn't like us to edit a dictionary while iterating over it.
	for str_species in d_spInfo.keys()[:]:
		path_genome          = d_spInfo[str_species]["genome"]
		path_gtfForTraining  = d_spInfo[str_species]["gtfForTraining"]
		path_augustusSpecies = d_spInfo[str_species]["augustusSpecies"]
		path_augSpeciesWDir  = makeIfAbsent(path_augWDir + "/" + str_species)
		if d_spInfo[str_species]["needsTraining"]:
			print("training augustus on " + str_species)
			async(pool, trainAugustusIndividual, args=(path_augustusSpecies, path_genome, path_gtfForTraining, path_augSpeciesWDir))

def makeGtfTrainingFile(path_inputGtf, path_outputGtf):
	"""For AUGUSTUS to train its algorithms correctly, we need to format
	   the gtf file in a certain way.
	"""
	#print("making training file " + path_outputGtf + " from  " + path_inputGtf + "...")
	#print("making training file " + path_outputGtf + "...")	
	path_tmp=path_outputGtf + ".tmp"
	path_bases=path_outputGtf + ".bases"
	getBases(path_inputGtf, path_bases)
	# Make sure there are no overlaps, by randomly choosing between overlapping entries, and sort the file.
	function="infile=\"" + path_inputGtf + "\"; basefile=\"" + path_bases + "\"; outfile=\"" + path_tmp + "\"; " + """
		td=`mktemp -d`
		rm -f $outfile

		#echo "Assuring transcripts..."
		infile_td=`mktemp $td/infile_tid.XXXXXX`
		sed -r  '/transcript_id/! s/gene_id([ =])\\"([^\\"]*)\\";?( ?)/gene_id\\1\\"\\2\\"; transcript_id\\1\\"\\2.t999\\";\\3/g' $infile > $infile_td

		#echo "Grouping into regions.."
		sort -k1,1V -k4,4n $basefile | bedtools merge -s -i - > $td/gtfmerged.bed.tmp

		cut -f1,2,3,4 $td/gtfmerged.bed.tmp | sed -r "s/\\t([^\\t]*)$/\\t.\\t.\\t\\1/g" > $td/gtfmerged.bed

		#echo "Intersecting..."
		bedtools intersect -a $td/gtfmerged.bed -b $infile_td -wa -wb > $td/gtfis.bed
		
		cat $td/gtfis.bed | shuf | sed -r  "s/(.*transcript_id[ =]\\")([^\\"]*)(\\".*)/\\2\\t\\1\\2\\3\\t\\2/g" | awk 'BEGIN {FS="\\t"} {if (a[$2"."$3"."$4"."$7] == "") { a[$2"."$3"."$4"."$7]=$1 } ; if (a[$2"."$3"."$4"."$7]==$1) {v[$2"."$3"."$4"."$7]=v[$2"."$3"."$4"."$7]"\\n"$0"\\t"$2"."$3"."$4"."$7 } } END { for ( i in a ) {print v[i] } } ' | awk 'NF' | cut -f8- | sed -r "s/.\tgene_id/.\tgene_id/g" | sed -r "s/\.\-/\.neg/g" | sed -r "s/\.\+/\.pos/g" > $td/tmp1
		awk -F "\\t" '{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t.\\t"$7"\\t.\\tgene_id \\""$11".gene\\"; transcript_id \\""$11".gene.t1\\";\t"$11}' $td/tmp1 | sort -u > $outfile
		
		rm $basefile
		rm -r $td"""
	callFunction(function)
	print("check st(art|op) codon consistency")
	# Check each gene has a start codon and a stop codon and that they're in the right place
	checkCdsHealth(path_tmp, path_outputGtf)
	os.remove(path_tmp)
	# Add exons as well as CDS
	function="infile=\"" + path_outputGtf + "\"; tmpfile=`mktemp`; tmpfile2=`mktemp`; grep -P \"\\tCDS\\t\" $infile | sed -r \"s/\\tCDS\\t/\\texon\\t/g\" | sed -r \"s/\\t[^\\t]*\\tgene_id/\\t\\.\\tgene_id/g\" > $tmpfile2; cat $infile $tmpfile2 | sort -u | sort -k1,1V -k4,4n > $tmpfile; mv $tmpfile $infile; rm $tmpfile2"
	callFunction(function)

def trainAugustusIndividual(str_augustusSpecies, path_genome, path_gtf, path_augSpeciesWDir):
	callFunctionQuiet("autoAugTrain.pl --useexisting -v -v -v --workingdir=" + \
		path_augSpeciesWDir + " --species=" + str_augustusSpecies + \
		" --trainingset=" + path_gtf + " --genome=" + path_genome +"; wait")

def runAndParseAugustus(path_goodHits, path_genome, path_augustusOut, path_augustusParsedOut, path_fastaOut, path_augustusSpeciesName, path_hintsFile, path_sourceGtf):
	#print("augustus out is " + path_augustusOut)
	#print("augustus parsed out is " + path_augustusParsedOut)
	runAugustus(path_goodHits, path_genome, path_augustusOut, path_augustusSpeciesName, path_hintsFile) #ql
	parseAugustusOutput(path_augustusOut, path_augustusParsedOut, path_fastaOut, path_sourceGtf)

def makeAugustusSpeciesName(species):
	return species+ ".orthofiller." + datetime.datetime.now().strftime("%y%m%d") + "." + \
				''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(9))

def goAugustus(d_spInfo, path_candidates, int_cores):
	jobs = []
	for str_species in d_spInfo:
		if not d_spInfo[str_species]["type"] == "target": continue
		path_proposedGenes = d_spInfo[str_species]["proposedgenes"]
		path_genome        = d_spInfo[str_species]["genome"]
		path_sourceGtf     = d_spInfo[str_species]["gtf"]
		str_augustusOutNameStub = path_candidates + "/" + str_species + ".proposedGenes"
		d_spInfo[str_species]["augustusoutput"]    = path_augustusOut       = str_augustusOutNameStub + ".AugustusModels.gtf"
		d_spInfo[str_species]["augustussequences"] = path_fastaOut          = str_augustusOutNameStub + ".AugustusParsed.sequences.fasta"
		d_spInfo[str_species]["augustusparsed"]    = path_augustusParsedOut = str_augustusOutNameStub + ".AugustusParsed.gtf"	
		d_spInfo[str_species]["hints"]             = path_hintsFile         = str_augustusOutNameStub + ".hints.gtf"
		print("Running Augustus on " + str_species)
		path_augustusSpeciesName = d_spInfo[str_species]["augustusSpecies"]
		jobs.append([runAndParseAugustus, (path_proposedGenes, path_genome, path_augustusOut, path_augustusParsedOut, path_fastaOut, path_augustusSpeciesName, path_hintsFile, path_sourceGtf)])
	runJobs(jobs, int_cores)

def runAugustus(path_goodHits, path_genome, path_augustusOut, path_augustusSpeciesName, path_hitsHintsGtf):
	print("making hints file....")	
	makeHintsFile(path_goodHits, path_hitsHintsGtf)
	print("running augustus, hints file: " + path_hitsHintsGtf + "; augustus species: " + path_augustusSpeciesName + "; genome: " + path_genome)
	callFunctionQuiet("augustus --genemodel=complete --hintsfile=" + path_hitsHintsGtf + \
			" --species=" + path_augustusSpeciesName + " " + path_genome + " > " + path_augustusOut)

def parseAugustusOutput(path_augustusOutput, path_outputGtf, path_outputFasta, path_sourceGtf):
	function = "infile=\"" + path_augustusOutput + "\"; outfile=\"" +\
			path_outputGtf + "\"; fastaout=\"" + path_outputFasta + "\"; sourceGtf=\"" + path_sourceGtf + "\";"+ \
			"""ot=`mktemp -d`; mkdir $ot/augsplit;
			#echo "parsing in $ot";
			awk -v RS="# start gene" -v ot="$ot" '{print "#"$0 > ot"/augsplit/augSplit."NR; close(ot"/augsplit/augSplit."NR) }' $infile
			mkdir $ot/success
			echo "" > $fastaout
			echo "" >  $outfile
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

			bedtools intersect -a $outfile.tmp.genes -b $sourceGtf -wa -wb | cut -f1-9 > $outfile.tmp.genes.remove
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

########################################################
########### Orthogroup test / final output #############
########################################################

def orthogroupTest(d_spInfo, str_species, silNew, oNew, sNew, orthogroups, d_sequenceInfoById):
	path_augustusParsed = d_spInfo[str_species]["augustusparsed_hintfiltered"]
	path_augustusParsedUniq = path_augustusParsed + ".uniq"
	callFunction("grep transcript_id " + path_augustusParsed + " | sed -r 's/.*transcript_id/transcript_id/g' | sort -u > " + path_augustusParsedUniq)
	str_newSpeciesName=d_spInfo[str_species]["newSpeciesName"]
	acceptedSequences=[]
	potentialSequences={}
	with open(path_augustusParsedUniq, "r") as csvfile:
		data = csv.reader(csvfile, delimiter="\t")
		for entry in data:
			sequenceName=re.sub(";g", "; g", re.sub("_id=", "_id ", entry[0]))
			potentialSequences[sequenceName] = re.split("[, ]+", re.sub("possibleOrthos=", "", entry[1]))
	for seqId in potentialSequences:
		seqIdAlt=seqId.replace("_id ", "_id=").replace(" ", "")
		uniqueId=[x for x in silNew if (silNew[x].species==str_newSpeciesName and compareOutputSequences(silNew[x].seqId, seqId))][0]
		list_newOrthogroup = [x for x in oNew if uniqueId in oNew[x]]
		if not list_newOrthogroup: continue
		newOrthogroup=list_newOrthogroup[0]
		oldOrthogroupsPotential=potentialSequences[seqId]
		success = 0
		for oldOrthogroup in oldOrthogroupsPotential:
			overlap = 0
			for oldOrthogroupSequenceId in orthogroups[oldOrthogroup]:
				oldOrthSeq = d_sequenceInfoById[oldOrthogroupSequenceId]
				correspondingSpecies = d_spInfo[oldOrthSeq.species]["newSpeciesName"]
				for newOrthogroupSequenceId in oNew[newOrthogroup]:
					newOrthSeq = silNew[newOrthogroupSequenceId]
					if newOrthSeq.species == correspondingSpecies and compareOutputSequences(newOrthSeq.seqId, oldOrthSeq.seqId):
						overlap = overlap + 1
			if overlap > 0:
				success = 1
				break
		if success == 1:
			acceptedSequences.append(seqId.replace(" ", "").replace("\"", ""))
	return acceptedSequences

def goOrthoGroupTest(path_newProteomesDir, d_spInfo, d_spInfo_modern, orthogroups, int_cores, d_sequenceInfoById):
	path_orthofinderOutputNew	= getOrthogroupsFile(path_newProteomesDir)
	path_orthofinderSingletonsNew	= getSingletonsFile(path_newProteomesDir)
	silNew, oNew, sNew = readOrthoFinderOutput(path_orthofinderOutputNew, \
									path_orthofinderSingletonsNew, d_spInfo_modern)
	jobs={}
	for str_species in d_spInfo:
		if not d_spInfo[str_species]["type"] == "target": continue
		print("Double-checking membership for species " + str_species)
		jobs[str_species] = [orthogroupTest, (d_spInfo, str_species, silNew, oNew, sNew, orthogroups, d_sequenceInfoById)]
	jobres = runAndTrackJobs(jobs, int_cores, ret = True, silent = True)
	return jobres

def assignNames(str_species, path_acceptedGtf, path_geneNameConversionTable, protSequencesAccepted, d_sequenceInfoById, path_augustusParsed, path_acceptedSequencesOut):
	originalNames = [d_sequenceInfoById[x].seqId for x in d_sequenceInfoById if d_sequenceInfoById[x].species == str_species]
	originalNamesStubs = [x.split(".")[0] for x in originalNames ]
	allNames= originalNames + originalNamesStubs
	# for each new gene, give it a nice name and check that it hasn't been used before.
	counter=1
	#print("updating names....")
	with open(path_geneNameConversionTable, "w") as f:
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
	SeqIO.write(protSequencesAccepted, path_acceptedSequencesOut, "fasta")
	blankFile(path_acceptedGtf)
	#print("writing out results....")
	callFunction("while read line ; do \
			sourceId=`echo \"$line\" | cut -f1`; \
			replacementId=`echo \"$line\" | cut -f2 | sed -r \"s/ //g\" | sed -r \"s/[^a-zA-Z0-9._]//g\"`;\
			tid=`echo $sourceId | sed -r \"s/.*transcript_id[= ]\\\"?([^\\\";]*)\\\"?;.*/\\1/g\"`;\
			gid=`echo $sourceId | sed -r \"s/.*gene_id[= ]\\\"?([^\\\";]*)\\\"?;.*/\\1/g\"`;\
			tnum=`echo $tid | sed -r \"s/^$gid\\.//g\"`;\
			grep -P \"(transcript_id[= ]\\\"$tid\\\")|(gene_id[= ]\\\"$gid\\\")|(transcript\\t.*\\t$tid$)|(gene\\t.*\\t$gid$)\" " + \
			path_augustusParsed + " | sed -r \"s/\\t$gid\\t/\\t$replacementId\\t/g\" \
						| sed -r \"s/\\t$tid\\t/\\t$replacementId\\.$tnum\\t/g\" \
						| sed -r \"s/transcript_id[= ]\\\"?$tid\\\"?; ?gene_id[= ]\\\"?$gid\\\"?;/transcript_id \\\"$replacementId\\.$tnum\\\"; gene_id \\\"$replacementId\\\";/g\" | \
								sed -r 's/(.*)\\tpossibleOrthos.*/\\1/g' >> " + path_acceptedGtf +";\
		done < " + path_geneNameConversionTable)

def finishUp(path_resultsDir, d_spInfo):
	print("Results directory is " + path_resultsDir)
	for str_species in d_spInfo:
		if d_spInfo[str_species]["type"] == "target":
			print("-------" + str(len(d_spInfo[str_species]["newGenes"])) + " new genes found for " + str_species)

def goRename(d_spInfo, jobres, path_resultsDir, path_newProteomesDir, d_sequenceInfoById, fullout):
	for str_species in d_spInfo:
		if not d_spInfo[str_species]["type"] == "target": continue
		acceptedSequences = jobres[str_species].get()
		path_acceptedSequencesOut = path_resultsDir + "/" + str_species + ".newSequences.fasta"
		path_newProteome = path_newProteomesDir + "/" + str_species + "newProteome.fasta"
		protSequencesAccepted = [x for x in SeqIO.parse(path_newProteome, "fasta") if x.description.replace(" ", "").replace("\"", "") in acceptedSequences]
		d_spInfo[str_species]["acceptedsequences"] = path_acceptedSequencesOut
		d_spInfo[str_species]["newGenes"] = protSequencesAccepted
		###########################################################
		# Have to make sure gene names are not being duplicated,
		# which can be a problem with iterated runs, for example.
		###########################################################
		# get the list of existing gene names
		print("Reassigning names for " + str_species + "...")
		path_acceptedGtf             = path_resultsDir + "/" + str_species + ".newGenes.gtf"
		path_augustusParsed          = d_spInfo[str_species]["augustusparsed_hintfiltered"]
		path_geneNameConversionTable = path_augustusParsed + ".geneNamesConversion.txt"
		assignNames(str_species, path_acceptedGtf, path_geneNameConversionTable, protSequencesAccepted, d_sequenceInfoById, path_augustusParsed, path_acceptedSequencesOut)
		# write everything out
		if fullout:
			path_resultsFasta = path_resultsDir + "/" + str_species + ".results.aa.fasta"
			path_resultsGtf   = path_resultsDir + "/" + str_species + ".results.gtf"
			d_spInfo[str_species]["resultsGtf"] = path_resultsGtf
			# write out the gtf results file
			blankFile(path_resultsGtf)
			appendFileToFile(d_spInfo[str_species]["gtf"], path_resultsGtf)
			appendFileToFile(path_acceptedGtf, path_resultsGtf)
			# write out the fasta results file
			blankFile(path_resultsFasta)
			appendFileToFile(d_spInfo[str_species]["protein"], path_resultsFasta)
			appendFileToFile(path_acceptedSequencesOut, path_resultsFasta)

	
########################################################
###################### Hit Filter ######################
########################################################

def checkReferenceHits(path_refFile):
	with open(path_refFile, "r") as f:
		data = csv.reader(f, delimiter="\t")
		match_good = 0; match_bad = 0;
		while match_good <= 500 or match_bad <= 500:
			a = next(data, "END")
			if a == "END": break
			if a[-1] == "match_good": match_good += 1
			elif a[-1] == "match_bad": match_bad += 1
		if match_good < 500 or match_bad < 500:
			exit("Error: hit quantity is not sufficient to classify potential new genes.")

def proposeNewGenes(d_spInfo, path_candidates, d_ogIntersectionFileNamesAnnotated, hitFilter, path_allHitsOgIntersectionFileNameAnnotated):
	path_pdir = makeIfAbsent(path_candidates+"/proposed_genes")
	jobs = []
	for str_species in d_spInfo:
		if d_spInfo[str_species]["type"] == "target":
			print "Fitting models to hit data for " + str_species + "..."
			d_spInfo[str_species]["proposedgenes"] = path_outFile = path_pdir + "/" + str_species + ".proposedGenes"
			path_hitsOgIntersectionFileNameAnnotated = d_ogIntersectionFileNamesAnnotated[str_species]
			jobs.append([proposeNewGenesIndividual, (path_hitsOgIntersectionFileNameAnnotated, path_allHitsOgIntersectionFileNameAnnotated, str_species, path_outFile, hitFilter)])
	runJobs(jobs, int_cores)

def unpackFitDistributionScript(path_scriptDestination):
	str_script='library("gamlss")\n\nargs <- commandArgs(TRUE)\n\n\nsourceF=args[1]\naltSourceF= args[2]\noutF=args[3]\n\nprint(paste("reading in source table: ", args[1], sep=""))\na <- read.table(sourceF, sep="\\t", header=FALSE)\n\nnames(a) <- c("hitChr", "hitStart", "hitEnd", "mystery1", "mystery2", "hitStrand", "eVal", "score", "bias", "hitSpecies", "hitOg", "targetChr", "mystery3", "mystery4", "targetStart", "targetEnd", "mystery5", "mystery6", "targetStrand", "geneLabel", "targetSpecies", "targetOg", "match")\n\na <- cbind(a, score_adj=a$score/(a$hitEnd - a$hitStart))\n\nh <- hist(a$score_adj, breaks=50, plot=FALSE)\nb <- h$breaks\n\na_none <- a[a$match=="match_none",]\na_good <- a[a$match=="match_good",]\n\n#Declare variables\ng_good = ""\ng_bad = ""\n\ng_prob_good = ""\ng_prob_bad = "" \n\n# Sampling 1000 data points makes the curve-fitting quicker and hardly affects the fit.\ngetGamlss <- function(theData) { print(theData$t); theData_s <- as.data.frame(sample(theData$t, 1000)); colnames(theData_s) <- c("t"); gamlss(t ~ 1, data=theData_s, family="ST1", method=RS(), gd.tol=10000000, c.cyc=0.001, control=gamlss.control(n.cyc=200)) } \n\nif(nrow(a_good) > 1000) {\n\tprint("source table is good, going ahead...")\n\ta_bad_og <- a[a$match=="match_bad",]\n\ta_bad_singleton <- a[a$match=="match_singleton",]\n\ta_bad <- rbind(a_bad_og, a_bad_singleton)\n\n\tprint("fitting good hits")\n\tgoodScores=as.data.frame(a_good$score_adj); colnames(goodScores) <- c("t")\n\t\n\tg_good <- getGamlss(goodScores)\n\tprint("fitting bad hits")\n\tbadScores=as.data.frame(a_bad$score_adj); colnames(badScores) <- c("t")\n\tg_bad <- getGamlss(badScores)\n\t\n\tg_prob_good =  nrow(a_good) / (nrow(a_bad) + nrow(a_good))\n\tg_prob_bad =  1 - g_prob_good\n\n} else {\n\tprint("source table too sparse, using aggregate distribution")\n\tz = read.table(altSourceF, sep="\\t", header=FALSE)\n\n\tnames(z) <- c("hitChr", "hitStart", "hitEnd", "mystery1", "mystery2", "hitStrand", "eVal", "score", "bias", "hitSpecies", "hitOg", "targetChr", "mystery3", "mystery4", "targetStart", "targetEnd", "mystery5", "mystery6", "targetStrand", "geneLabel", "targetSpecies", "targetOg", "match")\n\n\tz <- cbind(z, score_adj=z$score/(z$hitEnd - z$hitStart))\n\t\n\tz_good <- z[z$match=="match_good",]\n\tz_bad_og <- z[z$match=="match_bad",]\n     z_bad_singleton <- z[z$match=="match_singleton",]\n\tz_bad <- rbind(z_bad_og, z_bad_singleton)\n\n\tprint("fitting good hits")\t\n\tgoodScores=as.data.frame(z_good$score_adj); colnames(goodScores) <- c("t")\n        g_good <- getGamlss(goodScores)\n       \tprint("fitting bad hits")\n\tbadScores=as.data.frame(z_bad$score_adj); colnames(badScores) <- c("t")\n        g_bad <- getGamlss(badScores)\n\n       g_prob_good =  nrow(z_good) / (nrow(z_bad) + nrow(z_good))\n    g_prob_bad =  1 - g_prob_good\n\n}\n\nxg_fun <- function(x) { k=dST1(x, mu=g_good$mu.coefficients ,sigma=exp(g_good$sigma.coefficients), nu=g_good$nu.coefficients, tau=exp(g_good$tau.coefficients)) }\nxb_fun <- function(x) { k=dST1(x, mu=g_bad$mu.coefficients ,sigma=exp(g_bad$sigma.coefficients), nu=g_bad$nu.coefficients, tau=exp(g_bad$tau.coefficients)) }\n\ng_val <- function(x) { (xg_fun(x)*g_prob_good - xb_fun(x)*g_prob_bad) / (g_prob_good*xg_fun(x) + g_prob_bad*xb_fun(x)) }\n\nnone_g_scores = cbind(a_none, g_val=g_val(a_none$score_adj))\n\n#Not actually needed for the algorithm, but useful for other analyses.\ngood_g_scores = cbind(a_good, g_val=g_val(a_good$score_adj))\nbad_g_scores  = cbind(a_bad,  g_val=g_val(a_bad$score_adj))\n\n#Make double sure that very high scores win and very low scores lose\n#(Sometimes they can slip through due to the nature of the calculations)\n\nnone_win = none_g_scores[none_g_scores$g_val >= 0,]\nnone_lose = none_g_scores[none_g_scores$g_val < 0,]\n\nmaxwinscore = max(none_win$score_adj)\nminlosescore = min(none_lose$score_adj)\n\nrescued_win = none_g_scores[(none_g_scores$g_val >= 0 & none_g_scores$score_adj > minlosescore) | none_g_scores$score_adj > maxwinscore,]\nrescued_lose = none_g_scores[(none_g_scores$g_val < 0 & none_g_scores$score_adj < maxwinscore) | none_g_scores$score_adj < minlosescore,]\n\nrescued_win_bad = bad_g_scores[(bad_g_scores$g_val >= 0 & bad_g_scores$score_adj > minlosescore) | bad_g_scores$score_adj > maxwinscore,]\nrescued_lose_bad = bad_g_scores[(bad_g_scores$g_val < 0 & bad_g_scores$score_adj < maxwinscore) | bad_g_scores$score_adj < minlosescore,]\n\nrescued_win_good = good_g_scores[(good_g_scores$g_val >= 0 & good_g_scores$score_adj > minlosescore) | good_g_scores$score_adj > maxwinscore,]\nrescued_lose_good = good_g_scores[(good_g_scores$g_val < 0 & good_g_scores$score_adj < maxwinscore) | good_g_scores$score_adj < minlosescore,]\n\nprint(paste("writing to ", outF, sep=""))\n\nwrite.table(rescued_win, outF, quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\\t")\nwrite.table(rescued_lose, paste(outF, ".rejected", sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\\t")\n\nwrite.table(rescued_win_good, paste(outF, ".good", sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\\t")\nwrite.table(rescued_lose_good, paste(outF, ".good.rejected", sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\\t")\n\nwrite.table(rescued_win_bad, paste(outF, ".bad", sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\\t")\nwrite.table(rescued_lose_bad, paste(outF, ".bad.rejected", sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\\t")\n'
	f=open(path_scriptDestination, "w")
	f.write(str_script)

def unpackFitDistributionScript_noFilter(path_scriptDestination):
	str_script='library("gamlss")\n\nargs <- commandArgs(TRUE)\n\n\nsourceF=args[1]\naltSourceF= args[2]\noutF=args[3]\n\nprint(paste("reading in source table: ", args[1], sep=""))\na <- read.table(sourceF, sep="\\t", header=FALSE)\n\nnames(a) <- c("hitChr", "hitStart", "hitEnd", "mystery1", "mystery2", "hitStrand", "eVal", "score", "bias", "hitSpecies", "hitOg", "targetChr", "mystery3", "mystery4", "targetStart", "targetEnd", "mystery5", "mystery6", "targetStrand", "geneLabel", "targetSpecies", "targetOg", "match")\n\na <- cbind(a, score_adj=a$score/(a$hitEnd - a$hitStart))\n\nh <- hist(a$score_adj, breaks=50, plot=FALSE)\nb <- h$breaks\n\na_none <- a[a$match=="match_none",]\na_good <- a[a$match=="match_good",]\n\n#Declare variables\ng_good = ""\ng_bad = ""\n\ng_prob_good = ""\ng_prob_bad = "" \n\n# Sampling 1000 data points makes the curve-fitting quicker and hardly affects the fit.\n#getGamlss <- function(theData) { print(theData$t); theData_s <- as.data.frame(sample(theData$t, 1000)); colnames(theData_s) <- c("t"); gamlss(t ~ 1, data=theData_s, family="ST1", method=RS(), gd.tol=10000000, c.cyc=0.001, control=gamlss.control(n.cyc=200)) } \ngetGamlss <- function(theData) {}\n\nif(nrow(a_good) > 1000) {\n\tprint("source table is good, going ahead...")\n\ta_bad_og <- a[a$match=="match_bad",]\n\ta_bad_singleton <- a[a$match=="match_singleton",]\n\ta_bad <- rbind(a_bad_og, a_bad_singleton)\n\n\tprint("fitting good hits")\n\tgoodScores=as.data.frame(a_good$score_adj); colnames(goodScores) <- c("t")\n\t\n\tg_good <- getGamlss(goodScores)\n\tprint("fitting bad hits")\n\tbadScores=as.data.frame(a_bad$score_adj); colnames(badScores) <- c("t")\n\tg_bad <- getGamlss(badScores)\n\t\n\tg_prob_good =  nrow(a_good) / (nrow(a_bad) + nrow(a_good))\n\tg_prob_bad =  1 - g_prob_good\n\n} else {\n\tprint("source table too sparse, using aggregate distribution")\n\tz = read.table(altSourceF, sep="\\t", header=FALSE)\n\n\tnames(z) <- c("hitChr", "hitStart", "hitEnd", "mystery1", "mystery2", "hitStrand", "eVal", "score", "bias", "hitSpecies", "hitOg", "targetChr", "mystery3", "mystery4", "targetStart", "targetEnd", "mystery5", "mystery6", "targetStrand", "geneLabel", "targetSpecies", "targetOg", "match")\n\n\tz <- cbind(z, score_adj=z$score/(z$hitEnd - z$hitStart))\n\t\n\tz_good <- z[z$match=="match_good",]\n\tz_bad_og <- z[z$match=="match_bad",]\n     z_bad_singleton <- z[z$match=="match_singleton",]\n\tz_bad <- rbind(z_bad_og, z_bad_singleton)\n\n\tprint("fitting good hits")\t\n\tgoodScores=as.data.frame(z_good$score_adj); colnames(goodScores) <- c("t")\n        g_good <- getGamlss(goodScores)\n       \tprint("fitting bad hits")\n\tbadScores=as.data.frame(z_bad$score_adj); colnames(badScores) <- c("t")\n        g_bad <- getGamlss(badScores)\n\n       g_prob_good =  nrow(z_good) / (nrow(z_bad) + nrow(z_good))\n    g_prob_bad =  1 - g_prob_good\n\n}\n\nnone_good = a_none\n\nprint(paste("writing to ", outF, sep=""))\n\nwrite.table(none_good, outF, quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\\t")\nwrite.table(none_bad, paste(outF, ".rejected", sep=""), quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\\t")'
	f=open(path_scriptDestination, "w")
	f.write(str_script)

def proposeNewGenesIndividual(path_hitsOgIntersectionFileNameAnnotated, path_allHitsOgIntersectionFileNameAnnotated, str_species, path_candidatesFile, hitFilter):
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
	unpackFitDistributionScript(rfile) if hitFilter else unpackFitDistributionScript_noFilter(rfile)
	callFunctionQuiet("rfile=\""+rfile+"\"; Rscript $rfile " + path_hitsOgIntersectionFileNameAnnotated + " " + \
				path_allHitsOgIntersectionFileNameAnnotated + " "  + path_candidatesFile + "; rm $rfile")

def makeNewProteome(path_oldProteome, path_predictedProteinSequences, path_newProteome):
	oldSeqs = [a for a in SeqIO.parse(path_oldProteome, "fasta") if a.seq != ""]
	newSeqs = [a for a in SeqIO.parse(path_predictedProteinSequences, "fasta") if a.seq != ""]
	SeqIO.write(oldSeqs + newSeqs, path_newProteome, "fasta")
				
########################################################
##################### Hint Filter ######################
########################################################

def makeHintsFile(path_goodHits, path_hitsHintsGtf):
	"""Converts our hits output file into a gtf file which can be used to give hints to AUGUSTUS.
	"""
	callFunction("grep -v \"#\" " + path_goodHits + " | sed -r \"s/ +/\\t/g\" | perl -ne 'chomp;@l=split; printf \"%s\\tOrthoFiller\\texonpart\\t%s\\t%s\\t%s\\t%s\\t.\\torthogroup=%s;source=M\\n\", $l[0], $l[1], $l[2], $l[6], $l[5], $l[10]' | sed -r \"s/ +/\\t/g\"  > " + path_hitsHintsGtf)

def hintFscoreFilter(path_augustusParsed, path_hintFile, path_augustusParsedHintFiltered, num_threshold, path_augustusSequences, path_augustusSequencesHintFiltered):
	implementHintFscoreFilter(path_augustusParsed, path_hintFile, path_augustusParsedHintFiltered, num_threshold)
	extractFromFastaByName(path_augustusParsedHintFiltered, path_augustusSequences, path_augustusSequencesHintFiltered)

def goHintFscoreFilter(d_spInfo, hintFilter):
	jobs=[]
	for str_species in d_spInfo:
		if not d_spInfo[str_species]["type"] == "target": continue
		path_augustusParsed = d_spInfo[str_species]["augustusparsed"]
		path_hintFile       = d_spInfo[str_species]["hints"]
		path_augustusParsedHintFiltered = path_augustusParsed+".hintfiltered.gtf"
		d_spInfo[str_species]["augustusparsed_hintfiltered"] = path_augustusParsedHintFiltered
		num_threshold=0.8
		path_augustusSequences = d_spInfo[str_species]["augustussequences"]	
		path_augustusSequencesHintFiltered = path_augustusSequences + ".hintfiltered.fasta"
		d_spInfo[str_species]["augustussequences_hintfiltered"] = path_augustusSequencesHintFiltered	
		if hintFilter:
			jobs.append([hintFscoreFilter, (path_augustusParsed, path_hintFile, path_augustusParsedHintFiltered, num_threshold,  path_augustusSequences, path_augustusSequencesHintFiltered)])
		else:
			d_spInfo[str_species]["augustusparsed_hintfiltered"]    = path_augustusParsed
			d_spInfo[str_species]["augustussequences_hintfiltered"] = path_augustusSequences
	runJobs(jobs, int_cores)

def implementHintFscoreFilter(path_augustusParsed, path_hintFile, path_outFile, num_threshold):
	data = readCsv(path_augustusParsed)
	cds ={}; out =[]
	survivors = 0; 	nonSurvivors = 0
	for i in [a for a in data if a[2] == "CDS"]:
		key=re.sub(r".*transcript_id[ =]\"([^\"]*)\".*", r"\1", i[8])
		if not key in cds: cds[key] = []
		cds[key].append(i)
	data_entries = {}
	gene_entries = {}
	geneLookup = {}
	for i in data:
		if "transcript_id" in i[8]:
			tid=re.sub(r".*transcript_id[ =]\"([^\"]*)\".*", r"\1", i[8])
			gid=re.sub(r".*gene_id[ =]\"([^\"]*)\".*", r"\1", i[8])
			if not tid in data_entries: data_entries[tid] = []
			data_entries[tid].append(i)
			geneLookup[tid]=gid
		elif i[2] == "gene":
			gid = i[8]
			if not gid in gene_entries: gene_entries[gid]=[]
			gene_entries[gid].append(i)
		elif i[2] == "transcript":
			tid=i[8]
			if not tid in data_entries: data_entries[tid]=[]
			data_entries[tid].append(i)
	for tid in cds.keys():
		#print("checking hint scores for " + tid + " from " + path_augustusParsed)
		path_entry	= tempfile.mktemp()
		entry		= cds[tid]
		writeCsv(entry, path_entry)
		gL=sum(int(x[4]) + 1 for x in entry) - sum(int(x[3]) for x in entry)
		path_compatibleHints = tempfile.mktemp()
		path_hintIs = tempfile.mktemp()
		callFunction("bedtools intersect -wa -s -b " + path_entry + " -a " + path_hintFile +" | sort -u > " + path_compatibleHints)
		compatibleHints = readCsv(path_compatibleHints)
		deleteIfPresent(path_compatibleHints)
		hintGroups = {}
		for i in compatibleHints:
			key = i[8]
			if not key in hintGroups: hintGroups[key]=[]
			hintGroups[key].append(i)
		success=False
		for hg in hintGroups:
			path_hg = tempfile.mktemp()
			path_hintIs = tempfile.mktemp()
                        writeCsv(hintGroups[hg], path_hg)
                        callFunction("bedtools intersect -s -a " + path_hg + " -b " + path_entry + " > " + path_hintIs)
			intersection = readCsv(path_hintIs)
			deleteIfPresent(path_hg)
			deleteIfPresent(path_hintIs)
			hL=sum(int(x[4]) + 1 for x in hintGroups[hg]) - sum(int(x[3]) for x in hintGroups[hg])
			iL=sum(int(x[4]) + 1 for x in intersection) - sum(int(x[3]) for x in intersection)
			hR=iL/float(hL)
			hP=iL/float(gL)
			hFsc= 2 * hR * hP / (hR + hP)
			if hFsc >= num_threshold:
				success=True
				break
		if success:
			survivors +=1
			out += data_entries[tid]
			if tid in geneLookup and geneLookup[tid] in gene_entries:
				out += gene_entries[geneLookup[tid]]
		else:
			nonSurvivors += 1
	print("Finished implementing hint score filter for " + path_augustusParsed + ". " + str(survivors) + " survivors and " + str(nonSurvivors) + " non-survivors.")
	writeCsv(out, path_outFile)
	deleteIfPresent(path_entry)

def compareOutputSequences(seq1, seq2):
	return seq1.replace(" ", "").replace("\"", "") == seq2.replace(" ", "").replace("\"", "")

####################################
########### Entry code #############
####################################

def run(d_spInfo, d_sequenceInfoById, orthogroups, singletons, path_resultsDir, path_wDir, path_orthoFinderOutputFile, path_singletonsFile, fullout, int_cores=16, augOnly=False, hitFilter=True, hintFilter=True, splitByChr=False):
	"""Takes orthofinder output and a collection of genome info locations as input.
	   Processes orthogroups in parallel, filters hits, and generates gene models.
	"""
	#####################################################
	# Set off the Augustus training
	#####################################################
	trainingPool = multiprocessing.Pool(int_cores, init_worker)
	print("Training AUGUSTUS")
	trainAugustus(d_spInfo, path_wDir, trainingPool)
	######################################################
	# Prepare some working directories
	######################################################
	path_ogDir      = makeIfAbsent(path_wDir + "/orthogroups")
	path_ogGtfDir   = makeIfAbsent(path_ogDir + "/gtf")
	path_ogAlDir    = makeIfAbsent(path_ogDir + "/alignments")
	path_ogHmmDir   = makeIfAbsent(path_ogDir + "/hmm")
	path_ogHitsDir  = makeIfAbsent(path_ogDir + "/hits")
	path_hmmDbDir   = makeIfAbsent(path_wDir + "/hmmdb")
	path_candidates = makeIfAbsent(path_wDir + "/candidates")
	#####################################################
	# If we're on a second pass, where the first pass was
	# used to create training files, we don't need to 
	# calculate all the hmms and hits.
	#####################################################
	if not augOnly:
		######################################################
		# Set up an hmm database for each species
		######################################################
		stage("2.1. Preparing HMM databases")
		prepareHmmDbInfo(d_spInfo, path_hmmDbDir, splitByChr)
		prepareHmmDbs(d_spInfo, path_hmmDbDir, int_cores)#ql
		#####################################################
		# Produce gtf files for each orthogroup/species pair
		#####################################################
	        stage("2.2. Extracting orthogroup gtf files")
		gtfsForOrthoGroups(path_ogGtfDir, path_orthoFinderOutputFile, path_singletonsFile, d_spInfo, int_cores)#ql
		#####################################################
		# Process each individual orthogroup in parallel
		#####################################################
		proteinSequences = getProteinSequences(d_sequenceInfoById, d_spInfo)
		stage("2.3. Extracting protein fasta sequences")
		print("Getting protein fasta sets for all orthogroups...")
		getProteinFastaFiles(orthogroups, proteinSequences, d_sequenceInfoById, d_spInfo, path_ogAlDir, int_cores)
		stage("2.4. Aligning orthogroup sequences")
		print("Grabbing alignments...")
		getProteinAlignments(orthogroups, path_ogAlDir, int_cores)
		stage("2.5. Extracting nucleotide alignments")
		print("Threading nucleotides through alignments...")
		getNucleotideAlignments(orthogroups, path_ogAlDir, d_sequenceInfoById, d_spInfo, int_cores)
		stage("2.6. Building HMMs")
		print("Grabbing HMMs for each orthogroup...")
		buildHmms(orthogroups, path_ogAlDir, path_ogHmmDir, int_cores)
		stage("2.7. Running HMMs...")
		prepareHitDirs(orthogroups, d_spInfo, path_ogHmmDir, path_ogHitsDir, splitByChr)
		runHmms(orthogroups, d_spInfo, path_ogHmmDir, path_ogHitsDir, int_cores, splitByChr)
		####################################################
		# Start a new pool for processing the hmm outfiles.
		####################################################
		stage("3. Processing HMM output files")
		d_ogIntersectionFileNamesAnnotated = processHmmOutput(d_spInfo, path_candidates, path_ogHitsDir, path_ogGtfDir, splitByChr)
		print("Done processing HMM output files")
		####################################################
		# Concatenate all the files in case we need to 
		# use the aggregate distribution.
		####################################################
		print("Generating concatenated version of HMM output")
		path_allHitsOgIntersectionFileNameAnnotated = path_wDir + "/allSpecies.hitsIntersectionOrthogroups.annotated.bed"
		concatFiles(d_ogIntersectionFileNamesAnnotated.values(), path_allHitsOgIntersectionFileNameAnnotated)
		print("Checking quantity of reference hits...")
		checkReferenceHits(path_allHitsOgIntersectionFileNameAnnotated)
		####################################################
		# Fit a model for each individual species. If data
		# is insufficient, use aggregated data.
		####################################################
		proposeNewGenes(d_spInfo, path_candidates, d_ogIntersectionFileNamesAnnotated, hitFilter, path_allHitsOgIntersectionFileNameAnnotated)
	####################################################
	# Run Augustus. We need the training pool to have 
	# finished by this point. Parse output
	####################################################
#	finally:
	print("Waiting for training to finish before continuing...")
	trainingPool.close()
	trainingPool.join()
#	except KeyboardInterrupt:
#		trainingPool.terminate()
#		trainingPool.join()
#	else:
#		trainingPool.close()
#		trainingPool.join()
	print("Done training.")
        stage("4. Running Augustus")
	path_aug = makeIfAbsent(path_candidates + "/augustus_predictions")
	goAugustus(d_spInfo, path_aug, int_cores)
	####################################################
	# Get a hint f score for the new genes and abandon
	# those genes whose score is not adequate.
	####################################################
	if hintFilter: stage("5.1 Filtering by hint F-score")
	goHintFscoreFilter(d_spInfo, hintFilter)
	####################################################
	# Reinsert the sequences into the proteome and 
	# rerun OrthoFiller
	####################################################
	stage("5.2 Re-running OrthoFinder")
	d_spInfo_modern, path_newProteomesDir = rerunOrthoFinder(path_wDir, d_spInfo)
	####################################################
	# Check genes have ended up in the right orthogroup
	####################################################
	stage("5.3 Running orthogroup membership test")
	jobres = goOrthoGroupTest(path_newProteomesDir, d_spInfo, d_spInfo_modern, orthogroups, int_cores, d_sequenceInfoById)
	stage("5.4 Renaming new genes")
	goRename(d_spInfo, jobres, path_resultsDir, path_newProteomesDir, d_sequenceInfoById, fullout)
	###################################################
	# Leave a friendly goodbye message
	###################################################
	stage("6. Finishing up!")
	finishUp(path_resultsDir, d_spInfo)

def start(path_speciesInfoFile, path_referenceFile, path_orthogroups, path_singletonsFile, path_outDir, path_resultsDir, path_wDir, hitFilter, hintFilter, int_cores, splitByChr, fullout):
	######################################################
	# Read in the locations of the input files and the
	# orthofinder output.
	# MD-CC: will need to have a consistency check for this.:
	######################################################
	stage("1.1. Reading and checking orthofinder output and sequence info files")
	d_spInfo = readInputLocations(path_speciesInfoFile, path_referenceFile)
	d_sequenceInfoById, orthogroups, singletons = readOrthoFinderOutput(path_orthogroups, path_singletonsFile, d_spInfo)
	#####################################################
	# How many genes are there for each species? If any
	# species has less than 100, it needs special training.
	#####################################################
	stage("1.2. Preparing gtf files for Augustus training")
	path_trainingDir = makeIfAbsent(path_wDir + "/training")
	jobs=[]
	for str_species in d_spInfo:
		sequences = [ d_sequenceInfoById[x].seqId for x in d_sequenceInfoById if d_sequenceInfoById[x].species == str_species ]
		if len(sequences) < 100:
			exit("Species " + str_species + " had fewer than 100 genes annotated. For training purposes, all species must have at least 100 genes annotated.")
		path_gtf = d_spInfo[str_species]["gtf"]
		d_spInfo[str_species]["needsTraining"] = (d_spInfo[str_species]["type"] == "target")
		path_gtfForTraining = path_trainingDir + "/" + str_species + ".training.gtf"
#qr		d_spInfo[str_species]["augustusSpecies"]=commands.getstatusoutput("a=`find " + path_wDir+ "/augustus/"+str_species+"/autoAugTrain -name \"tmp_opt*\" -exec stat {} --printf=\"%y\\t%n\\n\" \\;  | sort -t\"-\" -k1,1n -k2,2n -k3,3n | head -n1  | cut -f2`; echo ${a##*/} | sed -r \"s/tmp_opt_//g\"")[1]
		d_spInfo[str_species]["augustusSpecies"] = makeAugustusSpeciesName(str_species)
		d_spInfo[str_species]["gtfForTraining"]  = path_gtfForTraining
		makeGtfTrainingFile(path_gtf, path_gtfForTraining)
		if d_spInfo[str_species]["needsTraining"]: jobs.append([makeGtfTrainingFile, (path_gtf, path_gtfForTraining)])#ql
	runJobs(jobs, int_cores)
	######################################################
	# Run it
	######################################################
	run(d_spInfo, d_sequenceInfoById, orthogroups, singletons, path_resultsDir, path_wDir, path_orthogroups, path_singletonsFile, fullout, int_cores, False, hitFilter, hintFilter, splitByChr)

if __name__ == '__main__':
	# Read in command-line arguments
	parser = argparse.ArgumentParser(description="Run OrthoFiller")
	parser.add_argument("--noHintFilter", action="store_true", dest="noHintFilter", default=False)
	parser.add_argument("--noHitFilter", action="store_true", dest="noHitFilter", default=False)
	parser.add_argument("-c", "--cores", metavar="cores", help="The maximum number of cores you wish to use", dest="CO", default=1)
	parser.add_argument("-o", "--outdir", metavar="outdir", help="The output directory", dest="OD", default="")
	parser.add_argument("-i", "--infoFiles", metavar="info", dest="IN")
	parser.add_argument("-r", "--referenceFiles", metavar="info", dest="RE", default="")
	parser.add_argument("--prep", help="Input data in pre-prepared form", dest="prep", action="store_true")
	parser.add_argument("-g", "--orthogroups", metavar="orthogroups", help="An orthofinder output file (orthogroups)", dest="OG")
	parser.add_argument("-s", "--singletons", metavar="singletons", help="An orthofinder output file (singles)", dest="SN")
	parser.add_argument("-t", "--translationtable", metavar="transtable", help="Which translation table to use", dest="TT")
	parser.add_argument("--checkPrograms", action="store_true", dest="checkOnly", default=False)
	parser.add_argument("--split", action="store_true", dest="SC", default=False)
	parser.add_argument("--fulloutput", action="store_true", dest="FO", default=False)
	args = parser.parse_args()
	prep = args.prep
	
	# Check all the required programsare installed.
	stage("0.1. Checking installed programs")
	checkShell()
	if args.checkOnly: exit("Finished checking installed programs. Everything looks fine.")

	#Check existence and non-confliction of arguments
	if args.IN == None:
		exit("Input file list -i required.")
	if prep:
		if (args.IN == None) | (args.OG == None) | (args.SN == None):
			exit("Option --prep requires options -i [input file list], -g [orthogroups file], and -s [singletons file].")
	else:
		if (args.OG != None) | (args.SN != None):
			exit("Options -g and -s can only be used with option --prep for pre-prepared data.")

	path_outDir = args.OD
	int_cores   = int(args.CO)
	splitByChr  = args.SC

	stage("0.2. Checking and unpacking input data")
	path_resultsDir, path_wDir = prepareOutputFolder(path_outDir)
	#If the data isn't pre-prepared, we must prepare it.
	#Else simply check each file exists. Later we will make sure every entry in the info file exists.
	if not prep:
		path_speciesInfoFile, path_orthoFinderOutputFile, path_singletonsFile, path_referenceFile = prepareFromScratch(args.IN, path_outDir, int_cores, args.RE)
	else:
		path_orthoFinderOutputFile = checkFileExists(args.OG)
		path_singletonsFile        = checkFileExists(args.SN)
		path_speciesInfoFile       = checkFileExists(args.IN)
		if not args.RE == "": path_referenceFile = checkFileExists(args.RE)
	start(path_speciesInfoFile, path_referenceFile, path_orthoFinderOutputFile, path_singletonsFile, path_outDir, path_resultsDir, path_wDir, not args.noHitFilter, not args.noHintFilter, int_cores, splitByChr, args.FO)

