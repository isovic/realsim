#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
import sys
import uuid
import shutil
import subprocess
import random

sys.path.append(SCRIPT_PATH + '')
sys.path.append(SCRIPT_PATH + '/wrappers')

import setup_realsim
import basicdefines
import wrapper_lastal
import wrapper_bwamem
import analyze_venn_of_mappings
import fastqparser
import utility_sam

# For adding quality values to LAST SAM file
from samfilter import add_quality_values

# working with CIGAR profiles
from cigarprofile import CIGARLine, CIGARProfile, loadCProfile, storeCProfile

# global variables for easier manipulation
output_path = basicdefines.INTERMEDIATE_PATH_ROOT_ABS
profiles_path = basicdefines.PROFILES_PATH_ABS


# Returns absolute path for a given profile, prints out error and exits
# if a given profile doesn't exist
def get_cprofile_path(cprofile):
	profile_path = os.path.join(basicdefines.PROFILES_PATH_ABS, cprofile + '.cpf')
	if not os.path.exists(profile_path):
		sys.stderr.write('\n\nCIGAR profile %s doesn\'t exist! Exiting ...' % cprofile)
		verbose_usage_and_exit()

	return profile_path


# Extracts error profile as CIGAR profile, from a SAM file
# Considers only the best alignment for each SAM line and
# considers only SAM lines whose chosen quality is greater then 0
def extract_cprofile_from_LAST_SAM(reference_path, sam_path_lastal, profile='cprofile'):

	sys.stderr.write('\n')

	qnames_with_multiple_alignments = {};
	sys.stderr.write('Loading SAM file into hash...\n');
	[sam_hash, sam_hash_num_lines, sam_hash_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_path_lastal, qnames_with_multiple_alignments);

	sys.stderr.write('Loading reference from fastq file...\n')
	reference_fastq = fastqparser.read_fastq(reference_path)
	reference_seq = reference_fastq[1][0]
	GCcontent = fastqparser.gc_content(reference_seq)

	# Initialize mutation count table
	utility_sam.init_mutCntTable()

	cprofile = CIGARProfile(profile, GCcontent)
	for qname in sam_hash.iterkeys():
		samline = sam_hash[qname]
		# TODO: Revise CalcExtendedCIGAR to receive only relevant part of reference
		cigar, gccontent = samline[0].CalcExtendedCIGARandGCContent(reference_seq)
		pos = samline[0].pos
		quality = samline[0].chosen_quality
		quals = samline[0].qual
		if quality > 0:
			cprofile.appendCLine(CIGARLine(cigar, pos, quality, gccontent, quals))

	cprofile.mutCntTable = utility_sam.mutCntTable

	cprofilefilepath = os.path.join(profiles_path, profile + '.cpf')
	sys.stderr.write('Saving CIGAR profile to:%s\n' % cprofilefilepath)
	storeCProfile(cprofilefilepath, cprofile)


def extract_cprofile_with_LAST(reference_path, reads_path, technology, profile='cprofile'):
	if (not os.path.exists(basicdefines.TOOLS_ROOT_ABS) or not os.path.exists(basicdefines.ALIGNERS_PATH_ROOT_ABS)):
		sys.stderr.write('ERROR: Required tools not set up. Rerun the script with the "setup" option.\nExiting.\n\n')
		exit(1)

	uuid_string = str(uuid.uuid4())			# Generate a random UUID so that multiple runs don't clash.

	# copy reference to output path with the same uuid
	ref_filename = "REF_" + uuid_string + ".fastq"
	new_ref_path = os.path.join(output_path, ref_filename)
	shutil.copy(reference_path, new_ref_path)

	if (not os.path.exists(output_path)):
		sys.stderr.write('Creating folder "%s".\n' % (output_path))
		os.makedirs(output_path)

	sam_path_lastal = wrapper_lastal.run(reads_path, reference_path, technology, output_path, uuid_string)

	basename = os.path.basename(sam_path_lastal)
	dirname = os.path.dirname(sam_path_lastal)

	# Adding Q in front of filename for sam file with quals
	sam_with_quals_path = os.path.join(dirname, 'Q' + basename)

	# Adding quality values to LAST SAM file
	add_quality_values(sam_path_lastal, reads_path, sam_with_quals_path)

	extract_cprofile_from_LAST_SAM(reference_path, sam_with_quals_path, profile)


def extract_cprofile_from_SAMs(reference_path, sam_path_bwamem, sam_path_lastal, profile='cprofile'):

	distance_threshold = -1			# This currently does nothing.
	# compare_two_sams returns three dicts:
	#	1. sam_hash1	-	hashed SAM file. Key of the dict is the qname of a SAM entry, and value is a list of SAMLine objects. List contains all alignments for the same qname
	#						and is composed of at least one element. The list is sorted in the descending order of the MAPQ parameter is available (i.e. if it is different from 255),
	#						otherwise it is sorted according to the AS value (again, if available).
	#	2. sam_hash2	-	same as the previous parameter, but for the second given SAM file.
	#	3. distance_to_qname_hash	-	Key of the dict is the distance between the mapping positions of the same read accross two SAM files. Value for each key contains
	#									a list of qnames at that particular distance. To retreieve concrete info on a specific alignment for the given qname, one can
	#									query the sam_hash1 and sam_hash2 objects.
	[sam_hash_bwamem, sam_hash_lastal, distance_to_qname_hash] = analyze_venn_of_mappings.compare_two_sams(sam_path_bwamem, sam_path_lastal, distance_threshold, out_summary_prefix='')

	extract_cprofile_from_dict(reference_path, sam_hash_bwamem, sam_hash_lastal, distance_to_qname_hash, profile)


# Extracts error profile as CIGAR profile, from dictionaries
def extract_cprofile_from_dict(reference_path, sam_hash_bwamem, sam_hash_lastal, distance_to_qname_hash, profile='cprofile'):
	# Load reference
	sys.stderr.write('\n')
	sys.stderr.write('Loading reference from fastq file...\n')
	reference_fastq = fastqparser.read_fastq(reference_path)
	reference_seq = reference_fastq[1][0]
	GCcontent = fastqparser.gc_content(reference_seq)

	# Initialize mutation count table
	utility_sam.init_mutCntTable()

	# Using only those qnames that BWA and LAST aligned to the same position
	# Using both CIGAR strings
	qnames = distance_to_qname_hash[0]
	cprofile = CIGARProfile(profile, GCcontent)
	for qname in qnames:
		lastsamline = sam_hash_lastal[qname]
		# TODO: Revise CalcExtendedCIGAR to receive only relevant part of reference
		cigar, gccontent = lastsamline[0].CalcExtendedCIGARandGCContent(reference_seq)
		pos = lastsamline[0].pos
		quality = lastsamline[0].chosen_quality
		quals = lastsamline[0].qual
		cprofile.appendCLine(CIGARLine(cigar, pos, quality, gccontent, quals))

		bwasamline = sam_hash_bwamem[qname]
		cigar, gccontent = bwasamline[0].CalcExtendedCIGARandGCContent(reference_seq)
		pos = bwasamline[0].pos
		quality = bwasamline[0].chosen_quality
		qual = bwasamline[0].qual
		cprofile.appendCLine(CIGARLine(cigar, pos, quality, gccontent, qual))

	cprofile.mutCntTable = utility_sam.mutCntTable

	cprofilefilepath = os.path.join(profiles_path, profile + '.cpf')
	sys.stderr.write('Saving CIGAR profile to:%s\n' % cprofilefilepath)
	storeCProfile(cprofilefilepath, cprofile)


# Not used. MSA was discarded as an idea
def extract_MSAprofile_from_dict(reference_path, sam_hash_bwamem, sam_hash_lastal, distance_to_qname_hash):
	# Load reference
	sys.stderr.write('\n')
	sys.stderr.write('Loading reference from fastq file...\n')
	reference_fastq = fastqparser.read_fastq(reference_path)
	reference_seq = reference_fastq[1][0]
	# print reference_fastq

	sys.stderr.write('Calculating maximum read size from SAM files ...\n')
	# Calculating largest read size, it seems enough to go through only one hash
	maxreadsize = 0
	minreadsize = 10000000		# at the moment this seems apropriate
	for (qname, samlines) in sam_hash_bwamem.iteritems():
		# looking only at the first samline
		samline = samlines[0]
		# print samline.seq
		if maxreadsize < len(samline.seq):
			maxreadsize = len(samline.seq)

		if minreadsize > len(samline.seq):
			minreadsize = len(samline.seq)

	windowsize = 2*maxreadsize

	print "Maxreadsize: " + str(maxreadsize)
	print "Minreadsize: " + str(minreadsize)

	# For lambda phage reads:
	# Maxreadsize: 21158
	# Minreadsize: 1296

	# For testing purposes taking a window of size 10k bp from the start of the reference genome
	# And taking all reads that both aligners completely placed inside isst
	# Trying to run MAFFT on it, just to see what i get
	# 458 lines, including reference
	# windowsize = 10000

	# for pbsim reads window size is 1000
	windowsize = 10000
	windowstart = 100000
	windowend = windowstart + windowsize-1
	seqcounter = 0

	# At the moment padding doesn't seem to make sense
	# padding = True
	padding = False

	# Constructing an input fasta file
	msainputpath = os.path.join(output_path, 'msainputfile.fasta')
	with open(msainputpath, 'w') as msainputfile:

		# Writing reference as the first sequence in fasta
		sys.stderr.write('Writing MSA input fasta file...\n')
		msainputfile.write('>Reference_' + reference_fastq[0][0] + '\n')
		msainputfile.write(reference_seq[windowstart:windowend] + '\n')

		sys.stderr.write('Writing reads...')
		# getting sam_lines that are within current window
		for (qname, samlines) in sam_hash_bwamem.iteritems():
			sys.stderr.write('\rWritting reads...')
			# looking only at the first samline
			samline = samlines[0]

			# if the line is within current window
			if (samline.clipped_pos >= windowstart) and (samline.clipped_pos + len(samline.seq) <= windowend):
				#check if last also places it within the current window
				if qname in sam_hash_lastal:
					samline2 = sam_hash_lastal[qname][0]
					if (samline2.clipped_pos >= windowstart) and (samline2.clipped_pos + len(samline2.seq) <= windowend):
						# write sequence in the file
						sys.stderr.write('\rWritting ' + qname)
						msainputfile.write('>' + qname + '\n')
						# extending each read to window length by adding Ns at the begining and at the end
						if padding:
							startpad = samline.clipped_pos - windowstart				# number of Ns in front of the sequence
							endpad = windowend - samline.clipped_pos - len(samline.seq) # number of Ns after the sequence
							seq = 'N'*startpad + samline.seq + 'N'*endpad
						else:
							seq = samline.seq
						msainputfile.write(seq + '\n')
						seqcounter += 1

	sys.stderr.write('\n\nMSA input contains %d sequences\n' % seqcounter)
	if seqcounter == 0:
		sys.stderr.write('Window contains no reads!\nExiting!\n')
		exit(0)

	# TODO:
	#	sort both sam_line lists according to clipped_pos
	#	move a window over reference genome and sam_line lists
	#	create an input fasta file for MAFFT
	#	run MAFFT
	#	format MAFFT output
	#	do something with the results

	# running MAFFT
	mafftoutputpath = os.path.join(output_path, 'mafftoutputfile.fasta')
	sys.stderr.write('\n\n')
	sys.stderr.write('Running MAFFT...\n')
	command = 'mafft --localpair --maxiterate 1000 %s > %s' % (msainputpath, mafftoutputpath)
	subprocess.call(command, shell=True)

	# combining MAFFT output with SAM data
	sys.stderr.write('\nCombining MAFFT output with SAM files...\n')
	mafftcombinedpath = os.path.join(output_path, 'mafftcombined.fasta')

	analyze_venn_of_mappings.combine_msa_output(mafftcombinedpath, mafftoutputpath, sam_hash_bwamem, sam_hash_lastal, windowstart)

	# # running MUSCLE
	# muscleoutputpath = os.path.join(output_path, 'muscleoutputfile.fasta')
	# sys.stderr.write('\n\n')
	# sys.stderr.write('Running MUSCLE...\n')
	# command = 'muscle -in %s -out %s' % (msainputpath, muscleoutputpath)
	# subprocess.call(command, shell=True)

	# # combining MUSCLE output with SAM data
	# sys.stderr.write('\nCombining MUSCLE output with SAM files...\n')
	# musclecombinedpath = os.path.join(output_path, 'musclecombined.fasta')

	# analyze_venn_of_mappings.combine_msa_output(musclecombinedpath, muscleoutputpath, sam_hash_bwamem, sam_hash_lastal, windowstart)

	# running CLUSTALO
	clustalooutputpath = os.path.join(output_path, 'clustalooutputfile.fasta')
	sys.stderr.write('\n\n')
	sys.stderr.write('Running clustalo...\n')
	command = 'clustalo -i %s -o %s -v' % (msainputpath, clustalooutputpath)
	subprocess.call(command, shell=True)


	# combining clustalo output with SAM data
	sys.stderr.write('\n\n')
	sys.stderr.write('Combining CLUSTALO output with SAM files...\n')
	clustalocombinedpath = os.path.join(output_path, 'clustalocombined.fasta')

	analyze_venn_of_mappings.combine_msa_output(clustalocombinedpath, clustalooutputpath, sam_hash_bwamem, sam_hash_lastal, windowstart)

	# Example of how to access reads at various distances between two SAM files:
	# sys.stdout.write('Histogram of alignment distances:\n');
	# i = 0;
	# for distance_key in sorted(distance_to_qname_hash.keys()):
	# 	qnames = distance_to_qname_hash[distance_key];
	# 	sys.stdout.write('[%d] Distance: %d, number of reads: %d\n' % (i, distance_key, len(qnames)));
	# 	i += 1;
	# 	if (distance_key > 10):
	# 		break;

	# Further example, working but commented out not to clog the stdout.
	# try:
	# 	qnames = distance_to_qname_hash[0];
	# 	sys.stdout.write('SAM line from BWA-MEM mapping:\n');
	# 	sys.stdout.write('%s\n\n' % sam_hash_bwamem[qnames[0]][0].original_line);

	# 	sys.stdout.write('SAM line from LAST mapping:\n');
	# 	sys.stdout.write('%s\n\n' % sam_hash_lastal[qnames[0]][0].original_line);
	# except Exception, e:
	# 	sys.stderr.write('%s\n' % e);



def extract_cprofile_from_path(reference_path, reads_path, technology, profile = 'cprofile'):
	if (not os.path.exists(basicdefines.TOOLS_ROOT_ABS) or not os.path.exists(basicdefines.ALIGNERS_PATH_ROOT_ABS)):
		sys.stderr.write('ERROR: Required tools not set up. Rerun the script with the "setup" option.\nExiting.\n\n')
		exit(1)

	# output_path = '%s/../intermediate' % SCRIPT_PATH;
	# output path is a global variable
	uuid_string = str(uuid.uuid4())			# Generate a random UUID so that multiple runs don't clash.

	# copy reference to output path with the same uuid
	ref_filename = "REF_" + uuid_string + ".fastq"
	new_ref_path = os.path.join(output_path, ref_filename)
	shutil.copy(reference_path, new_ref_path)

	if (not os.path.exists(output_path)):
		sys.stderr.write('Creating folder "%s".\n' % (output_path))
		os.makedirs(output_path)

	sam_path_bwamem = wrapper_bwamem.run(reads_path, reference_path, technology, output_path, uuid_string)
	sam_path_lastal = wrapper_lastal.run(reads_path, reference_path, technology, output_path, uuid_string)

	# Adding quality values to LAST SAM file
	basename = os.path.basename(sam_path_lastal)
	dirname = os.path.dirname(sam_path_lastal)

	# Adding Q in front of filename for sam file with quals
	sam_with_quals_path = os.path.join(dirname, 'Q' + basename)

	add_quality_values(sam_path_lastal, reads_path, sam_with_quals_path)

	extract_cprofile_from_SAMs(new_ref_path, sam_path_bwamem, sam_with_quals_path, profile)


def calculate_statistics(reference_path, sam_path_bwamem, sam_path_lastal):

	distance_threshold = -1;			# This currently does nothing.
	# compare_two_sams returns three dicts:
	#	1. sam_hash1	-	hashed SAM file. Key of the dict is the qname of a SAM entry, and value is a list of SAMLine objects. List contains all alignments for the same qname
	#						and is composed of at least one element. The list is sorted in the descending order of the MAPQ parameter is available (i.e. if it is different from 255),
	#						otherwise it is sorted according to the AS value (again, if available).
	#	2. sam_hash2	-	same as the previous parameter, but for the second given SAM file.
	#	3. distance_to_qname_hash	-	Key of the dict is the distance between the mapping positions of the same read accross two SAM files. Value for each key contains
	#									a list of qnames at that particular distance. To retreieve concrete info on a specific alignment for the given qname, one can
	#									query the sam_hash1 and sam_hash2 objects.
	[sam_hash_bwamem, sam_hash_lastal, distance_to_qname_hash] = analyze_venn_of_mappings.compare_two_sams(sam_path_bwamem, sam_path_lastal, distance_threshold, out_summary_prefix='')

	# output_path = '%s/../intermediate' % SCRIPT_PATH;
	# output path is a global variable

	# looking only at reads (queries) that both aligners placed at exactly the same position
	# looking only at the best alignment for each query
	cntdifflen = 0
	cntdiff = 0
	cnteq = 0
	for qname in distance_to_qname_hash[0]:
		bwa_samline = sam_hash_bwamem[qname][0]
		# extracting bwa cigar without clipping
		cigar = bwa_samline.cigar
		begin = 0
		while cigar[begin].isdigit():
			begin += 1
		if cigar[begin] != 'S':
			begin = 0
		else:
			begin += 1

		end = -1
		if (cigar[end] == 'S'):
			end = -2
			while cigar[end].isdigit():
				end -= 1
			end += 1

		else:
			end = 0

		if end == 0:
			bwa_cigar = cigar[begin:]
		else:
			bwa_cigar = cigar[begin:end]

		last_samline = sam_hash_lastal[qname][0]
		cigar = last_samline.cigar
		begin = 0
		while cigar[begin].isdigit():
			begin += 1
		if cigar[begin] != 'H':
			begin = 0
		else:
			begin += 1

		end = -1
		if (cigar[end] == 'H'):
			end = -2
			while cigar[end].isdigit():
				end -= 1
			end += 1

		else:
			end = 0

		if end == 0:
			last_cigar = cigar[begin:]
		else:
			last_cigar = cigar[begin:end]

		if len(bwa_samline.cigar) != len(last_samline.cigar):
			cntdifflen += 1
		if bwa_cigar != last_cigar:
			cntdiff += 1
		else:
			cnteq += 1

	sys.stdout.write('\n%d SAM lines placed at exactly the same position!' % len(distance_to_qname_hash[0]))
	sys.stdout.write('\n%d SAM lines have CIGARs of different length!' % cntdifflen)
	sys.stdout.write('\n%d SAM lines have different CIGARs!' % cntdiff)
	sys.stdout.write('\n%d SAM lines have the same CIGARs !' % cnteq)
	sys.stdout.write('\n')


# Cleans up intermediate files by deleting all files in the corresponding folder
def cleanup():
	sys.stdout.write('\nRemoving intermediate files%s!')
	for filename in os.listdir(basicdefines.INTERMEDIATE_PATH_ROOT_ABS):
		path = os.path.join(basicdefines.INTERMEDIATE_PATH_ROOT_ABS, filename)
		sys.stdout.write('\nRemoving: %s' % path)
		os.remove(path)


# Lists all available CIGAR profiles
# All profiles are stored in a separate folder in files ending with .cpf
def list_profiles():
	sys.stdout.write('\nAvailable profiles:')
	for filename in os.listdir(basicdefines.PROFILES_PATH_ABS):
		if filename.endswith('.cpf'):
			profilename = filename[:-4]
			sys.stdout.write('\n\t- %s' % profilename)

	sys.stdout.write('\n')



def simulate_with_cprofile(cprofile, reference_path, output_file, numreads = None, coverage = None, header = None, fastq = True):

	profile_path = get_cprofile_path(cprofile)

	# Load reference
	sys.stderr.write('\n')
	sys.stderr.write('Loading reference from fastq file...\n')
	reference_fastq = fastqparser.read_fastq(reference_path)
	reference_seq = reference_fastq[1][0]
	ref_header = reference_fastq[0][0]

	# Header is written into every header line in simulated fastq file
	# If not set, reference header is used
	# Replacing spaces with underscores
	if header == None:
		header = ref_header
	header = header.replace(' ', '_')

	# Load CProfile
	sys.stderr.write('Loading CIGAR profile...\n')
	cprofile = loadCProfile(profile_path)

	readspath = os.path.join(output_path, output_file)

	# Callculating number of reads from coverage if necessary
	if numreads is None:
		sumlength = 0
		clines = cprofile.clines
		for cline in clines:
			sumlength += cline.length()
		avglength = sumlength / len(clines)

		# numreads = len(reference)*coverage / avg_read_length
		numreads = len(reference_seq) * coverage / avglength

	if fastq:
		firstchar = '@'
	else:
		firstchar = '>'

	# Writing simulated reads
	sys.stderr.write('Simulating reads...\n')
	with open(readspath, 'w') as readsfile:
		for i in xrange(numreads):
			qname = '%s_simulated_read_%d' % (header, i)
			read, quals = cprofile.generateRandomRead(reference_seq)
			readsfile.write(firstchar + qname + '\n')
			readsfile.write(read + '\n')
			# If quals exist in the profile and files needs to be fastq, write them to the file as well
			if quals is not None and fastq:
				readsfile.write('+' + qname + '\n')
				readsfile.write(quals + '\n')


def spike_with_cprofile(cprofile, reference_path, reads_file,  output_file, numreads = None, coverage = None, header = None, fastq = True):

	profile_path = get_cprofile_path(cprofile)

	# Generate temp output path for simulating new reads
	uuid_string = str(uuid.uuid4())			# Generate a random UUID so that multiple runs don't clash.
	temp_filename = 'TEMP_' + uuid_string + os.path.splitext(output_file)[1]
	temp_output_file = os.path.join(basicdefines.INTERMEDIATE_PATH_ROOT_ABS, temp_filename)

	simulate_with_cprofile(profile_path, reference_path, temp_output_file, numreads, coverage, header, fastq)
	combine_files(temp_output_file, reads_file, output_file)



def calculate_statistics_for_CProfile(cprofile):

	profile_path = get_cprofile_path(cprofile)

	# Load CProfile
	sys.stderr.write('Loading CIGAR profile...\n')
	cprofile = loadCProfile(profile_path)

	numlines = len(cprofile.clines)
	sumlen = 0.0
	minlen = 1000000
	maxlen = 0
	numgoodquals = 0
	numbadquals = 0
	gcsum = 0.0
	gcmin = 1
	gcmax = 0
	for cline in cprofile.clines:
		length = cline.length()
		sumlen += length
		if minlen > length:
			minlen = length
		if maxlen < length:
			maxlen = length

		# LAST SAM file has '*' as quals, checking for length 5 just in case
		if len(cline.quals) < 5:
			numbadquals += 1
		else:
			numgoodquals += 1

		gcsum += cline.GCcontent
		if gcmin > cline.GCcontent:
			gcmin = cline.GCcontent
		if gcmax < cline.GCcontent:
			gcmax = cline.GCcontent

	sys.stdout.write('\n')
	sys.stdout.write('CIGAR profile: %s\n' % cprofile.name)
	sys.stdout.write('Profile GC content: %f\n' % cprofile.GCcontent)
	sys.stdout.write('Number of lines: %d\n' % len(cprofile.clines))
	sys.stdout.write('Profile read length (MIN|AVG|MAX): (%d | %d | %d)\n' % (minlen, sumlen/len(cprofile.clines), maxlen))
	sys.stdout.write('GC content per read (MIN|AVG|MAX): (%0.4f | %0.4f | %0.4f)\n' % (gcmin, gcsum/len(cprofile.clines), gcmax))
	sys.stdout.write('Good / bad quals: (%d / %d)\n' % (numgoodquals, numbadquals))


# Combines two fasta/fastq files by randomly selecting a file in each step.
# Used for spiking existing reads with simulated reads
# TODO: ATM loading everything first into memory. Should consider reading files line by line
#       to significantly reduce memory consumption
def combine_files(path1, path2, resultspath):
	# Loading first file
	sys.stderr.write('\nLoading first file: %s' % path1)
	[headers1, seqs1, quals1] = fastqparser.read_fastq(path1)

	# Loading second file
	sys.stderr.write('\nLoading second file: %s' % path2)
	[headers2, seqs2, quals2] = fastqparser.read_fastq(path2)

	if resultspath.endswith('.fasta') or resultspath.endswith('.fa'):
		filetype = 'fasta'
	elif resultspath.endswith('.fastq') or resultspath.endswith('.fq'):
		filetype = 'fastq'
	else:
		sys.stderr.write('\n\nInvalid results file')
		return

	token = ''
	if filetype == 'fasta':
		token = '>'
	if filetype == 'fastq':
		token = '@'

	numreads1 = len(headers1)
	numreads2 = len(headers2)
	totalreads = numreads1 + numreads2

	index1 = index2 = 0

	sys.stderr.write('\nMixing reads .... \n' )
	with open(resultspath, 'w') as resultsfile:
		while index1 < numreads1 and index2 < numreads2:
			rnd = random.randint(1, totalreads)
			if rnd <= numreads1:
				resultsfile.write(token + headers1[index1] + '\n')
				resultsfile.write(seqs1[index1] + '\n')
				if filetype == 'fastq':
					resultsfile.write('+' + headers1[index1] + '\n')
					resultsfile.write(quals1[index1] + '\n')
				index1 += 1
			else:
				resultsfile.write(token + headers2[index2] + '\n')
				resultsfile.write(seqs2[index2] + '\n')
				if filetype == 'fastq':
					resultsfile.write('+' + headers2[index2] + '\n')
					resultsfile.write(quals2[index2] + '\n')
				index2 += 1

		while index1 < numreads1:
			resultsfile.write(token + headers1[index1] + '\n')
			resultsfile.write(seqs1[index1] + '\n')
			if filetype == 'fastq':
				resultsfile.write('+' + headers1[index1] + '\n')
				resultsfile.write(quals1[index1] + '\n')
			index1 += 1

		while index2 < numreads2:
			resultsfile.write(token + headers2[index2] + '\n')
			resultsfile.write(seqs2[index2] + '\n')
			if filetype == 'fastq':
				resultsfile.write('+' + headers2[index2] + '\n')
				resultsfile.write(quals2[index2] + '\n')
			index2 += 1




def verbose_usage_and_exit():
	sys.stderr.write('RealSim - sequence simulator based on empirical observations.\n')
	sys.stderr.write('\n')
	sys.stderr.write('Usage:\n')
	sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
	sys.stderr.write('\n')
	sys.stderr.write('\tmode:\n')
	sys.stderr.write('\t\tsetup\n')
	sys.stderr.write('\t\tcleanup\n')
	sys.stderr.write('\t\textract\n')
	sys.stderr.write('\t\textract_with_LAST\n')
	sys.stderr.write('\t\textract_from_SAMs\n')
	sys.stderr.write('\t\textract_from_SAM\n')
	sys.stderr.write('\t\tsimulate_with_CProfile\n')
	sys.stderr.write('\t\tspike_with_CProfile\n')
	sys.stderr.write('\t\tcombine_files\n')
	sys.stderr.write('\t\tcombine_profiles\n')
	sys.stderr.write('\t\tlist_profiles\n')
	sys.stderr.write('\t\tstatistics_from_SAMs\n')
	sys.stderr.write('\t\tstatistics_for_CProfile\n')
	sys.stderr.write('\n')
	exit(0)

if __name__ == '__main__':
	if (len(sys.argv) < 2):
		verbose_usage_and_exit()

	mode = sys.argv[1]

	if (mode == 'setup'):
		if (len(sys.argv) != 2):
			sys.stderr.write('Setup the folder structures and install necessary tools.\n')
			sys.stderr.write('Requires no additional parameters to run.\n')
			sys.stderr.write('\n')
			exit(1)

		setup_realsim.setup_all()

	elif (mode == 'cleanup'):
		if (len(sys.argv) != 2):
			sys.stderr.write('Cleans up intermediate files.\n')
			sys.stderr.write('Requires no additional parameters to run.\n')
			sys.stderr.write('\n')
			exit(1)

		cleanup()

	elif (mode == 'extract'):
		if (len(sys.argv) != 6):
			sys.stderr.write('Extracts the error profile from a given set of reads and a reference genome.\n')
			sys.stderr.write('Error profile is stored in a file as CIGARProfile.\n')
			sys.stderr.write('\n')
			sys.stderr.write('Usage:\n')
			sys.stderr.write('\t%s %s <reference> <reads> technology <profile>' % (sys.argv[0], sys.argv[1]))
			sys.stderr.write('\n')
			sys.stderr.write('\ttechnology - "illumina", "pacbio" or "nanopore"\n')
			sys.stderr.write('\n')
			exit(1)

		reference_path = sys.argv[2]
		reads_path = sys.argv[3]
		technology = sys.argv[4]
		if not (technology in ['illumina', 'pacbio', 'nanopore']):
			sys.stderr.write('Unsuported sequencing technology!\n\n')
			verbose_usage_and_exit()
		profile = sys.argv[5]
		extract_cprofile_from_path(reference_path, reads_path, technology, profile)

	elif (mode == 'extract_with_LAST'):
		if (len(sys.argv) != 6):
			sys.stderr.write('Extracts the error profile using only LAST aligner.\n')
			sys.stderr.write('Error profile is stored in a file as CIGARProfile.\n')
			sys.stderr.write('\n')
			sys.stderr.write('Usage:\n')
			sys.stderr.write('\t%s %s <reference> <reads> technology <profile>' % (sys.argv[0], sys.argv[1]))
			sys.stderr.write('\n')
			sys.stderr.write('\ttechnology - "illumina", "pacbio" or "nanopore"\n')
			sys.stderr.write('\n')
			exit(1)

		reference_path = sys.argv[2]
		reads_path = sys.argv[3]
		technology = sys.argv[4]
		if not (technology in ['illumina', 'pacbio', 'nanopore']):
			sys.stderr.write('Unsuported sequencing technology!\n\n')
			verbose_usage_and_exit()
		profile = sys.argv[5]
		extract_cprofile_with_LAST(reference_path, reads_path, technology, profile)

	# So that aligners do not need to be run every time
	# Speeds up testing
	elif (mode == 'extract_from_SAMs'):
		if (len(sys.argv) != 7):
			sys.stderr.write('Extracts the error profile from a given set of SAM files (BWA and LAST).\n')
			sys.stderr.write('Only SAM lines aligned to the same place by both aligners are used.\n')
			sys.stderr.write('Error profile is stored in a file as CIGARProfile.\n')
			sys.stderr.write('\n')
			sys.stderr.write('Usage:\n')
			sys.stderr.write('\t%s %s <reference> <BWA_SAM> <LAST_SAM> technology <profile>' % (sys.argv[0], sys.argv[1]))
			sys.stderr.write('\n')
			sys.stderr.write('\ttechnology - "illumina", "pacbio" or "nanopore"\n')
			sys.stderr.write('\n')
			exit(1)

		reference_path = sys.argv[2]
		sam_path_bwamem = sys.argv[3]
		sam_path_lastal = sys.argv[4]
		technology = sys.argv[5]
		if not (technology in ['illumina', 'pacbio', 'nanopore']):
			sys.stderr.write('Unsuported sequencing technology!\n\n')
			verbose_usage_and_exit()
		profile = sys.argv[6]
		extract_cprofile_from_SAMs(reference_path, sam_path_bwamem, sam_path_lastal, profile)

	elif (mode == 'extract_from_SAM'):
		if (len(sys.argv) != 6):
			sys.stderr.write('Extracts the error profile using only one SAM file (intended for LAST).\n')
			sys.stderr.write('Error profile is stored in a file as CIGARProfile.\n')
			sys.stderr.write('\n')
			sys.stderr.write('Usage:\n')
			sys.stderr.write('\t%s %s <reference> <SAM_file> technology <profile>' % (sys.argv[0], sys.argv[1]))
			sys.stderr.write('\n')
			sys.stderr.write('\ttechnology - "illumina", "pacbio" or "nanopore"\n')
			sys.stderr.write('\n')
			exit(1)

		reference_path = sys.argv[2]
		sam_path_lastal = sys.argv[3]
		technology = sys.argv[4]
		if not (technology in ['illumina', 'pacbio', 'nanopore']):
			sys.stderr.write('Unsuported sequencing technology!\n\n')
			verbose_usage_and_exit()
		profile = sys.argv[5]
		extract_cprofile_from_LAST_SAM(reference_path, sam_path_lastal, profile)

	elif (mode == 'simulate_with_CProfile'):
		if (len(sys.argv) < 5):
			sys.stderr.write('Applies an extracted CIGAR profile to a given reference in order to simulate reads.\n')
			sys.stderr.write('Output is generated as a FASTQ file.\n')
			sys.stderr.write('\n')
			sys.stderr.write('Usage:\n')
			sys.stderr.write('\t%s %s <profile> <reference> <output file> options\n' % (sys.argv[0], sys.argv[1]))
			sys.stderr.write('\n')
			sys.stderr.write('\toptions:\n')
			sys.stderr.write('\t\t-h header\n')
			sys.stderr.write('\t\t-n number_of_reads\n')
			sys.stderr.write('\t\t-c coverage\n')
			sys.stderr.write('\t\t-t type (fasta or fatsq, fastq is default)\n')
			sys.stderr.write('\n')
			exit(1)

		profile_path = sys.argv[2]
		reference_path = sys.argv[3]
		output_file = sys.argv[4]

		# Look at extension of the output file to decide whether to generate reads in fasta or fastq format
		# Format can also be specified by parameter -t. Parameter is stronger then file extension.
		# If nothing is specified, fastq reads are generated
		fastq = True
		if output_file.endswith('.fa') or output_file.endswith('.fasta'):
			fastq = False

		# parsing options
		header = numreads = coverage = None
		for i in range(5, len(sys.argv), 2):
			if sys.argv[i] == '-h':
				header = sys.argv[i+1]
			elif sys.argv[i] == '-n':
				try:
					numreads = int(sys.argv[i+1])
				except:
					sys.stderr.write('Wrong value for parameter "number of reads"! Value needs to be an int.\n\n')
					exit()
			elif sys.argv[i] == '-c':
				try:
					coverage = int(sys.argv[i+1])
				except:
					sys.stderr.write('Wrong value for parameter "coverage"! Value needs to be an int.\n\n')
					exit()
			elif sys.argv[i] == '-t':
				ftype = sys.argv[i+1].lower()
				if ftype == 'fastq':
					fastq = True
				elif ftype == 'fasta':
					fastq = False
				else:
					sys.stderr.write('Parameter -t goes with \'fasta\' or \'fastq\'.\n\n')
					exit()
			else:
				sys.stderr.write('Invalid parameter!.\n\n')
				exit()

		if numreads is None and coverage is None:
			sys.stderr.write('Either "coverage" or "number of reads" parameter must be specified.\n\n')
			exit()

		simulate_with_cprofile(profile_path, reference_path, output_file, numreads, coverage, header, fastq)

	elif (mode == 'spike_with_CProfile'):
		if (len(sys.argv) != 6):
			sys.stderr.write('Applies an extracted CIGAR profile to a given reference in order to simulate reads.\n')
			sys.stderr.write('Simulated reads are inserted into a given reads file (file is spiked with simulated reads).\n')
			sys.stderr.write('\n')
			sys.stderr.write('Usage:\n')
			sys.stderr.write('\t%s %s <profile> <reference> <reads file> <output file> options\n' % (sys.argv[0], sys.argv[1]))
			sys.stderr.write('\n')
			sys.stderr.write('\toptions:\n')
			sys.stderr.write('\t\t-h header\n')
			sys.stderr.write('\t\t-n number_of_reads\n')
			sys.stderr.write('\t\t-c coverage\n')
			sys.stderr.write('\t\t-t type (fasta or fatsq, fastq is default)\n')
			sys.stderr.write('\n')
			exit(1)

		profile_path = sys.argv[2]
		reference_path = sys.argv[3]
		reads_file = sys.argv[4]
		output_file = sys.argv[5]

		# Look at extension of the output file to decide whether to generate reads in fasta or fastq format
		# Format can also be specified by parameter -t. Parameter is stronger then file extension.
		# If nothing is specified, fastq reads are generated
		fastq = True
		if output_file.endswith('.fa') or output_file.endswith('.fasta'):
			fastq = False

		# parsing options
		header = numreads = coverage = None
		for i in range(5, len(sys.argv), 2):
			if sys.argv[i] == '-h':
				header = sys.argv[i+1]
			elif sys.argv[i] == '-n':
				try:
					numreads = int(sys.argv[i+1])
				except:
					sys.stderr.write('Wrong value for parameter "number of reads"! Value needs to be an int.\n\n')
					exit()
			elif sys.argv[i] == '-c':
				try:
					coverage = int(sys.argv[i+1])
				except:
					sys.stderr.write('Wrong value for parameter "coverage"! Value needs to be an int.\n\n')
					exit()
			elif sys.argv[i] == '-t':
				ftype = sys.argv[i+1].lower()
				if ftype == 'fastq':
					fastq = True
				elif ftype == 'fasta':
					fastq = False
				else:
					sys.stderr.write('Parameter -t goes with \'fasta\' or \'fastq\'.\n\n')
					exit()
			else:
				sys.stderr.write('Invalid parameter!.\n\n')
				exit()

		if numreads is None and coverage is None:
			sys.stderr.write('Either "coverage" or "number of reads" parameter must be specified.\n\n')
			exit()

		spike_with_cprofile(profile_path, reference_path, reads_file, output_file, numreads, coverage, header, fastq)

	elif (mode == 'combine_files'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Combines reads from two fasta or fastq files.\n')
			sys.stderr.write('\n')
			sys.stderr.write('Usage:\n')
			sys.stderr.write('\t%s %s <file1> <file2> <output file>' % (sys.argv[0], sys.argv[1]))
			sys.stderr.write('\n')
			exit(1)

		path1 = sys.argv[2]
		path2 = sys.argv[3]
		resultspath = sys.argv[4]
		combine_files(path1, path2, resultspath)

	elif (mode == 'combine_profiles'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Combines two profiles.\n')
			sys.stderr.write('\n')
			sys.stderr.write('Usage:\n')
			sys.stderr.write('\t%s %s <profile1> <profile2> <combined profile>' % (sys.argv[0], sys.argv[1]))
			sys.stderr.write('\n')
			exit(1)

		profile1 = sys.argv[2]
		profile2 = sys.argv[3]
		combined_profile = sys.argv[4]

		# TODO: implement combining two CIGAR profiles
		sys.stderr.write('\n\nCombining profiles is not yet implemented!\n')

	elif (mode == 'list_profiles'):
		if (len(sys.argv) != 2):
			sys.stderr.write('Lists available profiles.\n')
			sys.stderr.write('Requires no additional parameters to run.\n')
			sys.stderr.write('\n')
			exit(1)

		list_profiles()

	elif (mode == 'statistics_from_SAMs'):
		if (len(sys.argv) != 6):
			sys.stderr.write('Calculates some relevant statistics from SAM files.\n')
			sys.stderr.write('\n')
			sys.stderr.write('Usage:\n')
			sys.stderr.write('\t%s %s <reference> <BWA_SAM> <LAST_SAM> technology' % (sys.argv[0], sys.argv[1]))
			sys.stderr.write('\n')
			sys.stderr.write('\ttechnology - "illumina", "pacbio" or "nanopore"\n')
			sys.stderr.write('\n')
			exit(1)

		reference_path = sys.argv[2]
		sam_path_bwamem = sys.argv[3]
		sam_path_lastal = sys.argv[4]
		technology = sys.argv[5]
		if not (technology in ['illumina', 'pacbio', 'nanopore']):
			sys.stderr.write('Unsuported sequencing technology!\n\n')
			verbose_usage_and_exit()
		calculate_statistics(reference_path, sam_path_bwamem, sam_path_lastal)

	elif (mode == 'statistics_for_CProfile'):
		if (len(sys.argv) != 3):
			sys.stderr.write('Calculates relevant statistics for a CProfile file.\n')
			sys.stderr.write('\n')
			sys.stderr.write('Usage:\n')
			sys.stderr.write('\t%s %s <CIGAR profile>' % (sys.argv[0], sys.argv[1]))
			sys.stderr.write('\n')
			exit(1)

		profile_path = sys.argv[2]
		calculate_statistics_for_CProfile(profile_path)

	else:
		sys.stderr.write('Unsupported mode parameter.\n\n')
		verbose_usage_and_exit()
