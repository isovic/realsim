#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
import sys
import uuid;
sys.path.append(SCRIPT_PATH + '');
sys.path.append(SCRIPT_PATH + '/wrappers');

import setup_realsim;
import basicdefines;
import wrapper_lastal;
import wrapper_bwamem;
import analyze_venn_of_mappings;

# TODO:
#   1. Cleanup of intermediate files



def extract_error_profile(reference_path, reads_path, technology):
	if (os.path.exists(basicdefines.TOOLS_ROOT_ABS) == False or os.path.exists(basicdefines.ALIGNERS_PATH_ROOT_ABS) == False):
		sys.stderr.write('ERROR: Required tools not set up. Rerun the script with the "setup" option.\nExiting.\n\n');
		exit(1);

	output_path = '%s/../intermediate' % SCRIPT_PATH;
	uuid_string = str(uuid.uuid4());			# Generate a random UUID so that multiple runs don't clash.

	if (not os.path.exists(output_path)):
		sys.stderr.write('Creating folder "%s".\n' % (output_path));
		os.makedirs(output_path);

	sam_path_bwamem = wrapper_bwamem.run(reads_path, reference_path, technology, output_path, uuid_string);
	sam_path_lastal = wrapper_lastal.run(reads_path, reference_path, technology, output_path, uuid_string);

	distance_threshold = -1;			# This currently does nothing.
	# compare_two_sams returns three dicts:
	#	1. sam_hash1	-	hashed SAM file. Key of the dict is the qname of a SAM entry, and value is a list of SAMLine objects. List contains all alignments for the same qname
	#						and is composed of at least one element. The list is sorted in the descending order of the MAPQ parameter is available (i.e. if it is different from 255),
	#						otherwise it is sorted according to the AS value (again, if available).
	#	2. sam_hash2	-	same as the previous parameter, but for the second given SAM file.
	#	3. distance_to_qname_hash	-	Key of the dict is the distance between the mapping positions of the same read accross two SAM files. Value for each key contains
	#									a list of qnames at that particular distance. To retreieve concrete info on a specific alignment for the given qname, one can
	#									query the sam_hash1 and sam_hash2 objects.
	[sam_hash_bwamem, sam_hash_lastal, distance_to_qname_hash] = analyze_venn_of_mappings.compare_two_sams(sam_path_bwamem, sam_path_lastal, distance_threshold, out_summary_prefix='');

	# Example of how to access reads at various distances between two SAM files:
	sys.stdout.write('Histogram of alignment distances:\n');
	i = 0;
	for distance_key in sorted(distance_to_qname_hash.keys()):
		qnames = distance_to_qname_hash[distance_key];
		sys.stdout.write('[%d] Distance: %d, number of reads: %d\n' % (i, distance_key, len(qnames)));
		i += 1;
		if (distance_key > 10):
			break;

	# Further example, working but commented out not to clog the stdout.
	# try:
	# 	qnames = distance_to_qname_hash[0];
	# 	sys.stdout.write('SAM line from BWA-MEM mapping:\n');
	# 	sys.stdout.write('%s\n\n' % sam_hash_bwamem[qnames[0]][0].original_line);

	# 	sys.stdout.write('SAM line from LAST mapping:\n');
	# 	sys.stdout.write('%s\n\n' % sam_hash_lastal[qnames[0]][0].original_line);
	# except Exception, e:
	# 	sys.stderr.write('%s\n' % e);



def simulate(profiles_path, technology, reference_path, coverage):
	sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name));
	exit(1);
	pass;



def verbose_usage_and_exit():
	sys.stderr.write('RealSim - sequence simulator based on empirical observations.\n');
	sys.stderr.write('\n');
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s [mode]\n' % sys.argv[0]);
	sys.stderr.write('\n');
	sys.stderr.write('\tmode:\n');
	sys.stderr.write('\t\tsetup\n');
	sys.stderr.write('\t\textract\n');
	sys.stderr.write('\t\tsimulate\n');
	sys.stderr.write('\n');
	exit(0);

if __name__ == '__main__':
	if (len(sys.argv) < 2):
		verbose_usage_and_exit();

	mode = sys.argv[1];

	if (mode == 'setup'):
		if (len(sys.argv) != 2):
			sys.stderr.write('Setup the folder structures and install necessary tools.\n');
			sys.stderr.write('Requires no additional parameters to run.\n');
			sys.stderr.write('\n');
			exit(1);
		setup_realsim.setup_all();

	elif (mode == 'extract'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Extracts the error profile from a given set of reads and a reference genome.\n');
			sys.stderr.write('Error profile is output to stdout in a tab delimited format.\n');
			sys.stderr.write('\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <reference> <reads> technology' % (sys.argv[0], sys.argv[1]));
			sys.stderr.write('\n');
			sys.stderr.write('\ttechnology - "illumina", "pacbio" or "nanopore"\n');
			sys.stderr.write('\n');
			exit(1);

		reference_path = sys.argv[2];
		reads_path = sys.argv[3];
		technology = sys.argv[4];
		if ((technology in ['illumina', 'pacbio', 'nanopore']) == False):
			sys.stderr.write('Unsuported sequencing technology!\n\n');
			verbose_usage_and_exit();
		extract_error_profile(reference_path, reads_path, technology);

	elif (mode == 'simulate'):
		if (len(sys.argv) != 6):
			sys.stderr.write('Applies an extracted error profile on a given reference in order to simulate reads.\n');
			sys.stderr.write('Output is generated in two files: a FASTQ file and a SAM file.\n');
			sys.stderr.write('\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <profiles> technology <reference> coverage\n' % (sys.argv[0], sys.argv[1]));
			sys.stderr.write('\n');
			sys.stderr.write('\ttechnology - "illumina", "pacbio" or "nanopore"\n');
			sys.stderr.write('\tcoverage - a numeric (int) value for sequencing depth.\n');
			sys.stderr.write('\n');
			exit(1);

		profiles_path = sys.argv[2];
		technology = sys.argv[3];
		reference_path = sys.argv[4];
		if ((technology in ['illumina', 'pacbio', 'nanopore']) == False):
			sys.stderr.write('Unsuported sequencing technology!\n\n');
			verbose_usage_and_exit();
		try:
			coverage = int(sys.argv[5]);
		except:
			sys.stderr.write('Wrong value for parameter "coverage"! Value needs to be an int.\n\n');
			verbose_usage_and_exit();
		simulate(profiles_path, technology, reference_path, coverage);

	else:
		sys.stderr.write('Unsupported mode parameter.\n\n');
		verbose_usage_and_exit();
