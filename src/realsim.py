#! /usr/bin/python

import sys;
import setup_realsim;



def extract_error_profile(reference_path, reads_path, technology):
	sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name));
	exit(1);
	pass;

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
			sys.stderr.write('\t%s %s <reference> <reads> <technology>' % (sys.argv[0], sys.argv[1]));
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
			sys.stderr.write('\t%s %s <profiles> technology <reference> coverage' % (sys.argv[0], sys.argv[1]));
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
