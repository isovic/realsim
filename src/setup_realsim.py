#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
import sys
sys.path.append(SCRIPT_PATH + '');
sys.path.append(SCRIPT_PATH + '/wrappers');
import subprocess;

import basicdefines;

def create_folders():
	sys.stderr.write('Generating folders...\n');

	if not os.path.exists(basicdefines.TOOLS_ROOT_ABS):
		sys.stderr.write('Creating folder "%s".\n' % basicdefines.TOOLS_ROOT_ABS);
		os.makedirs(basicdefines.TOOLS_ROOT_ABS);

	if not os.path.exists(basicdefines.INTERMEDIATE_PATH_ROOT_ABS):
		sys.stderr.write('Creating folder "%s".\n' % basicdefines.INTERMEDIATE_PATH_ROOT_ABS);
		os.makedirs(basicdefines.INTERMEDIATE_PATH_ROOT_ABS);

	if not os.path.exists(basicdefines.PROFILES_PATH_ABS):
		sys.stderr.write('Creating folder "%s".\n' % basicdefines.PROFILES_PATH_ABS);
		os.makedirs(basicdefines.PROFILES_PATH_ABS);

	sys.stderr.write('\n');

def download_aligners():
	sys.stderr.write('Installing alignment tools.\n');

	aligner_wrappers = basicdefines.find_files(basicdefines.WRAPPERS_PATH_ROOT_ABS, 'wrapper_*.py');

	for wrapper in aligner_wrappers:
		wrapper_basename = os.path.splitext(os.path.basename(wrapper))[0];
		command = 'import %s; %s.download_and_install()' % (wrapper_basename, wrapper_basename);
		exec(command);

def setup_tools():
	sys.stderr.write('Cloning Cgmemtime Git repo. Git needs to be installed.\n');
	command = 'cd %s; git clone https://github.com/isovic/cgmemtime.git' % (basicdefines.TOOLS_ROOT_ABS);
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

def setup_all():
	create_folders();
	setup_tools();
	download_aligners();



def verbose_usage_and_exit():
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s\n' % sys.argv[0]);
	sys.stderr.write('\n');
	exit(0);

if __name__ == '__main__':
	setup_all();

	# if (len(sys.argv) > 2):
	# 	verbose_usage_and_exit();
