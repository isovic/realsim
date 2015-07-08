#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
sys.path.append(SCRIPT_PATH + '/../src');

import subprocess;
import multiprocessing;

import basicdefines;

ALIGNER_URL = 'http://sourceforge.net/projects/bbmap/files/BBMap_35.10.tar.gz'
ALIGNER_PATH = os.path.join(basicdefines.ALIGNERS_PATH_ROOT_ABS, 'bbmap')
ALIGNER_BIN_SHORT = os.path.join(ALIGNER_PATH, 'bbmap.sh')
ALIGNER_BIN_LONG = os.path.join(ALIGNER_PATH, 'mapPacBio.sh')
ALIGNER_NAME = 'BBMAP'

ZIP_FILE = 'BBMap_35.10.tar.gz'
ZIP_PATH = os.path.join(basicdefines.ALIGNERS_PATH_ROOT_ABS, ZIP_FILE)



# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#    reads_file            Path to a FASTA/FASTQ file containing reads.
#    reference_file        Path to a reference genome FASTA file.
#    machine_name        A symbolic name to specify a set of parameters for a specific sequencing platform.
#    output_path            Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#    output_suffix        A custom suffix that can be added to the output filename.
def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
    parameters = ''
    num_threads = multiprocessing.cpu_count() / 2

    # This is intended to be used to set a parameter for max read length
    # However, this currently seems not to be working in BBMap
    if reads_file.endswith('.fa') or reads_file.endswith('.fasta'):
        type = 'fasta'
    elif reads_file.endswith('.fq') or reads_file.endswith('.fastq'):
        type = 'fastq'
    else:
        type = 'unknown'


    if ((machine_name.lower() == 'illumina') or (machine_name.lower() == 'roche')):
        parameters = '-t %s' % str(num_threads)
        aligner_bin = ALIGNER_BIN_SHORT
    elif ((machine_name.lower() == 'pacbio')):
        parameters = '-t %s -x pacbio' % str(num_threads)
        aligner_bin = ALIGNER_BIN_LONG
    elif ((machine_name.lower() == 'nanopore')):
        parameters = '-t %s -x ont2d' % str(num_threads);
        aligner_bin = ALIGNER_BIN_LONG
    elif ((machine_name.lower() == 'debug')):
        parameters = '-t %s' % str(num_threads);
        aligner_bin = ALIGNER_BIN_SHORT
    else:            # default
        parameters = '-t %s' % str(num_threads);
        aligner_bin = ALIGNER_BIN_SHORT


    if (output_suffix != ''):
        output_filename = '%s-%s' % (ALIGNER_NAME, output_suffix)
    else:
        output_filename = ALIGNER_NAME

    reads_basename = os.path.splitext(os.path.basename(reads_file))[0]
    sam_file = '%s/%s.sam' % (output_path, output_filename)
    memtime_file = '%s/%s.memtime' % (output_path, output_filename)
    memtime_file_index = '%s/%s-index.memtime' % (output_path, output_filename)

    # Index reference file, and measure execution time and memory.
    sys.stderr.write('[%s wrapper] Generating index...\n' % (ALIGNER_NAME))

    # command = '%s %s ref=%s' % (basicdefines.measure_command(memtime_file_index), ALIGNER_BIN, reference_file)
    # Atm doing it without measurements
    command = '%s ref=%s' % (ALIGNER_BIN_SHORT, reference_file)

    sys.stderr.write('[%s wrapper] %s\n' % (ALIGNER_NAME, command))
    subprocess.call(command, shell=True);

    # Run the alignment process, and measure execution time and memory.
    sys.stderr.write('[%s wrapper] Running %s...\n' % (ALIGNER_NAME, ALIGNER_NAME))

    # command = '%s %s in=%s out=%s' % (basicdefines.measure_command(memtime_file), ALIGNER_BIN, reads_file, sam_file)
    # Atm doing it without measurements
    command = '%s in=%s out=%s' % (aligner_bin, reads_file, sam_file)

    sys.stderr.write('[%s wrapper] %s\n' % (ALIGNER_NAME, command))
    subprocess.call(command, shell=True)
    sys.stderr.write('\n\n')

    sys.stderr.write('[%s wrapper] %s wrapper script finished processing.\n' % (ALIGNER_NAME, ALIGNER_NAME))

    return sam_file


# This is a standard interface for setting up the aligner. It should assume that the aligner
# is not present localy, but needs to be retrieved, unpacked, compiled and set-up, without requireing
# root privileges.
def download_and_install():
    if os.path.exists(ALIGNER_BIN_SHORT):
        sys.stderr.write('[%s wrapper] Bin found at %s. Skipping installation ...\n' % (ALIGNER_NAME, ALIGNER_BIN_SHORT))
    else:
        sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (ALIGNER_NAME, ALIGNER_NAME));

        if not os.path.exists(ZIP_PATH):
            sys.stderr.write('[%s wrapper] Downloading zip...\n' % (ALIGNER_NAME))
            command = 'cd %s; wget %s' % (basicdefines.ALIGNERS_PATH_ROOT_ABS, ALIGNER_URL)
            sys.stderr.write('[%s wrapper] %s\n' % (ALIGNER_NAME, command))
            subprocess.call(command, shell='True')

        # Decompress
        command = 'cd %s; tar -xzf %s' % (basicdefines.ALIGNERS_PATH_ROOT_ABS, ZIP_FILE)
        sys.stderr.write('[%s wrapper] %s\n' % (ALIGNER_NAME, command))
        subprocess.call(command, shell='True')

        # make
        # BBMap is precompiled and doesn't require make!

        sys.stderr.write('[%s wrapper] All instalation steps finished.\n' % (ALIGNER_NAME));
        sys.stderr.write('\n');



def verbose_usage_and_exit():
    sys.stderr.write('Usage:\n');
    sys.stderr.write('\t%s mode [<reads_file> <reference_file> <machine_name> <output_path> [<output_suffix>]]\n' % sys.argv[0]);
    sys.stderr.write('\n');
    sys.stderr.write('\t- mode          - either "run" or "install". If "install" other parameters can be ommitted.\n');
    sys.stderr.write('\t- machine_name  - "illumina", "roche", "pacbio", "nanopore" or "default".\n');
    sys.stderr.write('\t- output_suffix - suffix for the output filename.\n');

    exit(0);

if __name__ == "__main__":
    if (len(sys.argv) < 2 or len(sys.argv) > 7):
        verbose_usage_and_exit();

    if (sys.argv[1] == 'install'):
        download_and_install();
        exit(0);

    elif (sys.argv[1] == 'run'):
        if (len(sys.argv) < 6):
            verbose_usage_and_exit();

        reads_file = sys.argv[2];
        reference_file = sys.argv[3];
        machine_name = sys.argv[4];
        output_path = sys.argv[5];
        output_suffix = '';

        if (len(sys.argv) == 7):
            output_suffix = sys.argv[6];
        run(reads_file, reference_file, machine_name, output_path, output_suffix);

    else:
        verbose_usage_and_exit();
