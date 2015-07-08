#! /usr/bin/python

# The MIT License (MIT)

# Copyright (c) 2015 Ivan Sovic

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
sys.path.append(SCRIPT_PATH + '/../src');

import subprocess;
import multiprocessing;

import basicdefines;
import fastqparser;

ALIGNER_URL = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz';
ALIGNER_PATH = os.path.join(basicdefines.ALIGNERS_PATH_ROOT_ABS, 'ncbi-blast-2.2.30+/bin/')
BIN = 'blastn';
MAPPER_NAME = 'BLAST';

# This specifies the output format of the BLAST alignments, for tabulated output.
outfmt = 'qseqid qlen qstart qend sseqid slen sstart send sstrand evalue bitscore score length btop qseq sseq';



# TODO:
# Fix this error:
	# [BLAST wrapper] BLAST wrapper script finished processing.
	# [BLAST wrapper] Converting BLAST output to SAM file...
	# Current BLAST output line: 21903400Traceback (most recent call last):
	#   File "./run-alignment.py", line 80, in <module>
	#     wrapper_blast.run(reads_fasta, reference_file, machine_names[machine_num], output_path, machine_suffix);
	#   File "aligneval/wrappers/wrapper_blast.py", line 324, in run
	#     convert_blast_to_sam(reference_file, reads_file, out_file, sam_file);
	#   File "aligneval/wrappers/wrapper_blast.py", line 218, in convert_blast_to_sam
	#     sam_cigar = convert_btop_to_cigar(btop, num_clip_front, num_clip_back, sstrand);
	#   File "aligneval/wrappers/wrapper_blast.py", line 67, in convert_btop_to_cigar
	#     if (match_start >= 0 and match_end >= 0):
	# UnboundLocalError: local variable 'match_end' referenced before assignment

def convert_btop_to_cigar(btop, num_clip_front, num_clip_back, sstrand, use_extended_cigar=False):
	cigar_ops = [];

	# Add the front clipping.
	if (num_clip_front > 0):
		cigar_ops.append([num_clip_front, 'S']);

	# Convert the BTOP format into the CIGAR format.
	match_start = 0;
	match_end = -1;
	i = 0;
	while (i < len(btop)):
		is_digit = btop[i].isdigit();
		if (is_digit == True):
			if (match_start == -1):
				match_start = i;
			match_end = i;
		elif (is_digit == False):
			# Handle the match events.
			if (match_start >= 0 and match_end >= 0):
				if (use_extended_cigar == True):
					cigar_ops.append([int(btop[match_start:(match_end + 1)]), '=']);
				else:
					cigar_ops.append([int(btop[match_start:(match_end + 1)]), 'M']);
				match_start = -1;	match_end = -1;

			# Mismatch/insertion/deletion operations should be two-character events.
			# Check if we have enough characters left.
			if ((i + 1) >= len(btop)):
				sys.stderr.write('ERROR: BTOP string deformed!\n');
				exit(1);

			# Handle the mismatches, insertions and deletions.
			current_error = [];
			if (btop[i] != '-' and btop[i+1] != '-'):
				if (use_extended_cigar == True):
					current_error = [1, 'X'];
				else:
					current_error = [1, 'M'];
			elif (btop[i] == '-' and btop[i+1] != '-'):
				current_error = [1, 'D'];
			elif (btop[i] != '-' and btop[i+1] == '-'):
				current_error = [1, 'I'];
			else:
				sys.stderr.write('ERROR: BTOP string deformed! Invalid combination of characters "--".\n');
				exit(1);

			# Skip the next character.
			i += 1;

			if (len(cigar_ops) > 0 and cigar_ops[-1][1] == current_error[1]):
				cigar_ops[-1][0] += current_error[0];
			else:
				cigar_ops.append(current_error);

		i += 1;

	# Append the last CIGAR operation.
	if (match_start >= 0 and match_end >= 0):
		current_cigar = [int(btop[match_start:(match_end + 1)]), '='] if (use_extended_cigar == True) else [int(btop[match_start:(match_end + 1)]), 'M'];

		if (len(cigar_ops) > 0 and cigar_ops[-1][1] == current_cigar[1]):
			cigar_ops[-1][0] += current_cigar[0];
		else:
			cigar_ops.append(current_cigar);

	# Add the back clipping.
	if (num_clip_back > 0):
		cigar_ops.append([num_clip_back, 'S']);

	# Create the final CIGAR string.
	cigar = '';
	if (sstrand == 'plus'):
		cigar = ''.join(['%d%s' % (cigar_op[0], cigar_op[1]) for cigar_op in cigar_ops]);
	else:
		cigar = '';
		i = len(cigar_ops) - 1;
		while (i >= 0):
			cigar += '%d%s' % (cigar_ops[i][0], cigar_ops[i][1]);
			i -= 1;

	if (cigar == ''):
		cigar = '*';

	return cigar;



def convert_blast_to_sam(reference_file, reads_file, blast_out_file, sam_file):
	sys.stderr.write('[%s wrapper] Converting BLAST output to SAM file...\n' % (MAPPER_NAME));

	[ref_headers, ref_seqs, ref_quals] = fastqparser.read_fastq(reference_file);
	ref_header_hash = {};
	i = 0;
	while (i < len(ref_headers)):
		ref_header_hash[ref_headers[i]] = i;		ref_header_hash[ref_headers[i].split()[0]] = i;		i += 1;

	[read_headers, read_seqs, read_quals] = fastqparser.read_fastq(reads_file);
	read_header_hash = {};
	i = 0;
	while (i < len(read_headers)):
		read_header_hash[read_headers[i]] = i;		read_header_hash[read_headers[i].split()[0]] = i;		i += 1;



	try:
		fp_in = open(blast_out_file, 'r');
	except Exception, e:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % blast_out_file);
		exit(1);
	# blast_lines = fp.readlines();

	try:
		fp_out = open(sam_file, 'w');
	except Exception, e:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % sam_file);
		exit(1);

	# Write the SAM header.
	i = 0;
	while i < len(ref_headers):
		line = '@SQ\tSN:%s\tLN:%d\n' % (ref_headers[i], len(ref_seqs[i]));
		fp_out.write(line);
		i += 1;

	# sam_lines = [];

	current_line = 0;
	for line in fp_in:
		current_line += 1;
		if ((current_line % 100) == 0):
			sys.stderr.write('\rCurrent BLAST output line: %d' % current_line);

		line = line.strip();
		if (len(line) == 0):
			continue;
		split_line = line.split('\t');
		split_outfmt = outfmt.split();
		if (len(split_line) != len(split_outfmt)):
			continue;
		command = '[%s] = split_line;' % (', '.join(split_outfmt));
		exec(command);

		flag = 0;
		flag |= (0 if (sstrand == 'plus') else (1 << 4));

		pos = 0;
		cigar = '*';
		try:
			seq = read_seqs[read_header_hash[qseqid]];
		except:
			# In this case, something weird happened. Most likely the header got messed up.
			# Another option is that someone changed the reads file. In any case, if the original read
			# cannot be found, we will call this alignment unmapped.
			seq = '*';
			flag = 4;

		qual = '*';
		if (len(read_quals) > 0):
			try:
				qual = read_quals[read_header_hash[qseqid]];
			except:
				qual = '*';

		sam_start = (int(sstart) - 1) if (sstrand == 'plus') else (int(send) - 1);
		sam_end = (int(send) - 1) if (sstrand == 'plus') else (int(sstart) - 1);
		sam_seq = (seq) if (sstrand == 'plus' or seq == '*') else (fastqparser.revcomp_seq(seq));			# Reverse the seq field if necessary.
		sam_qual = (qual) if (sstrand == 'plus' or qual == '*') else (qual[::-1]);							# Reverse the quality values if necessary.
		num_clip_front = int(qstart) - 1;
		num_clip_back = int(qlen) - (int(qend));
		sam_cigar = convert_btop_to_cigar(btop, num_clip_front, num_clip_back, sstrand);

		sam_line = '';
		sam_line += '%s\t' % (qseqid);						# 1. qname
		sam_line += '%d\t' % (flag);						# 2. flag
		sam_line += '%s\t' % (sseqid);						# 3. rname
		sam_line += '%d\t' % (sam_start + 1);				# 4. pos
		sam_line += '255\t';								# 5. mapq
		sam_line += '%s\t' % (sam_cigar);					# 6. CIGAR
		sam_line += '*\t';									# 7. rnext
		sam_line += '0\t';									# 8. pnext
		sam_line += '0\t';									# 9. tlen
		sam_line += '%s\t' % (sam_seq);						# 10. seq
		sam_line += '%s\t' % (sam_qual);					# 11. qual
		sam_line += 'AS:i:%s\t' % (score.strip());			# AS, custom
		sam_line += 'ZE:f:%s\t' % (evalue.strip());			# custom, evalue
		sam_line += 'ZB:i:%s\t' % (bitscore.strip());		# custom, bitscore
		sam_line += 'ZA:i:%s' % (length.strip());			# custom, alignment length

		# Write the SAM lines.
		fp_out.write(sam_line + '\n');

	fp_in.close();
	fp_out.close();

	sys.stderr.write('\n');



# Function 'run' should provide a standard interface for running a mapper. Given input parameters, it should run the
# alignment process, and convert any custom output results to the SAM format. Function should return a string with the
# path to the output file.
#	reads_file			Path to a FASTA/FASTQ file containing reads.
#	reference_file		Path to a reference genome FASTA file.
#	machine_name		A symbolic name to specify a set of parameters for a specific sequencing platform.
#	output_path			Folder to which the output will be placed to. Filename will be automatically generated according to the name of the mapper being run.
#	output_suffix		A custom suffix that can be added to the output filename.
def run(reads_file, reference_file, machine_name, output_path, output_suffix=''):
	parameters = '';
	num_threads = multiprocessing.cpu_count() / 2;

	if ((machine_name.lower() == 'illumina') or (machine_name.lower() == 'roche')):
		parameters = '-num_threads %s' % str(num_threads);

	elif ((machine_name.lower() == 'pacbio')):
		parameters = '-reward 5 -penalty -4 -gapopen 8 -gapextend 6 -dust no';
		parameters += ' -num_threads %s' % str(num_threads);

	elif ((machine_name.lower() == 'nanopore')):
		# These parameters used in the paper: ""
		# These parameters used in the paper: "Oxford Nanopore Sequencing and de novo Assembly of a Eukaryotic Genome", Supplemental Notes and Figures
		# http://biorxiv.org/content/biorxiv/suppl/2015/01/06/013490.DC1/013490-1.pdf
		# Quote: "Overall accuracy was calculated by aligning the raw Oxford Nanopore reads to the W303 pacbio assembly using Blast version 2.2.27+ with the following parameters:"
		# parameters += ' -reward 5 -penalty -4 -gapopen 8 -gapextend 6 -dust no -evalue 1e-10';
		parameters = '-reward 5 -penalty -4 -gapopen 8 -gapextend 6 -dust no';
		parameters += ' -num_threads %s' % str(num_threads);

	elif ((machine_name.lower() == 'debug')):
		parameters = '-num_threads %s' % str(num_threads);
		# sys.stderr.write('ERROR: Debug parameters not implemented yet!\n');
		# exit(1);

	else:			# default
		parameters = '-num_threads %s' % str(num_threads);



	# http://www.kenkraaijeveld.nl/genomics/bioinformatics/
	# The first thing to do is to build your contig.fa file into a Blast database. Type:
	# $ makeblastdb -in [path to contigs.fa] -dbtype nucl -out [path to output directory]
	# You can now query this database with sequences that you want to find. For example:
	# $ blastn -query [path to file with sequence of interest] -task blastn -db [path to your database] -out [path to output directory] -num_threads 8

	if (output_suffix != ''):
		output_filename = '%s-%s' % (MAPPER_NAME, output_suffix);
	else:
		output_filename = MAPPER_NAME;

	reads_basename = os.path.splitext(os.path.basename(reads_file))[0];
	sam_file = '%s/%s.sam' % (output_path, output_filename);
	out_file = '%s/%s.out' % (output_path, output_filename);
	filtered_out_file = '%s/%s-filtered.out' % (output_path, output_filename);
	out_db_path = '%s-blastdb' % (reference_file);
	memtime_file = '%s/%s.memtime' % (output_path, output_filename);
	memtime_file_index = '%s/%s-index.memtime' % (output_path, output_filename);

	# Run the indexing process, and measure execution time and memory.
	if (not os.path.exists(out_db_path + '.nsq')):
		sys.stderr.write('[%s wrapper] Generating index...\n' % (MAPPER_NAME));

		# command = '%s %s/makeblastdb -in %s -dbtype nucl -out %s' % (basicdefines.measure_command(memtime_file_index), ALIGNER_PATH, reference_file, out_db_path);
		# ATM doiing it withot measurements
		command = '%s/makeblastdb -in %s -dbtype nucl -out %s' % (ALIGNER_PATH, reference_file, out_db_path);

		sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
		# subprocess.call(command, shell=True);
		sys.stderr.write('\n\n');
	else:
		sys.stderr.write('[%s wrapper] Reference index already exists. Continuing.\n' % (MAPPER_NAME));
		sys.stderr.flush();

	# Run the alignment process, and measure execution time and memory.
	sys.stderr.write('[%s wrapper] Running %s...\n' % (MAPPER_NAME, MAPPER_NAME));
	# command = '%s %s/%s -task blastn -db %s -query %s -out %s %s' % (basicdefines.measure_command(memtime_file), ALIGNER_PATH, BIN, out_db_path, reads_file, out_file, parameters);

	# command = '%s %s/%s -task blastn -db %s -query %s -out %s %s -outfmt "6 %s"' % (basicdefines.measure_command(memtime_file), ALIGNER_PATH, BIN, out_db_path, reads_file, out_file, parameters, outfmt);
	# ATM doiing it without measurements
	command = '%s/%s -task blastn -db %s -query %s -out %s %s -outfmt "6 %s"' % (ALIGNER_PATH, BIN, out_db_path, reads_file, out_file, parameters, outfmt)

	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	# subprocess.call(command, shell=True);
	sys.stderr.write('\n\n');

	# Filter the BLAST out file and extract only one alignment per read (the one with highest alignment score).
	sys.stderr.write('[%s wrapper] Filtering BLAST output...\n' % (MAPPER_NAME));
	command = '%s/filterblastout/bin/filterblastout %s > %s' % (SCRIPT_PATH, out_file, filtered_out_file);
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n\n');

	convert_blast_to_sam(reference_file, reads_file, filtered_out_file, sam_file);

	sys.stderr.write('[%s wrapper] %s wrapper script finished processing.\n' % (MAPPER_NAME, MAPPER_NAME));

	return sam_file



# This is a standard interface for setting up the aligner. It should assume that the aligner
# is not present localy, but needs to be retrieved, unpacked, compiled and set-up, without requireing
# root privileges.
def download_and_install():
	sys.stderr.write('[%s wrapper] Started installation of %s.\n' % (MAPPER_NAME, MAPPER_NAME));
	sys.stderr.write('[%s wrapper] Downloading the package.\n' % (MAPPER_NAME));
	command = 'cd %s; wget %s; tar -xzvf %s' % (basicdefines.ALIGNERS_PATH_ROOT_ABS, ALIGNER_URL, os.path.basename(ALIGNER_URL));
	subprocess.call(command, shell='True');
	sys.stderr.write('\n');

	# No make required - these are precompiled binaries.

	sys.stderr.write('[%s wrapper] All instalation steps finished.\n' % (MAPPER_NAME));
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



# http://blastedbio.blogspot.sg/2012/05/blast-tabular-missing-descriptions.html

#  -outfmt
#    alignment view options:
#      0 = pairwise,
#      1 = query-anchored showing identities,
#      2 = query-anchored no identities,
#      3 = flat query-anchored, show identities,
#      4 = flat query-anchored, no identities,
#      5 = XML Blast output,
#      6 = tabular,
#      7 = tabular with comment lines,
#      8 = Text ASN.1,
#      9 = Binary ASN.1,
#     10 = Comma-separated values,
#     11 = BLAST archive format (ASN.1)

#    Options 6, 7, and 10 can be additionally configured to produce
#    a custom format specified by space delimited format specifiers.
#    The supported format specifiers are:
#         qseqid means Query Seq-id
#            qgi means Query GI
#           qacc means Query accesion
#        qaccver means Query accesion.version
#           qlen means Query sequence length
#         sseqid means Subject Seq-id
#      sallseqid means All subject Seq-id(s), separated by a ';'
#            sgi means Subject GI
#         sallgi means All subject GIs
#           sacc means Subject accession
#        saccver means Subject accession.version
#        sallacc means All subject accessions
#           slen means Subject sequence length
#         qstart means Start of alignment in query
#           qend means End of alignment in query
#         sstart means Start of alignment in subject
#           send means End of alignment in subject
#           qseq means Aligned part of query sequence
#           sseq means Aligned part of subject sequence
#         evalue means Expect value
#       bitscore means Bit score
#          score means Raw score
#         length means Alignment length
#         pident means Percentage of identical matches
#         nident means Number of identical matches
#       mismatch means Number of mismatches
#       positive means Number of positive-scoring matches
#        gapopen means Number of gap openings
#           gaps means Total number of gaps
#           ppos means Percentage of positive-scoring matches
#         frames means Query and subject frames separated by a '/'
#         qframe means Query frame
#         sframe means Subject frame
#           btop means Blast traceback operations (BTOP)
#    When not provided, the default value is:
#    'qseqid sseqid pident length mismatch gapopen qstart qend sstart send
#    evalue bitscore', which is equivalent to the keyword 'std'
#    Default = `0'



# $ blastp -query spo_like.faa -db nr -max_target_seqs 1 -outfmt "6 qseqid sseqid"
# spo_like gi|384161333|ref|YP_005543406.1|



# http://www.ncbi.nlm.nih.gov/books/NBK279682/
# The "Blast trace-back operations" (BTOP) string describes the alignment produced by BLAST.
# This string is similar to the CIGAR string produced in SAM format, but there are important
# differences. BTOP is a more flexible format that lists not only the aligned region but also
# matches and mismatches. BTOP operations consist of 1.) a number with a count of matching letters,
# 2.) two letters showing a mismatch (e.g., "AG" means A was replaced by G), or 3.) a dash ("-") and a
# letter showing a gap. The box below shows a blastn run first with BTOP output and then the same
# run with the BLAST report showing the alignments.






# This is from a different source, it also has a parameter for the strand:
# *** Formatting options
#  -outfmt <String>
#    alignment view options:
#      0 = pairwise,
#      1 = query-anchored showing identities,
#      2 = query-anchored no identities,
#      3 = flat query-anchored, show identities,
#      4 = flat query-anchored, no identities,
#      5 = XML Blast output,
#      6 = tabular,
#      7 = tabular with comment lines,
#      8 = Text ASN.1,
#      9 = Binary ASN.1,
#     10 = Comma-separated values,
#     11 = BLAST archive format (ASN.1)

   # Options 6, 7, and 10 can be additionally configured to produce
   # a custom format specified by space delimited format specifiers.
   # The supported format specifiers are:
   #         qseqid means Query Seq-id
   #            qgi means Query GI
   #           qacc means Query accesion
   #        qaccver means Query accesion.version
   #           qlen means Query sequence length
   #         sseqid means Subject Seq-id
   #      sallseqid means All subject Seq-id(s), separated by a ';'
   #            sgi means Subject GI
   #         sallgi means All subject GIs
   #           sacc means Subject accession
   #        saccver means Subject accession.version
   #        sallacc means All subject accessions
   #           slen means Subject sequence length
   #         qstart means Start of alignment in query
   #           qend means End of alignment in query
   #         sstart means Start of alignment in subject
   #           send means End of alignment in subject
   #           qseq means Aligned part of query sequence
   #           sseq means Aligned part of subject sequence
   #         evalue means Expect value
   #       bitscore means Bit score
   #          score means Raw score
   #         length means Alignment length
   #         pident means Percentage of identical matches
   #         nident means Number of identical matches
   #       mismatch means Number of mismatches
   #       positive means Number of positive-scoring matches
   #        gapopen means Number of gap openings
   #           gaps means Total number of gaps
   #           ppos means Percentage of positive-scoring matches
   #         frames means Query and subject frames separated by a '/'
   #         qframe means Query frame
   #         sframe means Subject frame
   #           btop means Blast traceback operations (BTOP)
   #        staxids means Subject Taxonomy ID(s), separated by a ';'
   #      sscinames means Subject Scientific Name(s), separated by a ';'
   #      scomnames means Subject Common Name(s), separated by a ';'
   #     sblastnames means Subject Blast Name(s), separated by a ';'
   #              (in alphabetical order)
   #     sskingdoms means Subject Super Kingdom(s), separated by a ';'
   #              (in alphabetical order)
   #         stitle means Subject Title
   #     salltitles means All Subject Title(s), separated by a '<>'
   #        sstrand means Subject Strand
   #          qcovs means Query Coverage Per Subject
   #        qcovhsp means Query Coverage Per HSP
