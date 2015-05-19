#! /usr/bin/python

# Written by Ivan Sovic, March 2015.
#
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



import os;
import sys;
import utility_sam;
import fastqparser


def combine_msa_output(combinedoutputpath, msaoutputpath, sam_hash_bwamem, sam_hash_lastal, windowstart=0):

	with open(combinedoutputpath, 'w') as combinedoutputfile:
		msaoutput = fastqparser.read_fastq(msaoutputpath)
		headers = msaoutput[0]
		seqs = msaoutput[1]

		# if there is nothing to write just exit
		if len(headers) == 0:
			sys.stderr.out('\n\nNo output to combine!\n')
			return


		# first line (index 0) is for reference
		combinedoutputfile.write(headers[0] + '\n')
		combinedoutputfile.write(seqs[0] + '\n')

		# noting where MSA put inserts into reference sequence
		# each field in list refinserts contains the number of inserts in
		# reference sequence up to and including that point
		numins = 0
		refinserts = [0] * len(seqs[0])
		for i in xrange(len(refinserts)):
			if seqs[0][i] == '-':
				numins += 1
			refinserts[i] = numins


		for i in xrange(1, len(headers)):

			# writing sequence and header from MAFFT output
			combinedoutputfile.write(headers[i] + '\n')
			combinedoutputfile.write(seqs[i] + '\n')

			# getting sequence data from BWA SAM
			samline = sam_hash_bwamem[headers[i]][0]
			combinedoutputfile.write('BWA_' + headers[i] + '\n')
			# prefix = '-'*samline.clipped_pos							# to correctly position sequence from SAM
			# clipped pos didnt seem to work, so trying just pos
			# taking into account that reference sequence can have inserts, -1 is because pos is one-based
			# BWA uses soft clipping so I'm using clipped_pos
			prefix = '-'*(samline.clipped_pos - windowstart + refinserts[samline.clipped_pos - windowstart] - 1)
			combinedoutputfile.write(prefix + samline.seq + '\n')
			combinedoutputfile.write(prefix + samline.cigar + '\n')

			# getting sequence data from LAST SAM
			samline = sam_hash_lastal[headers[i]][0]
			combinedoutputfile.write('LAST_' + headers[i] + '\n')
			# prefix = '-'*samline.clipped_pos							# to correctly position sequence from SAM
			# clipped pos didnt seem to work, so trying just pos
			# taking into account that reference sequence can have inserts, -1 is because pos is one-based
			# LAST uses hard clipping so I'm using samline.pos
			prefix = '-'*(samline.pos - windowstart + refinserts[samline.pos - windowstart] - 1)
			combinedoutputfile.write(prefix + samline.seq + '\n')
			combinedoutputfile.write(prefix + samline.cigar + '\n')


def compare_two_sams(sam_file1, sam_file2, distance_threshold, out_summary_prefix=''):

	# sys.stderr.write('Loading first SAM file...');
	qnames_with_multiple_alignments = {};
	# print sam_file1, sam_file2, distance_threshold, out_summary_prefix
	sys.stderr.write('Loading the first SAM file into hash...\n');
	[sam_hash1, sam_hash1_num_lines, sam_hash1_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_file1, qnames_with_multiple_alignments);
	sys.stderr.write('Loading the second SAM file into hash...\n');
	[sam_hash2, sam_hash2_num_lines, sam_hash2_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_file2, qnames_with_multiple_alignments);

	not_in_sam_file1 = 0;
	not_in_sam_file2 = 0;

	num_different_reference = 0;
	num_different_orientation = 0;
	num_not_mapped_1 = 0;
	num_not_mapped_2 = 0;
	num_mapped_1 = 0;
	num_mapped_2 = 0;

	qname_to_distance_hash = {};
	distance_count_hash = {};
	distance_to_qname_hash = {};
	distance_to_sam_hash = {};

	num_processed = 0;

	for qname in sam_hash1.iterkeys():
		num_processed += 1;
		sys.stderr.write('\rProcessed %d alignments...' % num_processed);

		if (len(sam_hash1[qname]) > 0 and sam_hash1[qname][0].IsMapped() == True):
			num_mapped_1 += 1;

		# TODO: THIS NEEDS TO BE REMOVED OR IMPLEMENTED SOMEHOW DIFFERENTLY!!
		# The point of this was that, BLASR doesn't conform to the SAM standard, and makes it difficult to
		# uniformly evaluate the results!
		# if 'blasr' in sam_file1.lower():
		# 	qname = '/'.join(qname.split('/')[:-1]);

		sam_line_list1 = sam_hash1[qname];

		sam_line_list_2 = [];
		try:
			sam_line_list2 = sam_hash2[qname];
		except:
			not_in_sam_file2 += 1;
			continue;

		sorted_sam_line_list1 = sorted(sam_line_list1, key=lambda sam_line: (-sam_line.chosen_quality));
		sorted_sam_line_list2 = sorted(sam_line_list2, key=lambda sam_line: (-sam_line.chosen_quality));
		if (len(sorted_sam_line_list1) > 0 and len(sorted_sam_line_list2) > 0):

			if (sorted_sam_line_list1[0].IsMapped() == False):
				num_not_mapped_1 += 1;
			if (sorted_sam_line_list2[0].IsMapped() == False):
				num_not_mapped_2 += 1;
			if (sorted_sam_line_list1[0].IsMapped() == False or sorted_sam_line_list2[0].IsMapped() == False):
				continue;

			if (not ((sorted_sam_line_list1[0].rname in sorted_sam_line_list2[0].rname) or (sorted_sam_line_list2[0].rname in sorted_sam_line_list1[0].rname))):
				num_different_reference += 1;
				continue;
			if (sorted_sam_line_list1[0].IsReverse() != sorted_sam_line_list2[0].IsReverse()):
				num_different_orientation += 1;
				continue;

			distance = abs(sorted_sam_line_list1[0].clipped_pos - sorted_sam_line_list2[0].clipped_pos);
			qname_to_distance_hash[qname] = distance;
			if (distance in distance_count_hash):
				distance_count_hash[distance] += 1;
				distance_to_qname_hash[distance].append(qname);
				distance_to_sam_hash[distance].append(sorted_sam_line_list1[0]);
			else:
				distance_count_hash[distance] = 1;
				distance_to_qname_hash[distance] = [qname];
				distance_to_sam_hash[distance] = [sorted_sam_line_list1[0]];
		else:
			if (len(sorted_sam_line_list1) == 0):
				not_in_sam_file1 += 1;
			if (len(sorted_sam_line_list2) == 0):
				not_in_sam_file2 += 1;
			continue;

		# min_distance = -1;
		# i = 0;
		# for sam_line1 in sorted_sam_line_list1:
		# 	for sam_line2 in sorted_sam_line_list2:
		# 		distance = abs(sam_line1.clipped_pos - sam_line2.clipped_pos);
		# 		if (i == 0 or distance < min_distance):
		# 			min_distance = distance;
		# 		i += 1;
		# distance_hash[qname] = min_distance;



	sys.stderr.write('\n');
	sys.stderr.write('Counting qnames present in sam_file2 that are missing from sam_file1...\n');
	for qname in sam_hash2.iterkeys():
		if (len(sam_hash2[qname]) > 0):
			if (sam_hash2[qname][0].IsMapped() == True):
				num_mapped_2 += 1;
			if (qname in sam_hash1.iterkeys()):
				pass;
			else:
				not_in_sam_file1 += 1;

	fp_out = None;
	fp_out_lt0bp = None;
	fp_out_gt5000bp = None;
	out_file = out_summary_prefix + '.csv';
	out_file_lt0bp = out_summary_prefix + '_lt0bp.csv';
	out_file_gt5000bp = out_summary_prefix + '_gt5000bp.csv';

	if (out_summary_prefix != ''):
		try:
			fp_out = open(out_file, 'w');
			fp_out_lt0bp = open(out_file_lt0bp, 'w');
			fp_out_gt5000bp = open(out_file_gt5000bp, 'w');
		except IOError:
			sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_file));
			return;
			# exit(1);

	summary_line = '';
	summary_line += 'sam_file1 = %s\n' % sam_file1;
	summary_line += 'sam_file2 = %s\n' % sam_file2;
	summary_line += 'not_in_sam_file1 = %d\n' % (not_in_sam_file1);
	summary_line += 'not_in_sam_file2 = %d\n' % (not_in_sam_file2);
	summary_line += 'num_different_reference = %d\n' % (num_different_reference);
	summary_line += 'num_different_orientation = %d\n' % (num_different_orientation);
	summary_line += 'num_not_mapped_1 = %d\n' % (num_not_mapped_1);
	summary_line += 'num_not_mapped_2 = %d\n' % (num_not_mapped_2);
	summary_line += 'num_mapped_1 = %d\n' % (num_mapped_1);
	summary_line += 'num_mapped_2 = %d\n' % (num_mapped_2);
	summary_line += '\n';

	length_threshold = 9000;

	sys.stderr.write(summary_line);
	if (out_summary_prefix != ''):
		fp_out.write(summary_line);
	summary_line = '';

	summary_line_lt0bp = '';
	summary_line_gt5000bp = '';

	num_same_alignments = 0;
	i = 0;
	# while i < len(distance_to_qname_hash.iterkeys()
	# print distance_to_qname_hash;
	for distance in sorted(distance_to_qname_hash.iterkeys()):
		sorted_by_length = sorted(distance_to_sam_hash[distance], reverse=True, key=lambda sam_line: len(sam_line.seq));
		# sorted_qnames = ['%s <%d, %d>' % (single_sam_line.qname, len(single_sam_line.seq), single_sam_line.mapq) for single_sam_line in sorted_by_length];
		sorted_qnames = ['%s <%d>' % (single_sam_line.qname, len(single_sam_line.seq)) for single_sam_line in sorted_by_length];

		sorted_qnames_above_length = [('%s' % (single_sam_line.qname)) for single_sam_line in sorted_by_length if (len(single_sam_line.seq) > length_threshold)];
		if (distance == 0):
			summary_line_lt0bp = ' \\\n'.join(sorted_qnames_above_length);
		if (distance > 5000):
			if (len(summary_line_gt5000bp) > 0):
				summary_line_gt5000bp += ' \\\n';
			summary_line_gt5000bp += ' \\\n'.join(sorted_qnames_above_length);

		# sorted_qnames = [str(len(single_sam_line.seq)) for single_sam_line in sorted(distance_to_sam_hash[distance], reverse=True, key=lambda sam_line: len(sam_line.seq))];
		# summary_line = str(distance) + '\t' + str(len(distance_to_qname_hash[distance])) + '\t' + '\t'.join(distance_to_qname_hash[distance]) + '\n';
		summary_line = str(distance) + '\t' + str(len(distance_to_qname_hash[distance])) + '\t' + '\t'.join(sorted_qnames) + '\n';
		if (distance < distance_threshold):
			num_same_alignments += len(distance_to_qname_hash[distance]);

		# sys.stdout.write(summary_line);
		if (out_summary_prefix != ''):
			fp_out.write(summary_line);
		summary_line = '';

	summary_line = 'distance_threshold = %d\n' % distance_threshold;
	summary_line += 'num_same_alignments = %d\n' % num_same_alignments;
	summary_line += '\n';
	sys.stderr.write(summary_line);
	if (out_summary_prefix != ''):
		fp_out.write(summary_line);
		fp_out_lt0bp.write(summary_line_lt0bp);
		fp_out_gt5000bp.write(summary_line_gt5000bp);
		summary_line = '';
		summary_line_lt0bp = '';
		summary_line_gt5000bp = '';

	if (out_summary_prefix != ''):
		fp_out.close();
		fp_out_lt0bp.close();
		fp_out_gt5000bp.close();


	return [sam_hash1, sam_hash2, distance_to_qname_hash];

	# sys.stderr.write('Starting to count the number of correctly mapped bases in the tested SAM file!\n');

# src/sam_compare_cigars.py /home/ivan/work/eclipse-workspace/golden-bundle/alignments_for_testing/reads-simulated/PacBio-100k/escherichia_coli/graphmap-params_SSW_r-test-drive.sam /home/ivan/work/eclipse-workspace/golden-bundle/reads-simulated/PacBio-100k/escherichia_coli/reads.sam temp/cigcompare.temp
if __name__ == "__main__":
	if (len(sys.argv) < 3 or len(sys.argv) > 5):
		sys.stderr.write('Compares alignment positions between two SAM files.');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <input_sam_file_1> <input_sam_file_2> [distance_threshold] [<out_file_prefix>]\n' % sys.argv[0]);
		sys.stderr.write('\n');
		sys.stderr.write('\tdistance_threshold - default value is 100\n');
		exit(1);

	sam_file1 = sys.argv[1];
	sam_file2 = sys.argv[2];
	distance_threshold = 100;
	out_summary_prefix = '';

	if (len(sys.argv) >= 4):
		distance_threshold = int(sys.argv[3]);

	if (len(sys.argv) >= 5):
		out_summary_prefix = sys.argv[4];

	compare_two_sams(sam_file1, sam_file2, distance_threshold, out_summary_prefix);

	# print 'Percent correctly mapped bases: %.2f' % CompareCigars(sam_file, sam_reference_file, out_summary_prefix);
