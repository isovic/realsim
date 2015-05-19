#! /usr/bin/python

import os
import re
import sys
from collections import defaultdict

CIGAR_OPERATIONS_ALL = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'];
CIGAR_OPERATIONS_BASIC = ['M', 'I', 'D', 'S', 'H'];
CIGAR_OPERATIONS_EXTENDED = ['M', 'I', 'D', 'S', 'H', '=', 'X'];

# Global table containing mutation counts
mutCntTable = {letter: defaultdict(int) for letter in ('A', 'T', 'G', 'C')}

class SAMLine:
	def __init__(self, line='', sam_basename=''):
		if line == '':
			self.Clear();
		else:
			self.ParseLine(line, sam_basename);

	def Clear(self):
		# (I) These parameters are parsed from the SAM file.
		self.qname = '';
		self.flag = 0;
		self.rname = '';
		self.pos = 0;
		self.mapq = 0;
		self.cigar = '';
		self.mrnm = '';
		self.mpos = 0;
		self.isize = 0;
		self.seq = '';
		self.qual = '';
		self.optional = {};
		self.original_line = '';

		# (II) These parameters are parsed from the parameters parsed from the SAM file.
		self.clipped_pos = 0;			# If the read has clipping operation at the beginning, then self.clipped_pos = self.pos - num_clipped_bases . This is to allow easy comparison to the reference.
		self.clip_op_front = '';
		self.clip_count_front= 0;
		self.clip_op_back = '';
		self.clip_count_back = 0;
		self.chosen_quality = 0;			# In SAM, quality can be specified either as mapping quality (mapq), or as alignment score (AS optional parameter). For instance, LAST outputs AS, and not mapq. If mapq == 255, then it will be checked if there is an AS entry in the given optional parameters. If AS is not present, then self.chose_quality will be equal to self.mapq, that is, 255.

		# (III) These parameters are modified during evaluation of mapping.
		self.evaluated = 0;
		self.is_correct_ref_and_orient = 0;
		self.is_duplicate = 0;
		self.is_best_of_duplicates = 0;
		self.actual_ref_reverse = 0;
		self.is_filtered_out = False;
		self.num_occurances_in_sam_file = 0;

		self.is_header_deformed = 0;
		self.actual_ref_header = '';
		self.trimmed_header = '';
		self.actual_ref_pos = 0;
		self.mapped_pos_with_shift = 0;
		self.min_distance = -1;

		self.sam_basename = '';

		# For sanity checking - this is True if the SAM line contains all fields. It is false i.e. if file writing was interrupted before line was finished.
		self.line_fields_ok = False;

	def VerboseFormatLine(self):
		#line = 'qname = %s\tflag = %d\trname = %s\tpos = %d\tmapq = %d\tis_header_deformed = %d\tis_correct_ref_and_orient = %d\tis_duplicate = %d\tis_best_of_duplicates = %d\tactual_ref_pos = %d\tmapped_pos_with_shift = %d\tmin_distance = %d\tsecondary = %s\tactual_ref_reverse = %d\treverse = %d\tpaired = %s\t' % (
			#self.qname, self.flag, self.rname, self.pos, self. mapq, self.is_header_deformed, int(self.is_correct_ref_and_orient), int(self.is_duplicate), int(self.is_best_of_duplicates), self.actual_ref_pos, self.mapped_pos_with_shift, self.min_distance, str(self.IsSecondary()), self.actual_ref_reverse, self.IsReverse(), self.IsPaired());
		line = 'qname = %s\tactual_ref_pos = %d\tmapped_pos_with_shift = %d\tmin_distance = %d\tmapq = %d\tactual_ref_reverse = %d\treverse = %d\tis_correct_ref_and_orient = %d\tis_header_deformed = %d\tflag = %d\trname = %s\tpos = %d\tis_duplicate = %d\tis_filtered_out = %d\tis_best_of_duplicates = %d\tsecondary = %s\tpaired = %s\t' % (
			self.qname, self.actual_ref_pos, self.mapped_pos_with_shift, self.min_distance, self.mapq, self.actual_ref_reverse, self.IsReverse(), int(self.is_correct_ref_and_orient), self.is_header_deformed, self.flag, self.rname, self.pos, int(self.is_duplicate), int(self.is_filtered_out), int(self.is_best_of_duplicates), str(self.IsSecondary()), self.IsPaired());

		return line;

	def ParseLine(self, line, sam_basename=''):
		split_line = line.split('\t');

		if len(split_line) < 11:
			print 'ERROR: Line does not contain all mandatory SAM fields!';
			print 'Line: "%s"' % line;
			self.Clear();
			self.line_fields_ok = False;
			return;

		self.original_line = line;

		self.qname = split_line[0];
		self.flag = int(split_line[1]);
		self.rname = split_line[2];
		self.pos = int(split_line[3]);
		self.mapq = int(split_line[4]);
		self.cigar = split_line[5];
		self.mrnm = split_line[6];
		self.mpos = int(split_line[7]);
		self.isize = int(split_line[8]);
		self.seq = split_line[9];
		self.qual = split_line[10];

		self.optional = {};
		i = 11;
		while i < len(split_line):
			split_optional = split_line[i].split(':');
			if (len(split_optional) < 2):
				i += 1;
				continue;
			self.optional[split_optional[0].strip()] = split_optional[-1].strip();		# Example of an optional parameter: AS:i:717
			i += 1;

		self.chosen_quality = self.mapq + 0;
		if (self.chosen_quality == 255):
			try:
				self.chosen_quality = int(self.optional['AS']);
			except:
				pass;

		self.evaluated = 0;
		self.is_correct_ref_and_orient = 0;
		self.is_duplicate = 0;
		self.is_best_of_duplicates = 0;
		self.actual_ref_reverse = 0;
		self.is_filtered_out = False;

		self.is_header_deformed = 0;
		self.actual_ref_header = '';
		self.trimmed_header = '';
		self.actual_ref_pos = 0;
		self.mapped_pos_with_shift = 0;		# This is assigned only after evaluation of the mapping position. It is equal to self.clipped_pos if self.qname is found in the reference SAM file, otherwise it is equal to 0.
		self.min_distance = -1;

		self.clipped_pos = self.pos;		# The position of the read subtracted by the amount of clipped bases at the begining of the read.
		self.clip_op_front = '';
		self.clip_count_front= 0;
		self.clip_op_back = '';
		self.clip_count_back = 0;
		self.num_occurances_in_sam_file = 1;

		self.sam_basename = sam_basename;

		if (len(self.cigar) > 0 and self.cigar != '*'):
			m_front = re.match("^([\d]+)([SH])", self.cigar);
			m_back = re.match(".*[MIDNSHP=X]([\d]+)([SH])$", self.cigar);

			if (m_front):
				cigarcount = int(m_front.group(1));
				cigarop = m_front.group(2);
				self.clipped_pos -= cigarcount;
				self.clip_op_front = cigarop;
				self.clip_count_front = cigarcount;
				#print '(front) cigarcount = %d, cigarop = %s, sam_basename = %s, qname = %s' % (cigarcount, cigarop, self.sam_basename, self.qname);
			if (m_back):
				cigarcount = int(m_back.group(1));
				cigarop = m_back.group(2);
				self.clip_op_back = cigarop;
				self.clip_count_back = cigarcount;
				#print '(back) cigarcount = %d, cigarop = %s, sam_basename = %s, qname = %s' % (cigarcount, cigarop, self.sam_basename, self.qname);

		#if (self.pos != self.clipped_pos):
			#print 'Clipped position: %d (original: %d, clip_op = %s, clip_count = %d)' % (self.clipped_pos, self.pos, self.clip_op, self.clip_count);
		self.line_fields_ok = True;

	def Verbose(self):
		print 'qname = %s' % self.qname;
		print 'flag = %s' % self.flag;
		print 'rname = %s' % self.rname;
		print 'pos = %s' % self.pos;
		print 'mapq = %s' % self.mapq;
		print 'cigar = %s' % self.cigar;
		print 'mrnm = %s' % self.mrnm;
		print 'mpos = %s' % self.mpos;
		print 'isize = %s' % self.isize;
		print 'seq = %s' % self.seq;
		print 'qual = %s' % self.qual;
		print '(evaluated = %d)' % self.evaluated;
		print '(min_distance = %d)' % self.min_distance;
		print '(is_correct_ref_and_orient = %d)' % self.is_correct_ref_and_orient;
		print '(is_duplicate = %d)' % self.is_duplicate;
		print '(is_filtered_out = %d)' % self.is_filtered_out;
		print '(is_best_of_duplicates = %d)' % self.is_best_of_duplicates;
		print '(clipped_pos = %d)' % self.clipped_pos;

	def FormatAccuracy(self):
		line = '';

		#query.min_distance = ret_min_distance;
		#query.actual_ref_pos = ret_reference_pos;
		#query.actual_ref_reverse = ret_ref_reverse;
		#query.actual_ref_header = ret_ref_header;
		#query.mapped_pos_with_shift = ret_mapped_pos;

		line += 'distance = %d\t' % self.min_distance;
		line += 'header_hit = %s\t' % (str(self.rname == self.actual_ref_header));
		line += 'reverse_hit = %s\t' % (str(self.actual_ref_reverse == self.IsReverse()));
		line += 'qname = %s\t' % self.qname;
		line += 'mapped = %s\t' % str(self.IsMapped());
		#line = 'header_def = %s\t'
		line += 'ref_header = %s\t' % self.actual_ref_header;
		line += 'map_header = %s\t' % self.rname;
		line += 'ref_pos = %d\t' % self.actual_ref_pos;
		line += 'map_pos_clipped = %d\t' % self.clipped_pos;
		line += 'ref_reverse = %s\t' % str(self.actual_ref_reverse);
		line += 'map_reverse = %s\t' % str(self.IsReverse());
		#line += 'duplicate = %s\t' % (str(self.is_duplicate != 0));
		line += 'cigar_start = %s\t' % (self.cigar[0:5]);
		line += 'cigar_end = %s\t' % (self.cigar[(len(self.cigar)-5):]);
		line += 'chosen_quality = %d\t' % (self.chosen_quality);
		line += 'num_occ = %d' % (self.num_occurances_in_sam_file);

		return line;



	# Checks the SAM flag to see if the read is paired end.
	def IsPaired(self):
		return (True if (self.flag & 0x01) else False);

	# Checks the SAM flag to see if the alignment is reported as mapped.
	def IsMapped(self):
		return (False if (self.flag & 0x04) else True);
		#return (True if (self.flag & 0x04) else False);

	# Checks the SAM flag to see if the SEQ was reversed.
	def IsReverse(self):
		return (True if (self.flag & 0x10) else False);

	# Checks the SAM flag to see if this is a secondary alignment.
	def IsSecondary(self):
		return (True if (self.flag & 0x100) else False);

	# Splits the CIGAR string into individual operations, in the
	# same format as the original CIGAR string is in.
	# The result is returned as an array of a 2-element array, e.g.
	# [[12, 'M'], [3, 'D'], [4, 'M']].
	def SplitCigar(self):
		i = 0;
		cigarcount_string = '';
		cigar_operations = [];
		if (self.IsMapped() == False):
			return cigar_operations;
		while i < len(self.cigar):
			if (self.cigar[i] in CIGAR_OPERATIONS_EXTENDED):
				cigar_operations.append([int(cigarcount_string), self.cigar[i]]);
				cigarcount_string = '';
			else:
				cigarcount_string += self.cigar[i];
			i += 1;
		if (cigarcount_string != ''):
			print 'ERROR: Faulty CIGAR string!';
			print 'cigarcount_string: ', cigarcount_string
			print 'i = ', i;
			print self.original_line

			cigar_operations = [];
		return cigar_operations;

	# Splits the CIGAR string into individual operations. Unlike SplitCigar,
	# this function also converts the extended cigar format to the basic cigar
	# format. This includes counting the number of successive M, = and X operations.
	# The result is returned as an array of a 2-element array, e.g.
	# [[12, 'M'], [3, 'D'], [4, 'M']].
	def SplitCigarInBasicFormat(self):
		i = 0;
		cigarcount_string = '';
		cigar_operations = [];
		if (self.IsMapped() == False):
			return cigar_operations;
		while i < len(self.cigar):
			if (self.cigar[i] in CIGAR_OPERATIONS_EXTENDED):
				# Check if it is a match/mismatch operation and convert to M:
				if (self.cigar[i] in 'M=X'):
					# If it is a match/mismatch and the list is empty, simply add it.
					if (len(cigar_operations) == 0):
						cigar_operations.append([int(cigarcount_string), 'M']);
					# If the list is not empty, check if the previous operation was an 'M'.
					elif (cigar_operations[-1][1] == 'M'):
						cigar_operations[-1][0] += int(cigarcount_string);
					# If the previous operation was not an M, then simply add it again.
					else:
						cigar_operations.append([int(cigarcount_string), 'M']);
				else:
					cigar_operations.append([int(cigarcount_string), self.cigar[i]]);
				cigarcount_string = '';
			else:
				cigarcount_string += self.cigar[i];
			i += 1;
		if (cigarcount_string != ''):
			print 'ERROR: Faulty CIGAR string!';
			print 'cigarcount_string: ', cigarcount_string
			print 'i = ', i;
			print self.original_line

			cigar_operations = [];
		return cigar_operations;


	# Calculates and returns extended cigar string
	# = - match, X - mismatch
	# Compares all M operations areas with the reference genome to find mismaches
	def CalcExtendedCIGAR(self, reference):
		extcigar = ''
		operations = self.SplitCigar()
		refpos = self.pos-1	# start of alignment within referece, clipped bases are skipped
							# -1 is beacuse SAM position is one based!
		seqpos = 0			# startng position within original sequence

		for op in operations:
			opsize = op[0]
			optype = op[1]
			if optype == 'M':
				mcigar = getExtendedCIGAR(reference[refpos:refpos+opsize], self.seq[seqpos:seqpos+opsize])
				extcigar += mcigar
				refpos += opsize
				seqpos += opsize
			elif optype in 'IS':
				seqpos += opsize
				extcigar += '%s%s' % (opsize, optype)
			elif optype == 'D':
				refpos += opsize
				extcigar += '%s%s' % (opsize, optype)
			elif optype == 'H':
				extcigar += '%s%s' % (opsize, optype)
			else:
				print 'ERROR: Faulty basic CIGAR string operation!'
				print 'operation: %s%s' % (opsize, optype)
				print 'self.cigar'

		return extcigar


	# Calculates a GC content of a given read by looking at corresponding bases in the reference
	def CalcGCContent(self, reference):
		GCcontent = 0.0
		GCcount = 0
		length = self.CalcReferenceLengthFromCigar()
		refpos = self.pos-1

		for base in reference[pos:pos+length]:
			if base in 'GC':
				GCcount += 1

		GCcontent = float(GCcount)/length

		return GCcontent


	# Calculates both extended CIGAR and GC content at the same time
	# To pass reference as an argument only once
	def CalcExtendedCIGARandGCContent(self, reference):
		GCcontent = 0.0
		GCcount = 0
		length = self.CalcReferenceLengthFromCigar()
		extcigar = ''
		operations = self.SplitCigar()
		refpos = self.pos-1	# start of alignment within referece, clipped bases are skipped
							# -1 is beacuse SAM position is one based!
		seqpos = 0			# startng position within original sequence

		for base in reference[refpos:refpos+length]:
			if base in 'GC':
				GCcount += 1

		# Calculating GC content
		GCcontent = float(GCcount)/length

		# Calculating extended CIGAR
		for op in operations:
			opsize = op[0]
			optype = op[1]
			if optype == 'M':
				mcigar = getExtendedCIGAR(reference[refpos:refpos+opsize], self.seq[seqpos:seqpos+opsize])
				extcigar += mcigar
				refpos += opsize
				seqpos += opsize
			elif optype in 'IS':
				seqpos += opsize
				extcigar += '%s%s' % (opsize, optype)
			elif optype == 'D':
				refpos += opsize
				extcigar += '%s%s' % (opsize, optype)
			elif optype == 'H':
				extcigar += '%s%s' % (opsize, optype)
			else:
				print 'ERROR: Faulty basic CIGAR string operation!'
				print 'operation: %s%s' % (opsize, optype)
				print 'self.cigar'

		return extcigar, GCcontent



	# Sums the counts of M/I/S/=/X operations in the CIGAR string to determine the
	# length of SEQ. This is used for testing, debugging and sanity checking. If the
	# SAM line has been properly formatted, the result of this function should be
	# the same as len(self.seq).
	def CalcAlignmentLengthFromCigar(self):
		split_cigar = self.SplitCigar();
		alignment_length = 0;
		i = 0;
		while i < len(split_cigar):
			cigar_count = split_cigar[i][0];
			cigar_op = split_cigar[i][1];
			# From the SAM format specification:
			#     Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
			if (cigar_op in 'MIS=X'):
				alignment_length += cigar_count;
			i += 1;
		return alignment_length;

	# Sums the counts of M/D/=/X operations in the CIGAR string to determine the
	# length of the reference covered by the read.
	def CalcReferenceLengthFromCigar(self):
		split_cigar = self.SplitCigar();
		alignment_length = 0;
		i = 0;
		while i < len(split_cigar):
			cigar_count = split_cigar[i][0];
			cigar_op = split_cigar[i][1];
			# From the SAM format specification:
			#     Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
			if (cigar_op in 'MD=X'):
				alignment_length += cigar_count;
			i += 1;
		return alignment_length;

	# Given a base position within the read, find its position on the reference by inspecting CIGAR.
	def FindBasePositionOnReference(self, base_position_in_read):
		cigar_pos_list = self.CalcCigarStartingPositions(True);
		i = 0;
		while (i < len(cigar_pos_list)):
			[cigar_count, cigar_op, pos_on_reference, pos_on_read] = cigar_pos_list[i];
			if (cigar_op in 'M=X'):
				if (pos_on_read == base_position_in_read):
					return pos_on_reference;
				elif (pos_on_read < base_position_in_read and (pos_on_read + cigar_count) > base_position_in_read):
					return (pos_on_reference + (base_position_in_read - pos_on_read));
			i += 1;
		return -1;

	# Splits the CIGAR string into individual operations, and determines their starting positions
	# on the reference. I.e. this function identifies the positions of matches, mismatches and deletions.
	# This is used, for example, in comparison of alignment operations to the reference SAM file.
	# If separate_matches_in_individual_bases == True, then 'M', '=' and 'X' operations will
	# be split into individual operations of length 1 (i.e. if the input CIGAR operation was
	# 3M, this would generate 1M1M1M. This is a nice hack for quick checking coordinates of
	# individual matching/mismatching bases.)
	# Function returns an array of arrays:
	# [[cigar_count1, cigar_op1, pos_on_reference1, pos_on_read1], [[cigar_count2, cigar_op2, pos_on_reference2, pos_on_read2], ...]
	# Example usage:
	# import utility_sam;
	# [headers, sam_lines] = utility_sam.LoadSAM(sam_path);
	# for sam_line in sam_lines:
	#	cigar_pos_list = sam_line.CalcCigarStartingPositions();
	#	print cigar_pos_list;

	def CalcCigarStartingPositions(self, separate_matches_in_individual_bases=False):
		#cigar_list = self.SplitCigar();
		cigar_list = self.SplitCigarInBasicFormat();
		cigar_pos_list = [];
		pos_on_reference = self.clipped_pos + 0;
		pos_on_read = 0;

		i = 0;
		while (i < len(cigar_list)):
			cigar_count = cigar_list[i][0];
			cigar_op = cigar_list[i][1];

			if (separate_matches_in_individual_bases == False):
				cigar_pos_list.append([cigar_count + 0, cigar_op + '', pos_on_reference + 0]);
				#if (cigar_op in 'MDSH=X'):
					#pos_on_reference += cigar_count;

				# S and H are also used to walk down the reference, because pos_on_reference was initialized
				# with self.clipped_pos, which is calculated from the SAM pos field by subtracting the number
				# of clipped bases. Otherwise, S and H should not be used to increase pos_on_reference.
				if (cigar_op in 'MSH=X'):
					pos_on_reference += cigar_count;
					pos_on_read += cigar_count;
				elif (cigar_op == 'D'):
					pos_on_reference += cigar_count;
				elif (cigar_op == 'I'):
					pos_on_read += cigar_count;

			else:
				if ((cigar_op in 'M=X') == False):
					cigar_pos_list.append([cigar_count + 0, cigar_op + '', pos_on_reference + 0, pos_on_read + 0]);
					if (cigar_op in 'SH'):
						pos_on_reference += cigar_count;
						pos_on_read += cigar_count;
					elif (cigar_op == 'D'):
						pos_on_reference += cigar_count;
					elif (cigar_op == 'I'):
						pos_on_read += cigar_count;
				else:
					j = 0;
					while (j < cigar_count):
						cigar_pos_list.append([1, cigar_op, pos_on_reference + 0, pos_on_read + 0]);
						pos_on_reference += 1;
						pos_on_read += 1;
						j += 1;

			i += 1;

		return cigar_pos_list;

	def CalcNumMappedBases(self):
		num_mapped_bases = len(self.seq);
		if (self.clip_op_front == 'S'):
			num_mapped_bases -= self.clip_count_front;
		if (self.clip_op_back == 'S'):
			num_mapped_bases -= self.clip_count_back;
		return num_mapped_bases;

	def IsAlignmentSane(self):
		if (self.IsMapped() == False):
			return True;

		if (self.line_fields_ok == False):
			return False;

		seq_len = len(self.seq);

		split_cigar = self.SplitCigar();
		alignment_length_on_seq = 0;
		i = 0;
		while i < len(split_cigar):
			cigar_count = split_cigar[i][0];
			cigar_op = split_cigar[i][1];
			# From the SAM format specification:
			#     Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
			if (cigar_op in 'MIS=X'):
				alignment_length_on_seq += cigar_count;

			if ((cigar_op == 'D' or cigar_op == 'I') and cigar_count > 100):
				return False;
			if ((i + 1) < len(split_cigar) and (cigar_op in 'DI') and (split_cigar[i + 1][1] in 'DISH')):
				return False;

			i += 1;

		if (seq_len != alignment_length_on_seq):
			return False;

		return True;




# A helper function.
# Takes two string of the same length and return corresponding
# extended CIGAR string looking only at matches and mismatches
# = - match, X - mismatch
# Used to break up longer M operations in basic CIGAR into multiple = and X operations
# If attribute countMutations is set to True, function also counts different mutations
# considering string1 as reference and string2 as a read
def getExtendedCIGAR(string1, string2, countMutations = True):
	if len(string1) != len(string2):
		raise Exception('Cannot callculate extended CIGAR for strings of different length!')

	match = '='
	mismatch = 'X'
	skipped = 'N'

	extcigar = ''
	# checking first characters to determine first operation
	if string1[0] == 'N' or string2[0] == 'N':
		optype = skipped
	elif string1[0] == string2[0]:
		optype = match
	else:
		optype = mismatch
		if countMutations:
			mutCntTable[string1[0]][string2[0]] += 1
	opsize = 1

	# checking other characeters
	for i in xrange(1, len(string1)):
		if string1[i] == 'N' or string2[i] == 'N':
			newoptype = skipped
		elif string1[i] == string2[i]:
			newoptype = match
		else:
			newoptype = mismatch
			if countMutations:
				mutCntTable[string1[i]][string2[i]] += 1

		if newoptype == optype:
			opsize += 1
		else :
			extcigar += '%d%s' % (opsize, optype)
			optype = newoptype
			opsize = 1

	# Last operations was not noted in CIGAR string
	extcigar += '%d%s' % (opsize, optype)

	return extcigar


# Functions for working with global mutation count table
def init_mutCntTable():
	mut_cnt_table = {letter: defaultdict(int) for letter in ('A', 'T', 'G', 'C')}


def LoadSAM(sam_path, verbose=False):
	try:
		fp_reference = open(sam_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % sam_path;
		return [[], []];

	headers = [];
	sam_lines = [];

	i = 0;
	for line in fp_reference:
		i += 1;
		if (verbose == True):
			sys.stdout.write('\rParsing SAM line %d...' % (i));
		line = line.strip();
		if len(line) == 0:
			continue;
		if line[0] == '@':
			headers.append(line);
			continue;
		sam_line = SAMLine(line);
		sam_lines.append(sam_line);
	fp_reference.close();

	if (verbose == True):
		sys.stdout.write('done!\n');

	return [headers, sam_lines];

def LoadOnlySAMHeaders(sam_path, verbose=False):
	try:
		fp_reference = open(sam_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % sam_path;
		return [];

	headers = [];

	i = 0;
	for line in fp_reference:
		i += 1;
		if (verbose == True):
			sys.stdout.write('\rParsing SAM line %d...' % (i));
		line = line.strip();
		if len(line) == 0:
			continue;
		if line[0] == '@':
			headers.append(line);
			continue;
		else:
			break;
	fp_reference.close();

	if (verbose == True):
		sys.stdout.write('done!\n');

	return headers;

# Hashes the entries from a SAM file by their QNAME for faster lookup during comparison.
def HashSAM(sam_path):
	try:
		fp_reference = open(sam_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % sam_path;
		return [{}, 0];

	ret = {};

	num_unique_references = 0;
	num_references = 0;

	for line in fp_reference:
		line = line.strip();

		if len(line) == 0:
			continue;

		if line[0] == '@':
			continue;

		sam_line = SAMLine(line);

		modified_qname = sam_line.qname;
		#modified_qname = '/'.join(sam_line.qname.split('/')[:-1]);

		try:
			current_hash = ret[modified_qname];

			should_be_counted = True;

			# We have to check if this sequence is actually a paired read of another sequence, or there is something funny in the data.
			for existing_sam in current_hash:
				if (sam_line.IsSecondary() == True) or (sam_line.IsSecondary() == existing_sam.IsSecondary() and sam_line.IsReverse() == existing_sam.IsReverse()):
					# This case should not be counted. It means that, either the alignment is marked as secondary which means there should be a primary alignment as well, or that there is more than one primary alignment, and the orientation is the same, which means that an aligner is trying to artificially boost-up their statistics.
					should_be_counted = False;
					break;


			#if sam_line.qname == 'gi|48994873|gb|U00096.2|-463960':
				#print line;
				#for sam in ret[sam_line.qname]:
					#print 'is_secondary: %s, is_reverse = %s, num_unique_references: %d' % (str(sam.IsSecondary()), str(sam.IsReverse()), num_unique_references);

			if should_be_counted == True:
				num_unique_references += 1;	# Count only unique sequences.
				#print 'Tu sam 2!';
			#if sam_line.qname == 'gi|48994873|gb|U00096.2|-463960':
				#print '---';

			# At least one sequence with the same name has already been added.
			current_hash.append(sam_line);
		except:
			# This is a new sequence, unhashed before. Create a new list and count the sequence.
			ret[modified_qname] = [sam_line];
			if sam_line.IsSecondary() == False:
				num_unique_references += 1;	# Count only unique sequences, but not secondary ones.


			#if sam_line.qname == 'gi|48994873|gb|U00096.2|-463960':
				#print line;
				#if sam_line.IsSecondary() == False:
					#print 'Tu sam 1!';

				#print '---';

		#if sam_line.qname == 'gi|48994873|gb|U00096.2|-463960':
			#for sam in ret[sam_line.qname]:
				#print 'is_secondary: %s, is_reverse = %s, num_unique_references: %d' % (str(sam.IsSecondary()), str(sam.IsReverse()), num_unique_references);
			#print '---';

		num_references += 1;	# Count only unique sequences.

	fp_reference.close();

	for key in ret.keys():
		ret[key].sort(reverse=True, key=lambda sam_line: sam_line.chosen_quality);

	return [ret, num_references, num_unique_references];

# Hashes the entries from a SAM file by their QNAME for faster lookup during comparison.
def HashSAMWithFilter(sam_path, qname_hash_to_filter):
	try:
		fp_reference = open(sam_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % sam_path;
		return [{}, 0, 0];

	ret = {};

	num_unique_references = 0;
	num_references = 0;

	for line in fp_reference:
		line = line.strip();

		if len(line) == 0:
			continue;

		if line[0] == '@':
			continue;

		sam_line = SAMLine(line);

		modified_qname = sam_line.qname;
		#modified_qname = '/'.join(sam_line.qname.split('/')[:-1]);

		try:
			if (qname_hash_to_filter[modified_qname] == 1):
				sam_line.is_filtered_out = True;
		except:
				sam_line.is_filtered_out = False;

		try:
			current_hash = ret[modified_qname];

			should_be_counted = True;

			# We have to check if this sequence is actually a paired read of another sequence, or there is something funny in the data.
			for existing_sam in current_hash:
				if (sam_line.IsSecondary() == True) or (sam_line.IsSecondary() == existing_sam.IsSecondary() and sam_line.IsReverse() == existing_sam.IsReverse()):
					# This case should not be counted. It means that, either the alignment is marked as secondary which means there should be a primary alignment as well, or that there is more than one primary alignment, and the orientation is the same, which means that an aligner is trying to artificially boost-up their statistics.
					should_be_counted = False;
					break;

			if should_be_counted == True and sam_line.is_filtered_out == False:
				num_unique_references += 1;	# Count only unique sequences.

			# At least one sequence with the same name has already been added.
			current_hash.append(sam_line);
		except:
			# This is a new sequence, unhashed before. Create a new list and count the sequence.
			ret[modified_qname] = [sam_line];
			if sam_line.IsSecondary() == False and sam_line.is_filtered_out == False:
				num_unique_references += 1;	# Count only unique sequences, but not secondary ones.
		num_references += 1;	# Count only unique sequences.

	fp_reference.close();

	for key in ret.keys():
		ret[key].sort(reverse=True, key=lambda sam_line: sam_line.chosen_quality);

	return [ret, num_references, num_unique_references];

# This is deprecated, should be removed.
def CheckSamModified(sam_filename, out_path_prefix):
	sam_timestamp = str(os.path.getmtime(sam_filename));

	path_roc = out_path_prefix + '.roc';
	path_sum = out_path_prefix + '.sum';

	lines = [];

	try:
		fp_sum = open(path_sum, 'r');
		lines = fp_sum.readlines();
		fp_sum.close();

		fp_sum = open(path_roc, 'r');
		fp_sum.close();
	except IOError:
		return [True, sam_timestamp];

#	fp_roc.write('SAM timestamp: %s\n' % sam_timestamp);
	for line in lines:
		split_line = line.split(':');
		if split_line[0].strip() == 'SAM timestamp':
			if len(split_line) < 2:
				continue;
			last_sam_timestamp = split_line[1].strip();

			if last_sam_timestamp == sam_timestamp:
				return [False, sam_timestamp];

			break;

	return [True, sam_timestamp];

# Illumina and Roche reads that were generated using ART have different paired end notations in their headers. That
# makes sense if that simulates the original headers as would be produced by real technology.
# The problem is that, for some reason, ART trims the paired end notation (i.e. /1 and /2 for Illumina) from the qname
# field in the ground-truth SAM file. For this reason, function TrimHeader is used to remove the paired end notation from
# the Illumina and Roche headers.
def TrimHeader(header):
	ret = header;

	illumina_re = r'^(.*_[\d]+)-[12]';
	roche_re = r'^(.*-[\d]+)/[12]';
	match_header = None;
	if match_header == None:
		match_header = re.match(illumina_re, header);
	if match_header == None:
		match_header = re.match(roche_re, header);

	if match_header:
		ret = str(match_header.group(1));

	return ret;

def GetExecutionTime(sam_file):
	time_file = sam_file[0:-3] + 'time';

	try:
		fp_time = open(time_file, 'r');
		execution_time = fp_time.readline();
		fp_time.close();

		return execution_time;
	except IOError:
		return 'Execution time measurements not found!'

def GetExecutionStats(sam_file):
	time_file = sam_file[0:-3] + 'memtime';

	try:
		fp_time = open(time_file, 'r');
		lines = fp_time.readlines();
		fp_time.close();

	except IOError:
		return 'Execution time measurements not found!'

	#for line in lines:
		#line = line.strip();
		#line

	return ('\t' + '\n\t'.join([line.strip() for line in lines]));

def WriteSamLines(sam_lines, output_path):
	try:
		fp = open(output_path, 'w');
		i = 0;
		while i < len(sam_lines):
			sam_line = sam_lines[i];
#			temp_line = sam_line.VerboseFormatLine();
			temp_line = sam_line.FormatAccuracy();
			fp.write(temp_line + '\n');
			i += 1;
		fp.close();
	except IOError:
		print 'ERROR: Could not open file "%s" for writing!' % (output_path);
		return;

def GetBasicStats(sam_lines, allowed_distance=1):
	true_positive = 0;
	false_positive = 0;
	not_mapped = 0;

	for sam_line in sam_lines:
		if (sam_line.IsMapped() and sam_line.is_duplicate == 0):
			if (sam_line.is_correct_ref_and_orient == 1 and sam_line.min_distance <= allowed_distance):
				true_positive += 1;
			else:
				# print 'sam_line.min_distance > allowed_distance: %d > %d' % (sam_line.min_distance, allowed_distance);
				false_positive += 1;
		else:
			not_mapped = 0;

	return [true_positive, false_positive, not_mapped];

def GetDistanceHistogramStats(sam_lines, distance_limits):
	sorted_lines_by_distance = sorted(sam_lines, key=lambda sam_line: sam_line.min_distance);

	sorted_distance_limits = sorted(distance_limits);
	ret_distance_counts = [0 for distance_limit in sorted_distance_limits];

	current_distance_index = 0;

	total_correct = 0;

	i = 0;
	while i < len(sorted_lines_by_distance):
		sam_line = sorted_lines_by_distance[i];

		if (sam_line.IsMapped() and sam_line.is_duplicate == 0):
			if (sam_line.is_correct_ref_and_orient == 1):
				if (sam_line.min_distance < sorted_distance_limits[current_distance_index]):
					#ret_distance_counts[current_distance_index] += 1;
					total_correct += 1;
				else:
					ret_distance_counts[current_distance_index] = total_correct;
					current_distance_index += 1;
					if (current_distance_index >= len(sorted_distance_limits)):
						break;
					#ret_distance_counts[current_distance_index] = 0 + ret_distance_counts[current_distance_index - 1];
					continue;	# Do not increase i, because we want to check the current SAM line also.

		i += 1;

	if (current_distance_index > 0):
		ret_distance_counts[current_distance_index - 1] = total_correct;

	# If all the sam lines have been counted and the maximum observed distance is less than the maximum value of distance_limits parameter,
	# fill the rest of the array with the same value (the last value that was counted), to make the graph flatline.
	i = current_distance_index;
	while (i < len(sorted_distance_limits)):
		ret_distance_counts[i] = total_correct; # ret_distance_counts[current_distance_index];
		i += 1;

	# Nothing important, just format the X axis values for the counts.
	#min_sorted_distance_limit = sorted_distance_limits[0];
	#sorted_distance_limits_shifted = [(distance_limit - min_sorted_distance_limit) for distance_limit in sorted_distance_limits];
	#return [sorted_distance_limits_shifted, ret_distance_counts];

	return [sorted_distance_limits, ret_distance_counts];

def GetDistanceHistogramStatsScaleDuplicates(sam_lines, distance_limits, scale_by_num_occurances=True):
	sorted_lines_by_distance = sorted(sam_lines, key=lambda sam_line: sam_line.min_distance);

	sorted_distance_limits = sorted(distance_limits);
	ret_distance_counts = [0 for distance_limit in sorted_distance_limits];

	current_distance_index = 0;

	total_correct = 0.0;

	i = 0;
	while i < len(sorted_lines_by_distance):
		sam_line = sorted_lines_by_distance[i];

		if (sam_line.IsMapped()):
			#if (sam_line.is_correct_ref_and_orient == 1):
			if (sam_line.is_correct_ref_and_orient == 1 and sam_line.min_distance <= sorted_distance_limits[current_distance_index]):
				#ret_distance_counts[current_distance_index] += 1;
				if (scale_by_num_occurances == True):
					total_correct += (1.0 / float(sam_line.num_occurances_in_sam_file));
				else:
					#print 'total_correct = ', total_correct;
					total_correct += 1.0;

			else:
				if (sam_line.min_distance > sorted_distance_limits[current_distance_index]):
					ret_distance_counts[current_distance_index] = total_correct;
					current_distance_index += 1;
					if (current_distance_index >= len(sorted_distance_limits)):
						break;
					#ret_distance_counts[current_distance_index] = 0 + ret_distance_counts[current_distance_index - 1];
					continue;	# Do not increase i, because we want to check the current SAM line also.

		i += 1;

	if (current_distance_index > 0):
		ret_distance_counts[current_distance_index - 1] = total_correct;

	# If all the sam lines have been counted and the maximum observed distance is less than the maximum value of distance_limits parameter,
	# fill the rest of the array with the same value (the last value that was counted), to make the graph flatline.
	i = current_distance_index;
	while (i < len(sorted_distance_limits)):
		ret_distance_counts[i] = total_correct; # ret_distance_counts[current_distance_index];
		i += 1;

	# Nothing important, just format the X axis values for the counts.
	#min_sorted_distance_limit = sorted_distance_limits[0];
	#sorted_distance_limits_shifted = [(distance_limit - min_sorted_distance_limit) for distance_limit in sorted_distance_limits];
	#return [sorted_distance_limits_shifted, ret_distance_counts];

	return [sorted_distance_limits, ret_distance_counts];

#def GetDistanceHistogramStats(sam_lines, distance_limits):
	#sorted_lines_by_distance = sorted(sam_lines, key=lambda sam_line: sam_line.min_distance);

	#sorted_distance_limits = sorted(distance_limits);
	#ret_distance_counts = [0 for distance_limit in sorted_distance_limits];

	#current_distance_index = 0;

	##for sam_line in sorted_lines_by_distance:
	#i = 0;
	#while i < len(sorted_lines_by_distance):
		#sam_line = sorted_lines_by_distance[i];

		#if (sam_line.IsMapped()):
			#if (sam_line.is_correct_ref_and_orient == 1):
				#if (sam_line.min_distance < sorted_distance_limits[current_distance_index]):
					#ret_distance_counts[current_distance_index] += 1;
				#else:
					#current_distance_index += 1;
					#if (current_distance_index >= len(sorted_distance_limits)):
						#break;
					#ret_distance_counts[current_distance_index] = 0 + ret_distance_counts[current_distance_index - 1];
					#continue;	# Do not increase i, because we want to check the current SAM line also.

		#i += 1;

	#i = current_distance_index;
	#while (i < len(sorted_distance_limits)):
		#ret_distance_counts[i] = ret_distance_counts[current_distance_index];
		#i += 1;

	#min_sorted_distance_limit = sorted_distance_limits[0];
	#sorted_distance_limits_shifted = [(distance_limit - min_sorted_distance_limit) for distance_limit in sorted_distance_limits];

	#return [sorted_distance_limits_shifted, ret_distance_counts];

def GetMapqHistogramStats(sam_lines, mapq_limits, allowed_distance):
	sorted_lines_by_mapq = sorted(sam_lines, reverse=True, key=lambda sam_line: sam_line.mapq);
	sorted_mapq_limits = sorted(mapq_limits, reverse=True);
	ret_mapq_counts = [0 for mapq_limit in sorted_mapq_limits];

	current_mapq = 0;
	total_correct = 0;
	#for sam_line in sorted_lines_by_mapq:
	i = 0;
	while i < len(sorted_lines_by_mapq):
		sam_line = sorted_lines_by_mapq[i];

		if (sam_line.IsMapped() and sam_line.is_duplicate == 0):
			if (sam_line.is_correct_ref_and_orient == 1):
				if (sam_line.mapq >= sorted_mapq_limits[current_mapq]):
					if (sam_line.min_distance < allowed_distance):
						#ret_mapq_counts[current_mapq] += 1;
						total_correct += 1;
				else:
					#print '[%d, %d] %d, total_correct = %d' % (current_mapq, sorted_mapq_limits[current_mapq], sam_line.mapq, total_correct);
					ret_mapq_counts[current_mapq] = total_correct;
					current_mapq += 1;
					if (current_mapq >= len(sorted_mapq_limits)):
						break;

					continue;	# Do not increase i, because we want to check the current SAM line also.
					#total_correct += 1;
					#ret_mapq_counts[current_mapq] = 0 + ret_mapq_counts[current_mapq - 1];
		i += 1;

	ret_mapq_counts[current_mapq] = total_correct;

	i = current_mapq;
	while (i < len(sorted_mapq_limits)):
		ret_mapq_counts[i] = ret_mapq_counts[current_mapq];
		i += 1;

	return [allowed_distance, sorted_mapq_limits, ret_mapq_counts];

def FindMultipleQnameEntries(sam_files):
	#duplicate_list = [];
	duplicate_hash = {};

	for sam_file in sam_files:
		try:
			fp_sam = open(sam_file, 'r');
		except IOError:
			print 'ERROR: Could not open file "%s" for reading!' % sam_file;
			return [{}, 0];

		occurance_hash = {};
		for line in fp_sam:
			line = line.strip();
			if len(line) == 0:
				continue;
			if line[0] == '@':
				continue;
			sam_line = SAMLine(line);

			try:
				occurance_hash[sam_line.qname] += 1;
				duplicate_hash[sam_line.qname] = 1;
			except:
				occurance_hash[sam_line.qname] = 1;

		fp_sam.close();

	fp = open('test.multiple', 'w');
	fp.write('\n'.join(duplicate_hash.keys()));
	fp.close();

	return duplicate_hash;



def CountMappedReads(sam_file):
	fp = None;

	try:
		fp = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!\n' % ('CountMappedReads', sam_file));
		exit(1);

	num_alignments = 0;
	num_mapped_alignments = 0;
	sam_reads = {};
	sam_mapped_reads = {};
	num_mapped_bases = 0;	# Count of the number of non-clipped bases in an alignment. Only highest-scoring alignments are considered for each read,
	total_num_bases = 0;
	highest_scoring_alignment_for_read = {};	# Hash to keep track of the highest scoring alignment for a read.

	for line in fp:
		if (len(line) == 0 or line[0] == '@'):
			continue;

		sam_line = SAMLine(line.rstrip());

		# Count the number of occurances of each read in the SAM file (i.e. multiple mappings).
		try:
			sam_reads[sam_line.qname] += 1;
			#print sam_line.qname;
		except:
			sam_reads[sam_line.qname] = 1;

		if (sam_line.IsMapped() == True):
			# Count the occurances of mapped reads.
			try:
				sam_mapped_reads[sam_line.qname] += 1;
			except:
				sam_mapped_reads[sam_line.qname] = 1;

			try:
				if (highest_scoring_alignment_for_read[sam_line.qname].chosen_quality < sam_line.chosen_quality):
					highest_scoring_alignment_for_read[sam_line.qname] = sam_line;
			except:
				highest_scoring_alignment_for_read[sam_line.qname] = sam_line;

			num_mapped_alignments += 1;

		num_alignments += 1;

	fp.close();

	num_unique_reads = 0;
	for value in sam_reads.values():
		if (value == 1):
			num_unique_reads += 1;
	num_unique_mapped_reads = 0;
	for value in sam_mapped_reads.values():
		if (value == 1):
			num_unique_mapped_reads += 1;

	for best_read in highest_scoring_alignment_for_read.values():
		total_num_bases += len(sam_line.seq);
		num_mapped_bases += best_read.CalcNumMappedBases();

	num_mapped_reads = len(highest_scoring_alignment_for_read.keys());

	return [num_alignments, num_mapped_alignments, num_unique_reads, num_mapped_reads, num_mapped_bases];



if __name__ == "__main__":
	pass;
