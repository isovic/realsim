#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
import sys

from collections import defaultdict


# To generate random reads
import random

# For manipulating CIGAR strings
import re

# For working wit gzip files
# Compressing profile to reduce size and make it small enough for GitHub
import gzip

# Mainly for callculating a reverse complement of a sequence
# and other utility stuff
import fastqparser

# To be able to generate random reads in SAM format
import utility_sam


CIGAR_OPERATIONS_ALL = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'];
CIGAR_OPERATIONS_BASIC = ['M', 'I', 'D', 'S', 'H'];
CIGAR_OPERATIONS_EXTENDED = ['M', 'I', 'D', 'S', 'H', '=', 'X'];

# CIGAR profile version number
# Used to check if loaded profile is compatible with the current code
# Should be changed only when new code makes old profile unusable
PROFILE_VERSION = '1.1'


# CIGARProfile: Class containing information about a dataset
# The infromation is stored as read CIGAR strings, plus some additional stuff
# It is used to generate simulated reads applying stored CIGAR strings to a random part of a given reference
# TODO:
# - Decide how to generate quals
# - Decide how to generate Ns (lots of them in nanopore reads)
# - (least important) Decide what to do about short tandem repeats (STR), inserts and deletes both
# - decide what to do with clipping
#       - randomly generate bases
#       - take part of the reference

# CIGARLine: Class containg information about one CIGAR string
class CIGARLine:
    def __init__(self, cigar = '', position = -1, qname = '', rname = '', quality = -1, GCcontent = 0.0, quals = None):
        self.cigar = cigar
        self.position = position
        self.qname = qname
        self.rname = rname
        self.quality = quality
        self.GCcontent = GCcontent          # GC content of a read that generated a CIGAR string
        self.quals = quals                  # Quality values of a read that generated a CIGAR string

    def toTSVstring(self):
        tsvstring = ('%s\t%d\t%s\t%s\t%d\t%0.4f' % (self.cigar, self.position, self.qname, self.rname, self.quality, self.GCcontent))

        # Add quals if they are defined
        if self.quals is not None:
            tsvstring += '\t%s' % self.quals

        return tsvstring

    def fromTSVstring(self, line):
        # If present, removing \n from the end of line
        if line[-1] == '\n':
            line = line[:-1]
        elements = line.split('\t')
        self.cigar = elements[0]
        self.position = int(elements[1])
        self.qname = elements[2]
        self.rname  = elements[3]
        self.quality = int(elements[4])
        self.GCcontent = float(elements[5])

        # getting quals string, if it exists
        if len(elements) > 6:
            self.quals = elements[6]

    # Splits the CIGAR string into individual operations, in the
    # same format as the original CIGAR string is in.
    # The result is returned as an array of a 2-element tuples, e.g.
    # [(12, 'M'), (3, 'D'), (4, 'M')].
    def splitCigar(self):
        # Using regular expressions to find repeating digit and skipping one character after that
        pattern = '(\d+)(.)'
        operations = re.findall(pattern, self.cigar)

        return operations

    # Returns an aproximate length of the read that generated the cigar string
    def length(self):
        return sum([int(x[0]) for x in self.splitCigar()])


class CIGARProfile:
    def __init__(self, name, GCcontent = 0.0):
        self.name = name
        self.clines = []

        # GC content for the whole profile
        # For GC content to be relevan't, the profile should be constructed for a single genome
        self.GCcontent = GCcontent

        # Table containing count of different mutations
        # Used to generate different mutations with different probability
        self.mutCntTable = None

        # Placeholder that allows the definition of profile-wide distribution for quals and Ns
        self.qualDistrib = None
        self.NDistrib = None

        # Profile version, used to determine if loaded profile is compatible with the current code
        self.version = PROFILE_VERSION

    def appendCLine(self, cline):
        if cline is None:
            raise Exception('Trying to add nonexisting CIGAR line! (%d)' % len(self.clines))
        self.clines.append(cline)

    def getRandomCLine(self):
        rand = random.randint(0, len(self.clines)-1)
        return self.clines[rand]

    # ATM using all cigar lines to generate reads at random positions
    # TODO: Consider the situation where reference is too short for a randomly chosen CIGAR
    def generateRandomRead(self, reference):
        cline = self.getRandomCLine()
        reflen = len(reference)
        readlen = cline.length()

        # Generate random position
        randpos = random.randint(0, reflen-readlen)
        startread = reference[randpos:randpos+readlen]

        # 50% of the time generate reverse complement sequence
        complement = random.choice((True, False))
        if complement:
            startread = fastqparser.revcomp_seq(startread)

        read = ''
        pos = 0         # Where in the starting read are we currently positioned

        operations = cline.splitCigar()
        # applying cigar operations
        for op in operations:
            opsize = int(op[0])
            optype = op[1]
            if optype in ('H', 'S'):
                # Clipping: Generate random bases, atm just putting Hs
                bases = generate_randomseq(opsize)
                read += bases
                # read += 'H'*opsize
            elif optype == '=':
                # Copy from the starting read / reference
                read += startread[pos:pos+opsize]
                pos += opsize
            elif optype == 'X':
                for i in xrange(opsize):
                    read += mutate(startread[pos+i], self.mutCntTable)
                # read += 'X'*opsize
                pos += opsize
            elif optype == 'I':
                # Insert random bases, atm just putting Is
                # Could be joined together with clipping
                bases = generate_randomseq(opsize)
                read += bases
                # read += 'I'*opsize
            elif optype == 'D':
                # Skip bases in starting read
                pos += opsize
            elif optype == 'N':
                read += 'N'*opsize
                pos += opsize
            else:
                sys.stderr.write('\nInvalid CIGAR operation for generating random reads: %s\n' % optype)
                exit(1)

        # Return generated read in a form of a SAMline
        sline = utility_sam.SAMLine()
        sline.qname = ''
        sline.flag = 16 if complement else 0        # Setting forward/reverse flag
                                                    # Ignoring other possibilities, no secondary alignments and paired reads
        sline.rname = ''
        sline.pos = randpos
        sline.mapq = 255
        sline.cigar = cline.cigar
        sline.mrnm = ''
        sline.mpos = 0
        sline.isize = 0
        sline.seq = read
        sline.qual = cline.quals
        sline.original_line = ''

        return sline


    # Deprecated:
    # Separate function to generate a whole dataset faster and with less memory
    # (instead of calling generateRandomRead many times)
    def generateRandomReadsByNumber(self, reference, numreads):
        reads = []
        reflen = len(reference)

        # Generating reads
        for i in xrange(numreads):
            cline = self.getRandomCLine()
            readlen = cline.length()
            # Generate random position
            randpos = random.randint(0, reflen-readlen)
            startread = reference[randpos:randpos+readlen]

            # 50% of the time generate reverse complement sequence
            complement = random.choice((True, False))
            if complement:
                startread = fastqparser.revcomp_seq(startread)

            read = ''
            pos = 0         # Where in the starting read are we currently positioned

            operations = cline.splitCigar()
            # applying cigar operations
            for op in operations:
                opsize = int(op[0])
                optype = op[1]
                if optype in ('H', 'S'):
                    # Clipping: Generate random bases, atm just putting Hs
                    bases = generate_randomseq(opsize)
                    read += bases
                    # read += 'H'*opsize
                elif optype == '=':
                    # Copy from the starting read / reference
                    read += startread[pos:pos+opsize]
                    pos += opsize
                elif optype == 'X':
                    for i in xrange(opsize):
                        read += mutate(startread[pos+i],self.mutCntTable)
                    # read += 'X'*opsize
                    pos += opsize
                elif optype == 'I':
                    # Insert random bases, atm just putting Is
                    # Could be joined together with clipping
                    bases = generate_randomseq(opsize)
                    read += bases
                    # read += 'I'*opsize
                elif optype == 'D':
                    # Skip bases in starting read
                    pos += opsize
                elif optype == 'N':
                    read += 'N'*opsize
                    pos += opsize
                else:
                    sys.stderr.write('\nInvalid CIGAR operation for generating random reads: %s\n' % optype)
                    exit(1)

            reads.append((read, cline.quals))

        return reads


    # Deprecated:
    # Separate function to generate a whole dataset faster and with less memory
    # (instead of calling generateRandomRead many times)
    def generateRandomReadsByCoverage(self, reference, coverage):
        reads = []
        reflen = len(reference)

        numbases = 0
        numreads = 0

        while True:
            cline = self.getRandomCLine()
            readlen = cline.length()
            # Generate random position
            randpos = random.randint(0, reflen-readlen)
            startread = reference[randpos:randpos+readlen]

            # 50% of the time generate reverse complement sequence
            complement = random.choice((True, False))
            if complement:
                startread = fastqparser.revcomp_seq(startread)

            read = ''
            pos = 0         # Where in the starting read are we currently positioned

            operations = cline.splitCigar()
            # applying cigar operations
            for op in operations:
                opsize = int(op[0])
                optype = op[1]
                if optype in ('H', 'S'):
                    # Clipping: Generate random bases, atm just putting Hs
                    bases = generate_randomseq(opsize)
                    read += bases
                    # read += 'H'*opsize
                elif optype == '=':
                    # Copy from the starting read / reference
                    read += startread[pos:pos+opsize]
                    pos += opsize
                elif optype == 'X':
                    for i in xrange(opsize):
                        read += mutate(startread[pos+i], self.mutCntTable)
                    # read += 'X'*opsize
                    pos += opsize
                elif optype == 'I':
                    # Insert random bases, atm just putting Is
                    # Could be joined together with clipping
                    bases = generate_randomseq(opsize)
                    read += bases
                    # read += 'I'*opsize
                elif optype == 'D':
                    # Skip bases in starting read
                    pos += opsize
                elif optype == 'N':
                    read += 'N'*opsize
                    pos += opsize
                else:
                    sys.stderr.write('\nInvalid CIGAR operation for generating random reads: %s\n' % optype)
                    exit(1)

            reads.append((read, cline.quals))
            numbases += len(read)
            numreads += 1

            # Finish if a set coverage is reached
            if numbases/reflen > coverage:
                break

        return reads


def loadCProfile(filepath):

    cprofile = None
    try:
        if filepath.endswith('.cpf'):
            pfile = open(filepath, 'r')
        elif filepath.endswith('.cpf.gz'):
            pfile = gzip.GzipFile(filepath, 'r')
        else:
            sys.stderr.write('\n\nInvalide CIGAR profile file extension! Exiting....\n')
            exit(1)

        # loading profile name from the first line
        line = pfile.readline()
        # Removing '#CProfile: ' from the start and \n from the end of line
        name = line[11:-1]

        # Loading CProfile version number
        line = pfile.readline()
        # Removing '#CProfile version: ' from the start and \n from the end of line
        version = line[19:-1]
        if version != PROFILE_VERSION:
            raise Exception('Invalid profile version in %s(%s)! Expected %s' % (filepath, version, PROFILE_VERSION))

        # loading GC content
        line = pfile.readline()
        # Removing '#GC Content: ' from the start and \n from the end of line
        GCcontent = float(line[13:-1])
        cprofile = CIGARProfile(name, GCcontent)
        cprofile.version = version
        cprofile.mutCntTable = {}
        # Loading mutation count table
        line = pfile.readline()
        for base in ['A', 'C', 'G', 'T']:
            line = pfile.readline()
            elements = line.split('\t')
            if len(elements) != 7:
                sys.stderr.write('Invalid Mutation table row\n')
                return None
            mutCntDict = defaultdict(int)
            mutCntDict[elements[1]] = int(elements[2])
            mutCntDict[elements[3]] = int(elements[4])
            mutCntDict[elements[5]] = int(elements[6])
            cprofile.mutCntTable[elements[0]] = mutCntDict

        # Skipping line with headings
        line = pfile.readline()
        for line in pfile:
            cline = CIGARLine()
            cline.fromTSVstring(line)
            cprofile.appendCLine(cline)

        pfile.close()

    # TODO:
    # Load information about quals and Ns distribution
    except Exception:
        sys.stderr.write('\n\nError loading a profile!\n')
        exit(1)

    return cprofile


def storeCProfile(filepath, cprofile):

    pfile = None
    try:
        if filepath.endswith('.cpf'):
            pfile = open(filepath, 'w')
        elif filepath.endswith('.cpf.gz'):
            pfile = gzip.GzipFile(filepath, 'w')
        else:
            sys.stderr.write('\n\nInvalide CIGAR profile file extension! Exiting....\n')
            exit(1)

        # Writing profile name in the first line with # in front
        pfile.write('#CProfile: %s\n' % cprofile.name)
        # Writing profile version
        pfile.write('#CProfile version: %s\n' % cprofile.version)
        # Writing GCContent
        pfile.write('#GC Content: %0.4f\n' % cprofile.GCcontent)
        # Writing mutation count table
        if cprofile.mutCntTable is not None:
            pfile.write('#Mutation count table:\n')
            for base, mutCntDict in cprofile.mutCntTable.iteritems():
                pfile.write('%s' % base)
                for base, cnt in mutCntDict.iteritems():
                    pfile.write('\t%s\t%d' % (base, cnt))
                pfile.write('\n')
        pfile.write('#CIGAR lines: CIGAR pos mapq GCcontent quals\n')
        # Writing short CIGAR line, each in separate line
        for cline in cprofile.clines:
            pfile.write(cline.toTSVstring() + '\n')

        # TODO:
        # Writing information about qual and Ns distribution
        # ATM quals are stored within CLines
        pfile.close()
    except Exception:
        sys.stderr.write('\n\nError storing a profile!\n')


# Generate a sequence of random bases
# Used to generate clipped bases in a random read
def generate_randomseq(length):
    bases = 'ACGT'
    randomseq = ''.join(random.choice(bases) for i in xrange(length))
    return randomseq


# Mutate a given base
# Uses mutCntTable to calculate probabilities if available
def mutate(base, mutCntTable = None):
    # Shouldn't really happen, beacuse bases are taken from the reference
    # This first condition was added after the program tried to mutate only 'N' base in
    # Acinobacter Baumannii
    if base == 'N':
        return 'N'
    elif mutCntTable is None:
        bases = 'ACGT'
        # Remove a given base from the list, so its not generated
        bases.replace(base, '')

        # Return random base for a set (without the original base)
        return random.choice(bases)
    else:
        distrib = mutCntTable[base]
        sumcnt = sum(distrib.values())
        rand = random.randint(0, sumcnt-1)
        tmpsum = 0
        for base, cnt in distrib.items():
            tmpsum += cnt
            if tmpsum > rand:
                return base

        # Something is wrong
        sys.stderr.write('\nError generating mutation according to a given distribution!')
        return 'N'


if __name__ == "__main__":
    pass
