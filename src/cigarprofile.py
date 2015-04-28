#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
import sys

# To generate random reads
import random

# For manipulating CIGAR strings
import re

# Mainly for callculating a reverse complement of a sequence
# and other utility stuff
import fastqparser


CIGAR_OPERATIONS_ALL = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'];
CIGAR_OPERATIONS_BASIC = ['M', 'I', 'D', 'S', 'H'];
CIGAR_OPERATIONS_EXTENDED = ['M', 'I', 'D', 'S', 'H', '=', 'X'];


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
    def __init__(self, cigar = '', position = -1, quality = -1):
        self.cigar = cigar
        self.position = position
        self.quality = quality

    def toTSVstring(self):
        return ('%s\t%d\t%d' % (self.cigar, self.position, self.quality))

    def fromTSVstring(self, line):
        # If present, removing \n from the end of line
        if line[-1] == '\n':
            line = line[:-1]
        elements = line.split('\t')
        self.cigar = elements[0]
        self.position = int(elements[1])
        self.quality = int(elements[2])

    # Splits the CIGAR string into individual operations, in the
    # same format as the original CIGAR string is in.
    # The result is returned as an array of a 2-element tuples, e.g.
    # [(12, 'M'), (3, 'D'), (4, 'M')].
    def splitCigar(self):
        # Using regular expressions to find repeating digit and skipping one character after that
        pattern = '(\d+)(.)'
        operations = re.findall(pattern, self.cigar)

        return operations

    # Returns a length of the read that generated the cigar string
    def length(self):
        return sum([int(x[0]) for x in self.splitCigar()])


class CIGARProfile:
    def __init__(self, name):
        self.name = name
        self.clines = []

        # ATM infromation about quals Ns is stored for the whole profile
        self.qualDistrib = None
        self.NDistrib = None

    def appendCLine(self, cline):
        if cline is None:
            raise Exception('Trying to add nonexisting CIGAR line! (%d)' % len(self.clines))
        self.clines.append(cline)

    def getRandomCigar(self):
        rand = random.randint(0, len(self.clines)-1)
        return self.clines[rand]

    # ATM using all cigar lines to generate reads at random positions
    def generateRandomRead(self, reference):
        cline = self.getRandomCigar()
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
                    read += mutate(startread[pos+i])
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

        # For testing purposes returning a generated read together with the used CIGAR string
        return read, cline.cigar


    # Separate function to generate a whole dataset faster and with less memory
    # (instead of calling generateRandomRead many times)
    def generateRandomReadsByNumber(self, reference, numreads):
        reads = []
        reflen = len(reference)

        # Generating reads
        for i in xrange(numreads):
            cline = self.getRandomCigar()
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
                        read += mutate(startread[pos+i])
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

            reads.append(read)

        return reads


    # Separate function to generate a whole dataset faster and with less memory
    # (instead of calling generateRandomRead many times)
    def generateRandomReadsByCoverage(self, reference, coverage):
        reads = []
        reflen = len(reference)

        numbases = 0
        numreads = 0

        while True:
            cline = self.getRandomCigar()
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
                        read += mutate(startread[pos+i])
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

            reads.append(read)
            numbases += len(read)
            numreads += 1

            # Finish if a set coverage is reached
            if numbases/reflen > coverage:
                break

        return reads


def loadCProfile(filepath):

    cprofile = None
    with open(filepath, 'r') as pfile:
        # loading profile name from the first line
        line = pfile.readline()
        # Removing # from the start and \n from the end of line
        name = line[1:-1]
        cprofile = CIGARProfile(name)
        for line in pfile:
            cline = CIGARLine()
            cline.fromTSVstring(line)
            cprofile.appendCLine(cline)

    # TODO:
    # Load information about quals and Ns distribution

    return cprofile


def storeCProfile(filepath, cprofile):

    with open(filepath, 'w') as pfile:
        # Writing profile name in the first line with # in front
        pfile.write('#%s\n' % cprofile.name)
        # Writing short CIGAR line, each in separate line
        for cline in cprofile.clines:
            pfile.write(cline.toTSVstring() + '\n')


    # TODO:
    # Writing information about qual and Ns distribution


# Generate a sequence of random bases
# Used to generate clipped bases in a random read
def generate_randomseq(length):
    bases = 'ACGT'
    randomseq = ''.join(random.choice(bases) for i in xrange(length))
    return randomseq


# Mutate a given base
# ATM the probability is equal for all bases
def mutate(base):
    bases = 'ACGT'
    # Remove a given base from the list, so its not generated
    bases.replace(base, '')

    # Return random base for a set (without the original base)
    return random.choice(bases)


if __name__ == "__main__":
    pass
