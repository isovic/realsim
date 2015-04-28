#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
import sys

import random
import numpy

import fastqparser
from fqsprofile import FQSProfile
from utility_distrib import Probability_Distribution, avg_and_stdev

# Class containg statistic information gained from a dataset
# Use-ready for generating new reads
# Should be characteristic of a sequencing patform
class RSProfile:
    def __init__(self, name):
        self.name = name

        ### DELETEs
        # ATM delete size distribution is uniform across the whole read

        # Distribution of positions of deleted bases within a read, used to generate random delete position
        self.delPosDistrib = None

        # Average delete count and its standard deviation, used to generate a random delete count (number of deletes in a read)
        self.delCntAvg = 0
        self.delCntStDev = 0

        # Delete size distribution, used to generate random delete size
        self.delSizeDistrib = None

        ### INSERTs
        # ATM insert size distribution is uniform across the whole read

        # Distribution of positions of inserted bases within a read, used to generate random insert position
        self.insPosDistrib = None

        # Average insert count and its standard deviation, used to generate a random insert count (number of inserts in a read)
        self.insCntAvg = 0
        self.insCntStDev = 0

        # Insert size distribution, used to generate random insert size
        self.insSizeDistrib = None

        ### MUTATIONs

        # Distribution of positions of mutated bases withn a read, used to generate random mutation position
        # and, together with
        self.mutDistrib = None

        self.mutType = None

        # ATM not sure what to do with these:
        # self.posCount
        # self.primerCheck

        # Qual value distribution acros positions within a read
        # It should have average qual and standard deviation for each read position
        # TODO: Consider how to store this, atm not using this
        self.qualDistrib = None

        # Normalized read size distribution
        # keys should contain all values between 0 and maxReadLength
        self.readDistrib = None


    def import_fqsprofile(self, fqsprofile):
        self.name = fqsprofile.name

        # Distribution of positions of deletes in a read
        self.delPosDistrib = Probability_Distribution(fqsprofile.delCount)

        # Distribution of numbers of deletes in a read
        self.delCntAvg, self.delCntStDev = avg_and_stdev(fqsprofile.delsByRead)

        # Distribution of delete sizes
        # taking only the first delSize value (there are 2 in fqsprofile)
        tmpDelSize = {key:value[0] for (key,value) in fqsprofile.delSize.iteritems()}
        self.delSizeDistrib = Probability_Distribution(tmpDelSize)

        # Distribution of positions of inserts in a read
        self.insPosDistrib = Probability_Distribution(fqsprofile.insertCount)

        # Distribution of numbers of inserts in a read
        self.insCntAvg, self.insCntStDev = avg_and_stdev(fqsprofile.insertsByRead)

        # Distribution of insert sizes
        # taking only the first insertSize value (there are 2 in fqsprofile)
        tmpInsSize = {key:value[0] for (key,value) in fqsprofile.insertSize.iteritems()}
        self.insSizeDistrib = Probability_Distribution(tmpInsSize)

        # Distribution of positions of mutations in a read
        self.mutDistrib = Probability_Distribution(fqsprofile.mutationCount)

        # Mutation type, keepint the same structure are in FQSprofile
        self.mutType = fqsprofile.mutationType

        # Keeping this as it is ATM
        self.qualDistrib = fqsprofile.qualHist

        # Read size distribution
        self.readDistrib = Probability_Distribution(fqsprofile.readHist)


    def generate_quals(self, read_length):
        sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
        # exit(1)
        # pass

    def generate_deletes(self, read):
        # Generate a number of deletes according to Gaussian distribution using numpy
        numdels = numpy.random.normal(self.delCntAvg, self.delCntStDev)
        sys.stdout.write('\nGenerating %d deletes' % int(numdels))

        # For the testing purposes atm, instead of deleting, replacing corresponding bases with 'D'
        for i in xrange(int(numdels)):
            # generating delete position within read
            position = self.delPosDistrib.generateRandom()
            # generating delete size
            size = self.delSizeDistrib.generateRandom()

            sys.stdout.write('\nDelete %d: position: %d, size: %d' % (i, position, size))

            # Replacing bases to be deleted with 'D'
            replace = 'D'*size
            read = read[:position] + replace + read[position+size:]

        # TODO: Actually delete bases

        return read

    def generate_inserts(self, read):
        # Generate a number of inserts according to Gaussian distribution using numpy
        numinss = numpy.random.normal(self.insCntAvg, self.insCntStDev)
        sys.stdout.write('\nGenerating %d inserts' % int(numinss))

        # For the testing purposes atm, instead of actual bases, inserting 'I's
        # Also, have to sort insert positions to avoid one insert affecting another's position

        inserts = []

        for i in xrange(int(numinss)):
            # generating insert position within read
            position = self.insPosDistrib.generateRandom()
            # generating insert size
            size = self.insSizeDistrib.generateRandom()
            sys.stdout.write('\nInsert %d: position: %d, size: %d' % (i, position, size))
            inserts.append([position, size])

        # Sorting inserts according to position in descending order
        # TODO: check if the same position was generated more than once

        # Sorting function
        def keyfun(x):
            return x[1]

        inserts.sort(reverse=True, key=keyfun)

        for insert in inserts:
            position = insert[0]
            size = insert[1]
            insstring = size*'I'
            read = read[:position ] + insstring + read[position:]

        # TODO: Actually insert bases

        return read

    def generate_mutations(self, read):

        # Generating exactly one mutation for each read
        # Putting 'M' instead of mutated base

        position = self.mutDistrib.generateRandom()

        read = read[:position] + 'M' + read[position+1:]

        return read


    # DEPRECATED
    # Generate random starting read
    # Generates random position within the reference genome, generates random read length,
    # and reads from reference genome (double the size of read length to allow for latter manipulation)
    def generate_startread(self, reference):
        readlen = self.generate_readlength()
        # working read length is twice the actual read length, to allow for a lot of deletes
        tmplen = readlen*2
        reflen = len(reference)
        startpos = self.generate_startpos(reflen)

        # in case the tmplength is larger then or almost(TODO) as large as the genome
        if tmplen >= reflen:
            return reference

        count = 1
        while startpos > reflen - tmplen:
            startpos =  self.generate_startpos(len(reference))
            count += 1
            if count > 20:
                sys.stderr.write('\nFailed to generate starting read!\n')
                exit(1)

        return reference[startpos:startpos+tmplen], readlen

    # Generate initial position within reference genome from where a random read will be generated
    def generate_startpos(self, genome_length):
        return int(random.random()*genome_length)

    # Generate randomized read length
    def generate_readlength(self):
        return self.readDistrib.generateRandom()

    # Generate random read
    # Generates random position within the reference genome, generates random read length,
    # and reads from reference genome (double the size of read length to allow for latter manipulation)
    # Then applies deletes, inserts and mutations
    # at the end, shortens the read to the generated read length
    def generate_read(self, reference):
        readlen = self.generate_readlength()
        # working read length is twice the actual read length, to allow for a lot of deletes
        tmplen = readlen*2
        reflen = len(reference)
        startpos = self.generate_startpos(reflen)

        # in case the tmplength is larger then or almost(TODO) as large as the genome
        if tmplen >= reflen:
            read = reference

        else:
            count = 1
            while startpos > reflen - tmplen:
                startpos =  self.generate_startpos(len(reference))
                count += 1
                if count > 20:
                    sys.stderr.write('\nFailed to generate starting read!\n')
                    exit(1)
            read = reference[startpos:startpos+tmplen]

        read = self.generate_deletes(read)
        read = self.generate_inserts(read)
        read = self.generate_mutations(read)

        # Returning only a geneated read length
        return read[:readlen]


    # A function to test if rsprofile is loaded correctly
    def test(self):
        sys.stdout.write('\nTesting profile %s!\n' % self.name)

        if (self.delPosDistrib is not None):
            sys.stdout.write('Delete probability distribution loaded. Contains %d elements.\n' % self.delPosDistrib.numPoints())

        sys.stdout.write('Distribution of number of deletes in a read. Avg=%f, StDev=%f\n' % (self.delCntAvg, self.delCntStDev))

        if (self.delSizeDistrib is not None):
            sys.stdout.write('Delete size probability distribution loaded. Contains %d elements.\n' % self.delSizeDistrib.numPoints())


        if (self.insPosDistrib is not None):
            sys.stdout.write('Insert probability distribution loaded. Contains %d elements.\n' % self.insPosDistrib.numPoints())

        sys.stdout.write('Distribution of number of inserts in a read. Avg=%f, StDev=%f\n' % (self.insCntAvg, self.insCntStDev))

        if (self.insSizeDistrib is not None):
            sys.stdout.write('Insert size probability distribution loaded. Contains %d elements.\n' % self.insSizeDistrib.numPoints())


        if (self.mutDistrib is not None):
            sys.stdout.write('Mutation probability distribution loaded. Contains %d elements.\n' % self.mutDistrib.numPoints())

        if (self.mutType is not None):
            sys.stdout.write('Mutation type probability distribution loaded. Contains %d elements.\n' % len(self.mutType))


        if (self.qualDistrib is not None):
            sys.stdout.write('Qual distribution loaded. Contains %d elements.\n' % len(self.qualDistrib))


        if (self.readDistrib is not None):
            sys.stdout.write('Read size distribution loaded. Contains %d elements.\n' % self.readDistrib.numPoints())

        for i in xrange(10):
            print self.generate_readlength()
        # print self.readDistrib.distrib


def verbose_usage_and_exit():
    sys.stderr.write('RSProfile - loading FASTQSim profile stats.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tload_FQSprofile\n')
    sys.stderr.write('\t\tload_FQSfolder\n')
    sys.stderr.write('\t\tgenerate_read\n')
    sys.stderr.write('\t\tgenerate_dataset\n')
    sys.stderr.write('\n')
    exit(0)


if __name__ == "__main__":
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'load_FQSprofile'):
        if (len(sys.argv) != 3):
            sys.stderr.write('Load a single FASTQSim profile from a given folder.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <profile folder>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        profilefolder = sys.argv[2]

        # extracting inner-most folder name
        name = os.path.basename(os.path.normpath(profilefolder))
        rootfolder = os.path.dirname(os.path.normpath(profilefolder))

        fqsprofile = FQSProfile(name)
        fqsprofile.load_all(rootfolder)

        rsprofile = RSProfile('')
        rsprofile.import_fqsprofile(fqsprofile)
        rsprofile.test()

    elif (mode == 'load_FQSfolder'):
        if (len(sys.argv) != 3):
            sys.stderr.write('Load multiple FASTQSim profiles from a root profile folder.\n')
            sys.stderr.write('Each profile has a separate folder within root folder.\n')
            sys.stderr.write('Profile name is defined by folder name.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <root profile folder>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

            sys.stderr.write('\nOption load_FQSfolder not completeley implemented yet!\n')


    elif (mode == 'generate_read'):
        if (len(sys.argv) != 4):
            sys.stderr.write('Generate a single read using a FASTQSim profile.\n')
            sys.stderr.write('FASTQSim profile is defined by folder name.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <profile folder> <reference>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        profilefolder = sys.argv[2]
        referencepath = sys.argv[3]
        referencefastq = fastqparser.read_fastq(referencepath)
    	referenceseq = referencefastq[1][0]

        # extracting inner-most folder name
        name = os.path.basename(os.path.normpath(profilefolder))
        rootfolder = os.path.dirname(os.path.normpath(profilefolder))

        fqsprofile = FQSProfile(name)
        fqsprofile.load_all(rootfolder)

        rsprofile = RSProfile('')
        rsprofile.import_fqsprofile(fqsprofile)

        read = rsprofile.generate_read(referenceseq)
        sys.stdout.write(read)
        sys.stdout.write('\n')

        # sys.stderr.write('\nOption generate_read not completeley implemented yet!\n')

    elif (mode == 'generate_dataset'):
        if (len(sys.argv) != 3):
            sys.stderr.write('Generate a complete dataset using a FASTQSim profile.\n')
            sys.stderr.write('FASTQSim profile is defined by folder name.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <profile folder> <reference> <number of reads>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        profilefolder = sys.argv[2]

        # extracting inner-most folder name
        name = os.path.basename(os.path.normpath(profilefolder))
        rootfolder = os.path.dirname(os.path.normpath(profilefolder))

        fqsprofile = FQSProfile(name)
        fqsprofile.load_all(rootfolder)

        rsprofile = RSProfile('')
        rsprofile.import_fqsprofile(fqsprofile)

        sys.stderr.write('\nOption generate_dataset not completeley implemented yet!\n')


    else:
        sys.stderr.write('Unsupported mode parameter.\n\n')
        verbose_usage_and_exit()
