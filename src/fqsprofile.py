#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
import sys


# Class containg statistic information gained from a dataset
# Should be characteristic of a sequencing patform
# This class is used from importing FQSTQSim profiles
class FQSProfile:
    def __init__(self, name):
        self.name = name
        self.delCount = {}
        self.delsByRead = []
        self.delSize = {}
        self.insertCount = {}
        self.insertsByRead = []
        self.insertSize = {}
        self.mutationCount = {}
        self.mutationType = {}
        self.posCount = {}
        self.primerCheck = {}
        self.qualHist = {}
        self.readHist = {}


    # load delCount from .csv file
    def load_delCount(self, filepath):
        self.delCount = {}

        with open(filepath) as f:
            for line in f:
                # splitting each line using ',', there should be exactly two elements
                elements = line[:-1].split(',')         # disregarding \n at the end
                if len(elements) != 2:
                    sys.stderr.write('\n\nInvalid line in %s delCount file\n' % self.name)
                    exit(1)
                if int(elements[0]) in self.delCount.iterkeys():
                    sys.stderr.write('\n\nDuplicate element in %s delCount file\n' % self.name)
                    exit(1)

                self.delCount[int(elements[0])] = int(elements[1])

        return self.delCount


    def load_delsByRead(self, filepath):
        self.delsByRead = []

        with open(filepath) as f:
            # there should be only one line in the file, so I'm looking only at the first line
            line = f.readline()[:-1]     # removing \n from the end

            # splitting each line using ','
            elements = line.split(',')
            for element in elements:
                self.delsByRead.append(int(element))

        return self.delsByRead


    def load_delSize(self, filepath):
        self.delSize = {}

        with open(filepath) as f:
            for line in f:
                # splitting each line using ',', there should be exactly three elements
                elements = line.split(',')
                if len(elements) != 3:
                    sys.stderr.write('\n\nInvalid line in %s delSize file\n' % self.name)
                    import pdb
                    pdb.set_trace()
                    exit(1)
                if int(elements[0]) in self.delSize.iterkeys():
                    sys.stderr.write('\n\nDuplicate element in %s delSize file\n' % self.name)
                    exit(1)

                self.delSize[int(elements[0])] = [int(elements[1]), int(elements[2])]

        return self.delSize


    def load_insertCount(self, filepath):
        self.insertCount = {}

        with open(filepath) as f:
            for line in f:
                # splitting each line using ',', there should be exactly two elements
                elements = line.split(',')
                if len(elements) != 2:
                    sys.stderr.write('\n\nInvalid line in %s insertCount file\n' % self.name)
                    exit(1)
                if int(elements[0]) in self.insertCount.iterkeys():
                    sys.stderr.write('\n\nDuplicate element in %s insertCount file\n' % self.name)
                    exit(1)

                self.insertCount[int(elements[0])] = int(elements[1])

        return self.insertCount


    def load_insertsByRead(self, filepath):
        self.insertsByRead = []

        with open(filepath) as f:
            # there should be only one line in the file, so I'm looking only at the first line
            line = f.readline()[:-1]     # removing \n from the end

            # splitting each line using ','
            elements = line.split(',')
            for element in elements:
                self.insertsByRead.append(int(element))

        return self.insertsByRead


    def load_insertSize(self, filepath):
        self.insertSize = {}

        with open(filepath) as f:
            for line in f:
                # splitting each line using ',', there should be exactly three elements
                elements = line.split(',')
                if len(elements) != 3:
                    sys.stderr.write('\n\nInvalid line in %s insertSize file\n' % self.name)
                    exit(1)
                if int(elements[0]) in self.insertSize.iterkeys():
                    sys.stderr.write('\n\nDuplicate element in %s insertSize file\n' % self.name)
                    exit(1)

                self.insertSize[int(elements[0])] = [int(elements[1]), int(elements[2])]

        return self.insertSize


    def load_mutationCount(self, filepath):
        self.mutationCount = {}

        with open(filepath) as f:
            for line in f:
                # splitting each line using ',', there should be exactly two elements
                elements = line.split(',')
                if len(elements) != 2:
                    sys.stderr.write('\n\nInvalid line in %s mutationCount file\n' % self.name)
                    exit(1)
                if int(elements[0]) in self.mutationCount.iterkeys():
                    sys.stderr.write('\n\nDuplicate element in %s insertCount file\n' % self.name)
                    exit(1)

                self.mutationCount[int(elements[0])] = int(elements[1])

        return self.mutationCount


    # mutation type dictionary defines how likely is it for each base
    # to mutate to another base
    # this is stored in a shallow hierary of dictionaries, code shouldn't be too complex
    def load_mutationType(self, filepath):
        self.mutationType = {}

        with open(filepath) as f:
            # there should be exactly 5 lines in the file
            for line in f:
                # splitting each line using ',', there should be exactly two elements
                elements = line.split(',')
                # each line should contain exactly 9 elements, first defines starting base
                # and the following 4 pairs define final base and likelyhood of corresponding mutation
                if len(elements) != 9:
                    sys.stderr.write('\n\nInvalid line in %s mutationType file\n' % self.name)
                    exit(1)
                if elements[0] in self.mutationCount.iterkeys():
                    sys.stderr.write('\n\nDuplicate element in %s mutationType file\n' % self.name)
                    exit(1)

                mutationDict = {}
                mutationDict[elements[1]] = int(elements[2])
                mutationDict[elements[3]] = int(elements[4])
                mutationDict[elements[5]] = int(elements[6])
                mutationDict[elements[7]] = int(elements[8])
                self.mutationType[elements[0]] = mutationDict

        return self.mutationType


    def load_posCount(self, filepath):
        self.posCount = {}

        with open(filepath) as f:
            for line in f:
                # splitting each line using ',', there should be exactly two elements
                elements = line.split(',')
                if len(elements) != 2:
                    sys.stderr.write('\n\nInvalid line in %s posCount file\n' % self.name)
                    exit(1)
                if int(elements[0]) in self.posCount.iterkeys():
                    sys.stderr.write('\n\nDuplicate element in %s posCount file\n' % self.name)
                    exit(1)

                self.posCount[int(elements[0])] = int(elements[1])

        return self.posCount


    def load_primerCheck(self, filepath):
        self.primerCheck = {}

        with open(filepath) as f:
            for line in f:
                # splitting each line using ',', there should be exactly two elements
                elements = line.split(',')
                if len(elements) != 2:
                    sys.stderr.write('\n\nInvalid line in %s primerCheck file\n' % self.name)
                    exit(1)
                if elements[0] in self.primerCheck.iterkeys():
                    sys.stderr.write('\n\nDuplicate element in %s primerCheck file\n' % self.name)
                    exit(1)

                self.primerCheck[elements[0]] = int(elements[1])

        return self.primerCheck


    def load_qualHist(self, filepath):
        self.qualHist = {}

        with open(filepath) as f:
            for line in f:
                # splitting each line using ',', there is a variable number of elements
                # atm I'm not sure what each of them means
                # first element is the dictionary key, the reast are the value
                elements = line.split(',')
                if len(elements) < 2:
                    sys.stderr.write('\n\nInvalid line in %s qualHist file\n' % self.name)
                    exit(1)
                if int(elements[0]) in self.qualHist.iterkeys():
                    sys.stderr.write('\n\nDuplicate element in %s qualHist file\n' % self.name)
                    exit(1)

                self.qualHist[int(elements[0])] = elements[1:]

        return self.qualHist



    def load_readHist(self, filepath):
        self.readHist = {}

        with open(filepath) as f:
            for line in f:
                # splitting each line using ',', there should be exactly two elements
                elements = line.split(',')
                if len(elements) != 2:
                    sys.stderr.write('\n\nInvalid line in %s readHist file\n' % self.name)
                    exit(1)
                if int(elements[0]) in self.readHist.iterkeys():
                    sys.stderr.write('\n\nDuplicate element in %s readHist file\n' % self.name)
                    exit(1)

                self.readHist[int(elements[0])] = int(elements[1])

        return self.readHist

    def printStatistics(self, title = ''):

        sys.stdout.write('\n\n')
        if title != '':
            sys.stdout.write(title + '\n')
        sys.stdout.write('Printing statistics for %s sequencing platform:\n' % self.name)
        sys.stdout.write('DelCount dictionary contains %d enteries!\n' % len(self.delCount))
        sys.stdout.write('DelsByRead list contains %d enteries!\n' % len(self.delsByRead))
        sys.stdout.write('DelSize dictionary contains %d enteries!\n' % len(self.delSize))
        sys.stdout.write('InsertCount dictionary contains %d enteries!\n' % len(self.insertCount))
        sys.stdout.write('InsertsByRead list contains %d enteries!\n' % len(self.insertsByRead))
        sys.stdout.write('InsertSize dictionary contains %d enteries!\n' % len(self.insertSize))
        sys.stdout.write('MutationCount dictionary contains %d enteries!\n' % len(self.mutationCount))
        sys.stdout.write('MutationType dictionary contains %d enteries!\n' % len(self.mutationType))
        sys.stdout.write('PosCount dictionary contains %d enteries!\n' % len(self.posCount))
        sys.stdout.write('PrimerCheck dictionary contains %d enteries!\n' % len(self.primerCheck))
        sys.stdout.write('QualHist dictionary contains %d enteries!\n' % len(self.qualHist))
        sys.stdout.write('ReadHist dictionary contains %d enteries!\n' % len(self.readHist))


    # Load all statistics data from a folder with corresponding .CSV files
    def load_all(self, rootfolder, name=''):

        if name != '':
            self.name = name
        else:
            name = self.name

        # delCount
        filepath = os.path.join(rootfolder, name, name + '_delCount.csv')
        self.load_delCount(filepath)

        # delsByRead
        filepath = os.path.join(rootfolder, name, name + '_delsByRead.csv')
        self.load_delsByRead(filepath)

        # delSize
        filepath = os.path.join(rootfolder, name, name + '_delSize.csv')
        self.load_delSize(filepath)

        # insertCount
        filepath = os.path.join(rootfolder, name, name + '_insertCount.csv')
        self.load_insertCount(filepath)

        # insertsByRead
        filepath = os.path.join(rootfolder, name, name + '_insertsByRead.csv')
        self.load_insertsByRead(filepath)

        # insertSize
        filepath = os.path.join(rootfolder, name, name + '_insertSize.csv')
        self.load_insertSize(filepath)

        # mutationCount
        filepath = os.path.join(rootfolder, name, name + '_mutationCount.csv')
        self.load_mutationCount(filepath)

        # mutationType
        filepath = os.path.join(rootfolder, name, name + '_mutationType.csv')
        self.load_mutationType(filepath)

        # posCount
        filepath = os.path.join(rootfolder, name, name + '_posCount.csv')
        self.load_posCount(filepath)

        # primerCheck
        filepath = os.path.join(rootfolder, name, name + '_primerCheck.csv')
        self.load_primerCheck(filepath)

        # qualHist
        filepath = os.path.join(rootfolder, name, name + '_qualHist.csv')
        self.load_qualHist(filepath)

        # readHist
        filepath = os.path.join(rootfolder, name, name + '_readHist.csv')
        self.load_readHist(filepath)


    def generate_readsize(self):
        total_count = 0

        for value in self.readHist.itervalues():
            total_count += value




    # Generates a single read for the current profile
    def generate_read(self, reference):
        sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
        exit(1)
        pass


    # Generates a number of reads for the current profile
    def generate_reads_bynumber(self, number, reference):
        sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
        exit(1)
        pass


    # Generates a number of reads for the current profile, resulting in a given coverage
    def generate_reads_bycoverage(self, coverage, reference):
        sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
        exit(1)
        pass



def load_profile(profilefolder):

    # extracting last folder name
    name = os.path.basename(os.path.normpath(profilefolder))
    rootfolder = os.path.dirname(os.path.normpath(profilefolder))

    profile = FQSProfile(name)
    profile.load_all(rootfolder)
    # profile.printStatistics('After loading')

    return profile


def load_folder(profilesfolder):
    profiles = []

    filenames = os.listdir(profilesfolder)
    sys.stdout.write('\n\n')
    sys.stdout.write('Loading profiles from %s\n' % profilesfolder)
    for filename in filenames:
        filepath = os.path.join(profilesfolder, filename)
        if os.path.isdir(filepath):
            sys.stdout.write('Loading profile %s\n' % filename)
            profile = load_profile(filepath)
            profiles.append(profile)





def verbose_usage_and_exit():
    sys.stderr.write('FQSProfile - loading FASTQSim profile stats.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tload_profile\n')
    sys.stderr.write('\t\tload_folder\n')
    sys.stderr.write('\n')
    exit(0)


if __name__ == "__main__":
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'load_profile'):
        if (len(sys.argv) != 3):
            sys.stderr.write('Load a single FASTQSim profile from a given folder.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <profile folder>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        profilefolder = sys.argv[2]
        load_profile(profilefolder)

    elif (mode == 'load_folder'):
        if (len(sys.argv) != 3):
            sys.stderr.write('Load multiple FASTQSim profiles from a root profile folder.\n')
            sys.stderr.write('Each profile has a separate folder within root folder.\n')
            sys.stderr.write('Profile name is defined by folder name.\n')
            sys.stderr.write('\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <root profile folder>' % (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            exit(1)

        rootfolder = sys.argv[2]
        load_folder(rootfolder)

    else:
        sys.stderr.write('Unsupported mode parameter.\n\n')
        verbose_usage_and_exit()
