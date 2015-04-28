#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
import sys

# Using random number generator
import random

import numpy


# A class representing a discrete probability distribution function
# Main element is a dictionary whose keys represent values, and dictionary
# values represent probabilities that a random value is <= given value (dict key)
# Keys range from minValue (inclusive) to maxValue (exclusive)
# All keys in the range must be present.

class Probability_Distribution:
    def __init__(self, density = None, Gauss = False, Poisson = False):
        if density is not None:
            self.createFromDensity(density)
        elif Gauss:
            self.createGauss(0, 100, 50, 10)
        elif Poisson:
            self.createPoisson(0, 100, 50, 10)
        else:
            self.distrib = {}
            self.minValue = 0
            self.maxValue = 0
            self.normalized = False

    # Number of enteries in distribution dictionary
    def numPoints(self):
        return len(self.distrib)

    # Creates distribution function from density function
    # Density must be a dictionary, dict keys are values and dict keys are probabilities or number of occurences for a given value
    def createFromDensity(self, density):
        self.distrib = {}
        initial = True
        self.minValue = 0
        self.maxValue = 0
        self.normalized = True

        # Determining a minimum and maximum value
        for key in density.iterkeys():
            if initial:
                initial = False
                self.minValue = key
                self.maxValue = key
            else:
                if key > self.maxValue:
                    self.maxValue = key
                if key < self.minValue:
                    self.minValue = key

        # Calculating the distribution
        lastValue = self.minValue -1
        lastProb = 0
        for key in sorted(density.keys()):
            prob = density[key] + lastProb          # Adding last prob to make it a distribution function
            valuerange = key - lastValue
            probrange = prob - lastProb
            for i in xrange(1, valuerange+1):
                self.distrib[lastValue+i] = lastProb + i*float(probrange)/valuerange
            lastValue = key
            lastProb = prob

        self.normalize()


    def createGauss(self, min, max, avg, stdev):
        self.distrib = {}
        # to supress warnings
        tsum = min + max + avg + stdev
        sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
        # exit(1)
        # pass


    def createPoisson(self, min, max, avg, stdev):
        self.distrib = {}
        # to supress warnings
        tsum = min + max + avg + stdev
        sys.stderr.write('Function %s not implemented yet!\n\n' % (sys._getframe().f_code.co_name))
        # exit(1)
        # pass


    def printStatistics(self, title = ''):
        sys.stdout.write('\n\nPrinting statistics for distribution %s\n' % title)
        sys.stdout.write('Number of enteries: %d \n' % len(self.distrib))
        sys.stdout.write('Valuerange: [%d, %d] \n' % (self.minValue, self.maxValue))


    def normalize(self):
        maxProb = 0

        for prob in self.distrib.itervalues():
            if prob > maxProb:
                maxProb = prob

        for (key, value) in self.distrib.iteritems():
            self.distrib[key] /= maxProb

        self.normalized = True


    def generateRandom(self):
        if not self.normalized:
            sys.stderr.write('\nDistribution not normalized! Normalize before generating random numbers!\n')
        else:
            # lastProb = 0
            randnum = random.random()
            for key in sorted(self.distrib.keys()):
                # if self.distrib[key] > randnum and lastProb < randnum:
                if self.distrib[key] > randnum:
                    return key
                # lastProb = self.distrib[key]



# A function that calculates average and standard deviation of a list of numbers
def avg_and_stdev(numlist):
    avg = stdev = 0

    # Avg
    avg = numpy.mean(numlist)

    # StDev
    stdev = numpy.std(numlist)

    return avg, stdev



if __name__ == "__main__":
    pass
