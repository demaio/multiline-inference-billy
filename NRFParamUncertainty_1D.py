#!/usr/bin/python
# -*- coding: utf-8 -*-

# 1D NRF Parameter Uncertainty Calculator - Python 2.7
# Ruaridh Macdonald, MIT, 2016
# rmacd@mit.edu

# Uses NRFGamma class from:
# StandaloneNRFLineCalculator
# Jayson Vavrek, MIT, 2015
# jvavrek@mit.edu

print "1D NRF Parameter Uncertainty Calculator - Python 2.7"
print "Note: T_Debye still not properly implemented! Doppler-broadened widths will be off."

import sys
import numpy as np

import NRFmultiLine # Import file with functions to use here
from NRFGamma import NRFGamma # Import NRFGamma class for storing the line properties easily
# Takes (_z, _a, _Elevel, _Egamma, _Width, _prob, _GSprob, _J0, _Jr, _TDebye, _nDens, _sigmaNRLevel, _sigmaNRGamma, _counter)
# and returns the same as well as: (sigmaInt, sigmaDmax, sigmaNRLevel, sigmaNRGamma, alpha, Delta)

# Timing for debugging and comparisons
import time
startTime = time.time()

# Provide command line option functionality
import argparse
parser = argparse.ArgumentParser(description='Processes a standalone database of NRF gammas and returns lists of branced and neighbouring NRF lines which leak the most information, sorted by attenuation coefficient. \nExample:\
    .......................................................................\
    $./multilineSearch.py -neighE=0.01 -matList=matList.txt -branchOut=5 -neighOut=5 -EMin=1 -EMax=2.5 -bremsMin=2 -bremsMax=9 \
    .......................................................................\
    Produces a lists of the 5 branched and neighbouring pairs with the higher attenuation coefficients, with Elevel between 2 and 9 MeV and EGamma between 1 and 2.5 MeV ')
    
# -h and --help options exist by default
parser.add_argument('-neighE', help='Energy gap to qualify as ''neighbours'' (MeV), default = 1KeV ', type=float, default=0.001)
parser.add_argument('-matList', help='File name of list of isotopes to check and their number densities \n Format: A Z numDen*A', type=str)
parser.add_argument('-EMin', help='Detector minimum energy (MeV), default = 0 ', type=float, default=0)
parser.add_argument('-EMax', help='Detector maximum energy (MeV), default = 20 ', type=float, default=20)
parser.add_argument('-bremsMin', help='Photon sourve minimum energy (MeV), default = 0 ', type=float, default=0)
parser.add_argument('-bremsMax', help='Photon source maximum energy (MeV), default = 20 ', type=float, default=20)
parser.add_argument('-neighOut', help='Lists -neighOut worst neighbouring pairs, ranked by cross section, default = 0 ', type=int)
parser.add_argument('-branchOut', help='Lists -branchOut worst branched pairs, ranked by cross section, default = 0 ', type=int)
parser.add_argument('-plotOn', help='Do you want to plot the results, Yes = 1 ', default = 0, type=int)
args = parser.parse_args()

plotOn = args.plotOn

# Check that material list was given by user
if args.matList == None: sys.exit("User must specify a material file describing the foil and warhead isotopic content")
# Check that definition of 'neighbour' is positive
if args.neighE<= 0 : 
    sys.exit("Neighbour energy gap must be >= 0")
else: neighE = args.neighE
# Check that energy ranges make sense
if args.EMax<=args.EMin or args.EMax<=0: sys.exit("\nMaximum detector energy must be greater than minimum detector energy and positive")
if args.bremsMax<=args.bremsMin or args.bremsMax<=0: sys.exit("\nMaximum Bremsstrahlung energy must be greater than minimum Bremsstrahlung enery and positive")

# ------ ------ ------
# Create list of isotopes to work with and the correspondign number densities
matList, nDensList, thickList = NRFmultiLine.parse_materials(args.matList)
    
# Error checking on material lists (though this will normally be caught during reading)
if len(matList) != len(nDensList): sys.exit("Material input incorrect: different numbers of isotopes and number densities in file")
if len(matList) == 0 : sys.exit("Material input incorrect: must contain at least one material")
    
# Build a matrix of all the non-resonant cross sections up front for use in loops later
NRData = np.loadtxt('nonResonantAttenuation.txt',dtype=float,delimiter='|') # This file has to be carefully formatted
NRData = np.split(NRData,[1],axis=1)
NREnergy = NRData[0] # Vector of energies that the database uses
NRData = NRData[1]   # Non-resonant attenuation data by isotope, z = 1:100
    
# ------ ------ ------
# Import data (read-only) from standalone.dat database of gammas
f = open('standalone.dat','r')
print f
NRfile = open('nonResonantAttenuation.txt','r')
print NRfile

# ------ ------ ------
# Print go statement    
print "\nLooking for significant pairs of lines emitted between %s and %s MeV \nUsing Bremsstrahlung source from %s to %s MeV" %("{:6.3f}".format(args.EMin), "{:6.3f}".format(args.EMax),"{:6.3f}".format(args.bremsMin), "{:6.3f}".format(args.bremsMax))

# ------ ------ ------
# Begin main thread 

# List to contain the NRF lines
emitList = []
counter = -1 # Counter for managing indexing in the final list of NRF lines

f.seek(0)
for line in f: # Running through each line of the NRF database to check if it matches our parameters
    line = line.split(" ")

    # Grab data and convert strings to ints, floats
    [z, a] = map(int, [line[0], line[1]]);
    [Elevel, Egamma, Width, prob, GSprob, J0, Jr, TDebye] = map(float, [line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9]])
    
    # Exclude gammas without valid data, outside the energy ranges or not on the material list
    if Width>0 and ~np.isinf(Width) and prob>0 and [z,a] in matList and Egamma>args.EMin and Egamma<args.EMax and Elevel>args.bremsMin and Egamma<args.bremsMax:
        counter += 1 # Add another accepted line to the count
        
        nDens = nDensList[matList.index([z,a])]       # Calculate number density * 1e-24 
        
        # Search and interpolate for the non-resonant cross section    
        # Data file has one column of energy labels and then isotopes 1:100
        levelIndex = NRFmultiLine.find_nearestE(Elevel,NREnergy)    # Find energy index of the non-resonant cross sections
        gammaIndex = NRFmultiLine.find_nearestE(Egamma,NREnergy)    # Find energy index of the non-resonant cross sections
        
        sigmaNRLevelWarhead = 0
        sigmaNRGammaWarhead = 0
        sigmaNRLevelFoil = 0
        sigmaNRGammaFoil = 0
        for i in range(len(matList)) :
            sigmaNRLevelFoil += NRData[levelIndex,matList[i][0]]*nDensList[i][1]    # Non-resonant attenuaton from all isotopes at resonance energy
            sigmaNRGammaFoil += NRData[gammaIndex,matList[i][0]]*nDensList[i][1]    # Non-resonant attenuaton from all isotopes at emitted gamma energy
            sigmaNRLevelWarhead += NRData[levelIndex,matList[i][0]]*nDensList[i][0] # Non-resonant attenuaton from all isotopes at resonance energy
            sigmaNRGammaWarhead += NRData[gammaIndex,matList[i][0]]*nDensList[i][0] # Non-resonant attenuaton from all isotopes at emitted gamma energy   
        sigmaNRLevel = [sigmaNRLevelWarhead,sigmaNRLevelFoil]
        sigmaNRGamma = [sigmaNRGammaWarhead,sigmaNRGammaFoil]
        
        # Build NRFGamma instance
        x = NRFGamma(z,a,Elevel,Egamma,Width,prob,GSprob,J0,Jr,TDebye,nDens,thickList[matList.index([z,a])],sigmaNRLevel,sigmaNRGamma,counter)          
        emitList.append(x)

# ------ ------ ------
# Calculate counts at detector given the top lines
# Currently we find the highest intensity NRF peak for each isotope and set the smallest one to be 1e4
# This is then used to scale all the other lines

# Find the largest line for each isotope
maxLineCount = [0.0]*len(matList)
maxLineIndx = [0]*len(matList)
for line in emitList:
    if line.counts > maxLineCount[matList.index([line.z,line.a])]:
        maxLineIndx[matList.index([line.z,line.a])] = line.index
        maxLineCount[matList.index([line.z,line.a])] = line.counts
        
# Find smallest of these max lines and set source strength so that this line has 1e4 counts
minMaxLineIndex = maxLineCount.index(min(maxLineCount))
sourceStrength = 1e4/maxLineCount[minMaxLineIndex]

# Update the count values with the source strength
for line in emitList:
    line.counts = line.counts*sourceStrength
    
# For each line, the number of counts is calculated in the NRFGamma class, assuming Bremms flux = 1 but this is fine for our comparisons (but it makes it trickier to order them)

# ------ ------ ------
# Compare branched and neighbouring pairs vs resonance energy, minimum sigma_NRF and total sigma_NRF
branchData, neighData = NRFmultiLine.branchNeigh_compare(emitList,neighE,args.branchOut,args.neighOut,args.plotOn)

# The data lists returned by branchNeigh_compare contain:
#    [0] : The indexes of the branched / neighbouring pairs from emitList, grouped as pairs
#    [1] : The alpha values for the pairs, grouped as pairs
#    [2] : Resonance energies for the lines in the pairs, grouped as pairs for neigh and as a single for branch
#    [3] : Minimum number of counts for the items in a pair
#    [4] : Ordered index of the entries by number of counts, largest -> smallest

# ------ ------ ------
# Use numerical methods to estimate properties of interest

# ------ ------ ------
# Close files and tidy up

f.close()
NRfile.close()

# ------ ------ ------
print "Took: %s seconds" %"{:5.3}".format(time.time()-startTime)