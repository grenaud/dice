#!/usr/bin/env python

import subprocess
from optparse import OptionParser
import sys, os, random
import numpy as np
from scipy import stats
from collections import Counter,defaultdict
from operator import itemgetter 
from random import randint

parser = OptionParser("$prog [options]")
parser.add_option("-t", "--theta", dest="theta", help="Theta", default=20, type="float")
parser.add_option("-s", "--timesplit", dest="timesplit", help="Split time in 2N_0 generations", default=0.05, type="float")
parser.add_option("-f", "--destfolder", dest="destfolder", help="Destination folder", default=None, type="string")
parser.add_option("-n", "--numsim", dest="numsim", help="Number of simulations", default=100, type="int")
parser.add_option("-c", "--contrate", dest="contrate", help="Contamination rate", default=0.0, type="float")
parser.add_option("-e", "--errorrate", dest="errorrate", help="Error rate", default=0.0, type="float")
parser.add_option("-m", "--meancoverage", dest="meancoverage", help="Mean coverage", default=30.0, type="float")
parser.add_option("-o", "--outfile", dest="outfile", help="Output file", default=None, type="string")
parser.add_option("-a", "--admixrate", dest="admixrate",help="Admixture rate", default=0, type="float")
parser.add_option("-b", "--admixtime", dest="admixtime", help="Admixture time", default=0.0045, type="float")
parser.add_option("-j", "--numsamphumA", dest="numsamphumA", help="Number of sampled high-coverage humans A", default=6, type="int")
parser.add_option("-k", "--numsamphumB", dest="numsamphumB", help="Number of sampled high-coverage humans B", default=2, type="int")
(options,args) = parser.parse_args()

theta = options.theta
timesplit = options.timesplit
destfolder = options.destfolder
numsim = options.numsim
meancoverage = options.meancoverage
errorrate = options.errorrate
contrate = options.contrate
admixrate = options.admixrate
admixtime = options.admixtime
outfile = open(options.outfile,"w")

# Convert to units of 4N_0 generations
timesplitms = timesplit / 2.0
firstsplit = float(3200) / float(40000)
numsamphumA = options.numsamphumA
numsamphumB = options.numsamphumB
numtotalhumA = 100
numtotalhumB = 100
numarch = 2

puredict = defaultdict(int)
i=1
#if True:
try:
	while i < (numsim+1):

		# Create simulation file
		infile_name = "/home/fernando_racimo/TwoPopCont/"+destfolder+"/simul_"+str(i)+"_"+str(randint(0,10000))+".txt"
	        if admixrate > 0:
			commname = "/home/fernando_racimo/bin/msms/bin/msms "+str(numtotalhumA+numtotalhumB+numarch)+" 1 -t "+str(theta)+" -I 3 "+str(numtotalhumA)+" "+str(numtotalhumB)+" "+str(numarch)
			commname = commname +" -es "+str(admixtime)+" 1 "+str(1 - admixrate)+" -ej "+str(admixtime+0.00001)+" 4 3"+" -ej "+str(firstsplit)+" 1 2"+" -ej "+str(timesplitms)+" 2 3 | tail -n+7"
		else:
			commname = "/home/fernando_racimo/bin/msms/bin/msms "+str(numtotalhumA+numtotalhumB+numarch)+" 1 -t "+str(theta)+" -I 3 "+str(numtotalhumA)+" "+str(numtotalhumB)+" "+str(numarch)
			commname = commname +" -ej "+str(firstsplit)+" 1 2"+" -ej "+str(timesplitms)+" 2 3 | tail -n+7"
	#	print commname
	#	i += 1
		commname = commname + " | sed 's/1/1,/g' | sed 's/0/0,/g' | sed 's/,$//'"
		
		commname = commname + " > "+infile_name
		v = subprocess.Popen(commname,shell=True)
		v.communicate()
		
		
		try:
			
			print "simul number "+str(i)
			
			# Load and sum matrices
			sitemat = np.loadtxt(infile_name,delimiter=",",dtype="int")

			# Store site data
			for column in sitemat.transpose():
				modsampA = np.sum(column[0:numsamphumA])
				modfreqA = float(np.sum(column[0:numtotalhumA])) / float(numtotalhumA)
				modsampB = np.sum(column[numtotalhumA:(numtotalhumA+numsamphumB)])
				modfreqB = float(np.sum(column[numtotalhumA:(numtotalhumA+numtotalhumB)])) / float(numtotalhumB)
				arch = np.sum(column[(numtotalhumA+numtotalhumB):(numtotalhumA+numtotalhumB+numarch)])

			
				# Require the site to be segregating in sample
				if ( (modsampA+modsampB) > 0 and (modsampA+modsampB) < (numsamphumA+numsamphumB) ):
		       	                #if True:
					puredict[(float(modsampA)/numsamphumA,float(modsampB)/numsamphumB,arch)] += 1
			i += 1
					
		except:
			continue

		os.remove(infile_name)

except:
	print "Script interrupted"




countdict = defaultdict(int)
i = 1
totalsites = sum(puredict.values())
#totalsites = sum(puredict.keys())
print >>outfile, "Anc\tDer\tPanelFreqA\tPanelFreqB\tNum"
for key in puredict.keys():
	panelfreqB = float(key[1])
	panelfreqA =  float(key[0])
	archgeno = int(key[2])
	panelfreq = panelfreqA
#	print "Key: "+str(key)
#	print "Number: "+str(puredict[key])
	coveragevec = np.random.poisson(meancoverage, puredict[key])
#	print coveragevec
	for sitecoverage in coveragevec:
		if(archgeno == 2):
			q = contrate*panelfreq*(1-errorrate) + contrate*(1-panelfreq)*errorrate + (1-contrate)*(1-errorrate)
		elif(archgeno == 1):
			q = contrate*panelfreq*(1-errorrate) + contrate*(1-panelfreq)*errorrate + (1-contrate)*(1-errorrate)/2.0 + (1-contrate)*errorrate/2.0
		elif(archgeno == 0):
			q = contrate*panelfreq*(1-errorrate)+contrate*(1-panelfreq)*errorrate + (1-contrate)*errorrate
		numder = np.random.binomial(sitecoverage,q,1)[0]
		numanc = sitecoverage - numder

		idx = "\t".join([str(numanc),str(numder),str(panelfreqA),str(panelfreqB)])
		countdict[idx] += 1

		print str(i)+" out of "+str(totalsites)
		i += 1

sortedkeys = countdict.keys()
sortedkeys.sort()
for entry in sortedkeys:
	print >>outfile, entry+"\t"+str(countdict[entry])

print "Total sites: "+str(totalsites)

