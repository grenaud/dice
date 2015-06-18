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
parser.add_option("--admixrateold", dest="admixrateold", help="Admixture rate old", default=0.0, type="float")
parser.add_option("--admixtimeold", dest="admixtimeold", help="Admixture time old", default=0.0255, type="float")
parser.add_option("-u", "--numtotalhum", dest="numtotalhum",help="Number of total humans sampled", default=100, type="int")
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
admixrateold = options.admixrateold
admixtimeold = options.admixtimeold

outfile = open(options.outfile,"w")


# Convert to units of 4N_0 generations
timesplitms = timesplit / 2.0

numsamphum = 6
numtotalhum = options.numtotalhum
numarch = 2

puredict = defaultdict(int)
i=1
try:
	while i < (numsim+1):

                # Create simulation file                                                                                                                                                                                                    
                infile_name = ""+destfolder+"/simul_"+str(i)+"_"+str(randint(0,10000))+".txt"
                if admixrate > 0:
                        commname = "msms "+str(numtotalhum+numarch)+" 1 -t "+str(theta)+" -I 2 "+str(numtotalhum)+" "+str(numarch)
                        #commname = commname +" -es "+str(admixtime)+" 1 "+str(1 - admixrate)+" -ej "+str(admixtime+0.00001)+" 3 2"+" -ej "+str(timesplitms)+" 1 2 | tail -n+7"
			commname = commname +" -es "+str(admixtime)+" 1 "+str(1 - admixrate)+" -ej "+str(admixtime+0.00001)+" 3 2"+" -ej "+str(timesplitms)+" 1 2 ";
			commname = commname +" -es "+str(admixtimeold)+" 2 "+str(1 - admixrateold)+ " ";  #splits archaic lineage into 2, creates lineage #4
			commname = commname +" -ej "+str(admixtimeold+0.00001)+" 4 1 "+" "; #join #4 to humans (#1)
			commname = commname +"| tail -n+7";
                else:
                        commname = "msms "+str(numtotalhum+numarch)+" 1 -t "+str(theta)+" -I 2 "+str(numtotalhum)+" "+str(numarch)+" -ej "+str(timesplitms)+" 1 2 | tail -n+7"
        #       print commname                                                                                                                                                                                                              
        #       i += 1                                                                                                                                                                                                                      
                commname = commname + " | sed 's/1/1,/g' | sed 's/0/0,/g' | sed 's/,$//'"
		
		commname = commname + " > "+infile_name
		print commname;
		v = subprocess.Popen(commname,shell=True)
		#v.communicate()
		v.wait()
		
		#try:
		if True:	
			print "simul number "+str(i)
			
			# Load and sum matrices
			sitemat = np.loadtxt(infile_name,delimiter=",",dtype="int")

			#print sitemat.transpose().shape
			#print sitemat.transpose()
			matnumind = sitemat.shape[0]
			matsize = sitemat.size
			print matsize / 102
			# Check if it's a vector
			if matnumind == matsize and matsize > 0:
				column = sitemat
				modfreq = float(np.sum(column[0:numtotalhum])) / float(numtotalhum)
				arch = np.sum(column[numtotalhum:(numtotalhum+numarch)])
				if (modfreq > 0 and modfreq < 1):
					puredict[(modfreq,arch)] += 1				
					# Advance if we've collected one site
					i += 1
				
			# Check if it's a matrix
			elif matsize > 0:
				itermat = sitemat.transpose()
				atleastone = False
				for column in itermat:
					modfreq = float(np.sum(column[0:numtotalhum])) / float(numtotalhum)
					arch = np.sum(column[numtotalhum:(numtotalhum+numarch)])
					# Require the site to be segregating
					if (modfreq > 0 and modfreq < 1):
						puredict[(modfreq,arch)] += 1
						atleastone = True
						
				# Advance if we've collected at least one site
				if atleastone == True:
					i += 1


		os.remove(infile_name)

except:
      print "Script interrupted"




countdict = defaultdict(int)
i = 1
totalsites = sum(puredict.values())
#totalsites = sum(puredict.keys())
print >>outfile, "Anc\tDer\tPanelFreq\tNum"
for key in puredict.keys():
	#panelsamp = int(key[1])
	panelfreq =  float(key[0])
	archgeno = int(key[1])
#	print "Key: "+str(key)
#	print "Number: "+str(puredict[key])
	coveragevec = np.random.poisson(meancoverage, puredict[key])

	for sitecoverage in coveragevec:
		if sitecoverage > 0:
			if(archgeno == 2):
				q = contrate*panelfreq*(1-errorrate) + contrate*(1-panelfreq)*errorrate + (1-contrate)*(1-errorrate)
			elif(archgeno == 1):
				q = contrate*panelfreq*(1-errorrate) + contrate*(1-panelfreq)*errorrate + (1-contrate)*(1-errorrate)/2.0 + (1-contrate)*errorrate/2.0
			elif(archgeno == 0):
				q = contrate*panelfreq*(1-errorrate)+contrate*(1-panelfreq)*errorrate + (1-contrate)*errorrate

       			numder = np.random.binomial(sitecoverage,q,1)[0]
	       		numanc = sitecoverage - numder

		       	idx = "\t".join([str(numanc),str(numder),str(panelfreq)])
			countdict[idx] += 1
			
			i += 1
		else:
			totalsites = totalsites - 1

sortedkeys = countdict.keys()
sortedkeys.sort()
for key in sortedkeys:
#	print key+"\t"+str(countdict[key])
	print >>outfile, key+"\t"+str(countdict[key])

print "Total sites with coverage > 0: "+str(totalsites)

